use serde::{Deserialize, Serialize};
use std::{
    cmp::max,
    collections::HashMap,
    hash::{DefaultHasher, Hasher},
};

pub trait MapLike: FromIterator<(Sequence, usize)> {
    fn get(&self, key: &Sequence) -> Option<usize>;
}

impl MapLike for HashMap<Sequence, usize> {
    fn get(&self, key: &Sequence) -> Option<usize> {
        self.get(key).copied()
    }
}

#[inline(always)]
pub fn index_in_acgt(c: u8) -> usize {
    // Map A, C, G, T to 0, 1, 2, 3
    // The ilog2 compiles into a single `bsr` instruction, which is pretty neat
    debug_assert!(c == b'A' || c == b'C' || c == b'G' || c == b'T');
    (c - 0x3f).checked_ilog2().unwrap_or(1) as usize - 1
}

#[inline(always)]
pub fn murmur_hash_step(state: u32, c: u32) -> u32 {
    let mut acc = state;
    acc ^= c;
    acc = acc.wrapping_mul(0x5bd1e99);
    acc ^= acc >> 15;
    acc
}

pub fn murmur_hash(key: &[u8], salt: u32) -> u32 {
    let mut acc = salt.wrapping_mul(0x5bd1e99).wrapping_add(0xc613fc15);
    acc ^= acc >> 15;
    key.iter()
        .fold(acc, |acc, &c| murmur_hash_step(acc, c as u32))
}

pub fn murmur_hash_64(key: &[u8]) -> u64 {
    let mut acc = 0x749e3e6989df617u64;
    let (chunks, remainder) = key.as_chunks::<8>();
    for &c in chunks {
        acc ^= u64::from_le_bytes(c);
        acc = acc.wrapping_mul(0x5bd1e9955bd1e995u64);
        acc ^= acc >> 47;
    }
    for (i, &c) in remainder.iter().enumerate() {
        acc ^= (c as u64) << (8 * i);
    }
    acc = acc.wrapping_mul(0x5bd1e9955bd1e995u64);
    acc ^= acc >> 47;
    acc
}

// 2-bit encoded sequence, packed to 4 chars per byte
// First byte of vector encodes length mod 4
#[derive(Clone, Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub struct Sequence(Vec<u8>);

impl Sequence {
    pub fn new() -> Self {
        Self(vec![0])
    }

    pub fn with_capacity(length: usize) -> Self {
        let mut inner = Vec::with_capacity(1 + length.div_ceil(8));
        inner.push(0);
        Self(inner)
    }

    #[inline(always)]
    fn get_bit_offset(&self) -> u8 {
        // first element should always exist
        debug_assert!(self.0.len() > 0);
        unsafe { *self.0.get_unchecked(0) }
    }

    pub fn len(&self) -> usize {
        (self.0.len() - 1) * 4 - ((4 - self.get_bit_offset()) % 4) as usize
    }

    // Extend sequence with a 2-bit-encoded byte array
    pub fn extend_from_2bit_rep<S: AsRef<[u8]>>(&mut self, s: S) {
        let s = s.as_ref();
        let length = s.len();
        let remainder = self.get_bit_offset() as usize;
        let mut pre = 0;

        // update length
        self.0[0] = (self.get_bit_offset() + length as u8) % 4;

        if remainder != 0 {
            pre = 4 - remainder;
            // complete current byte

            debug_assert!(self.0.len() > 1);
            let last_byte = unsafe { self.0.last_mut().unwrap_unchecked() };

            if length <= pre {
                *last_byte |= s
                    .into_iter()
                    .enumerate()
                    .fold(0, |acc, (i, &c)| acc | (c << (2 * (i + remainder))));
                return;
            }
            *last_byte |= s[..pre]
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * (i + remainder))));
        }

        let (chunks, remainder) = s[pre..].as_chunks::<4>();
        for &[a, b, c, d] in chunks {
            self.0.push(a + 4 * b + 16 * c + 64 * d);
        }
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * i))) as u8;
            self.0.push(b);
        }
    }

    // Check if ACGT-encoded sequence exists at the given position
    pub fn check_genome<S: AsRef<[u8]>>(&self, s: S, pos: usize) -> bool {
        let s = s.as_ref();
        let length = s.len();
        if pos + length > self.len() {
            return false;
        }
        let (mut byte_pos, bit_pos) = (pos / 4 + 1, pos % 4);

        let mut pre = 0;
        if bit_pos != 0 {
            pre = 4 - bit_pos;
            let cur = self.0[byte_pos] >> (bit_pos * 2);
            if length <= pre {
                let mask = (1 << (2 * length)) - 1;
                let b = s
                    .into_iter()
                    .enumerate()
                    .fold(0, |acc, (i, &c)| acc | (index_in_acgt(c) << (2 * i)) as u8);

                return b == (cur & mask);
            }

            let b = s[..pre]
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (index_in_acgt(c) << (2 * i)) as u8);
            if b != cur {
                return false;
            }

            byte_pos += 1;
        }

        let (chunks, remainder) = s[pre..].as_chunks::<4>();
        for &[a, b, c, d] in chunks {
            let b = (index_in_acgt(a)
                + 4 * index_in_acgt(b)
                + 16 * index_in_acgt(c)
                + 64 * index_in_acgt(d)) as u8;
            if b != self.0[byte_pos] {
                return false;
            }

            byte_pos += 1;
        }

        if remainder.len() > 0 {
            let mask = (1 << (2 * remainder.len())) - 1;
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (index_in_acgt(c) << (2 * i)))
                as u8;
            return b == self.0[byte_pos] & mask;
        }

        true
    }

    // Create sequence from ACGT-encoded byte array
    pub fn from_genome<S: AsRef<[u8]>>(s: S) -> Self {
        let s = s.as_ref();
        let (chunks, remainder) = s.as_chunks::<4>();
        let mut sequence: Vec<u8> = std::iter::once(remainder.len() as u8)
            .chain(chunks.into_iter().map(|&[a, b, c, d]| {
                (index_in_acgt(a)
                    + 4 * index_in_acgt(b)
                    + 16 * index_in_acgt(c)
                    + 64 * index_in_acgt(d)) as u8
            }))
            .collect();
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (index_in_acgt(c) << (2 * i)))
                as u8;
            sequence.push(b);
        }
        Self(sequence)
    }

    // Create sequence from 2-bit-encoded byte array
    pub fn from_2bc<S: AsRef<[u8]>>(s: S) -> Self {
        let s = s.as_ref();
        let (chunks, remainder) = s.as_chunks::<4>();
        let mut sequence: Vec<u8> = std::iter::once(remainder.len() as u8)
            .chain(
                chunks
                    .into_iter()
                    .map(|&[a, b, c, d]| a + 4 * b + 16 * c + 64 * d),
            )
            .collect();
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * i))) as u8;
            sequence.push(b);
        }
        Self(sequence)
    }

    #[inline(always)]
    pub fn murmur_hash_n(&self, salt: u32, n: usize) -> usize {
        // If our hash function is good enough, this trick avoids an expensive modulo operation
        // We ignore sequence length since we will only hash sequences of the same length
        murmur_hash(&self.0[1..], salt) as usize * n >> 32
    }
}

// HyperLogLog sketch: https://en.wikipedia.org/wiki/HyperLogLog
pub struct HyperLogLog {
    log_m: u8,
    mask: u64,
    alpha: f64,
    counters: Vec<u8>,
}

impl HyperLogLog {
    // Create new HyperLogLog sketch with 2^log_m counters
    pub fn new(log_m: u8) -> Self {
        let mask: u64 = (1 << log_m) - 1;
        let alpha = match log_m {
            4 => 0.673,
            5 => 0.697,
            6 => 0.709,
            log_m => 0.7213 / (1f64 + 1.079 / ((1u64 << log_m) as f64)),
        };
        Self {
            log_m,
            mask,
            alpha,
            counters: vec![0; 1 << log_m],
        }
    }

    // Insert an item
    pub fn insert<S: AsRef<[u8]>>(&mut self, s: S) {
        // We use DefaultHasher because hash quality matters more than performance here
        let mut hasher = DefaultHasher::new();
        hasher.write(s.as_ref());
        let h: u64 = hasher.finish();
        let j = (h & self.mask) as usize;

        let remain = h & !self.mask;
        let rho = remain.trailing_zeros() as u8 - self.log_m + 1;

        self.counters[j] = max(self.counters[j], rho);
    }

    // Estimate number of distinct items inserted
    pub fn count(&self) -> usize {
        let n_empty = self.counters.iter().filter(|&&c| c == 0).count();
        if n_empty > 0 {
            // small range correction
            let m = (1u64 << self.log_m) as f64;
            return (m * (m / n_empty as f64).ln()) as usize;
        }
        let counter_max = 65 - self.log_m;
        let denominator: u64 = self
            .counters
            .iter()
            .map(|&c| 1u64 << (counter_max - c))
            .sum();
        let harmonic_mean = (1u64 << counter_max) as f64 / denominator as f64;
        let estimate = self.alpha * (1u64 << (2 * self.log_m)) as f64 * harmonic_mean;
        estimate.round() as usize
    }

    // Standard error of count estimate
    pub fn standard_error(&self) -> f64 {
        let m = 1u64 << self.log_m;
        1.04f64 / (m as f64).sqrt()
    }

    // Compute upper bound on count using stddev
    // e.g., stddevs = 2 is greater than actual count about 97.7% of the time
    pub fn upper_bound(&self, stddevs: f64) -> usize {
        (self.count() as f64 * (1.0f64 + self.standard_error() * stddevs)) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_in_acgt() {
        assert_eq!(index_in_acgt(b'A'), 0);
        assert_eq!(index_in_acgt(b'C'), 1);
        assert_eq!(index_in_acgt(b'G'), 2);
        assert_eq!(index_in_acgt(b'T'), 3);
    }

    #[test]
    fn test_sequence() {
        let s = b"TAACGTCTTCAAGGCG";
        let seq = Sequence::from_genome(s);
        for w in 1..s.len() {
            for (i, k) in s.windows(w).enumerate() {
                assert!(seq.check_genome(k, i), "{} {:?} {:?}", i, k, seq);
            }
        }
    }

    #[test]
    fn test_hyperloglog() {
        let mut sketch = HyperLogLog::new(11);
        let mut check_counts: Vec<usize> = vec![100, 500, 1000, 5_000, 10_000, 100_000, 1_000_000];

        if !cfg!(debug_assertions) {
            // we can go higher in release mode
            check_counts.push(10_000_000);
            check_counts.push(100_000_000);
        }

        for i in 1usize..=*check_counts.iter().max().unwrap() {
            sketch.insert(i.to_le_bytes());

            if check_counts.contains(&i) {
                let count = sketch.count();
                let error = count as i64 - i as i64;

                // error should probably be within 3.5 standard deviations
                let stddev = sketch.standard_error() * i as f64;
                let z = error as f64 / stddev;
                eprintln!("estimated {}, true count was {} (z = {})", count, i, z);
                assert!(z < 3.5_f64);
            }
        }
    }
}
