use serde::{Deserialize, Serialize};
use std::collections::HashMap;

pub trait MapLike: FromIterator<(Sequence, usize)> {
    fn get(&self, key: &Sequence) -> Option<usize>;
}

impl MapLike for HashMap<Sequence, usize> {
    fn get(&self, key: &Sequence) -> Option<usize> {
        self.get(key).copied()
    }
}

pub fn index_in_acgt(c: u8) -> usize {
    // Map A, C, G, T to 0, 1, 2, 3
    // The ilog2 compiles into a single `bsr` instruction, which is pretty neat
    debug_assert!(c == b'A' || c == b'C' || c == b'G' || c == b'T');
    (c - 0x3f).checked_ilog2().unwrap_or(1) as usize - 1
}

pub fn murmur_hash(key: &[u8], salt: u32) -> u32 {
    let mut acc = salt.wrapping_mul(0x5bd1e99).wrapping_add(0xc613fc15);
    acc ^= acc >> 15;
    for &c in key {
        acc ^= c as u32;
        acc = acc.wrapping_mul(0x5bd1e99);
        acc ^= acc >> 15;
    }
    acc
}

pub fn murmur_hash_n(key: &[u8], salt: u32, n: usize) -> usize {
    // If our hash function is good enough, this trick avoids an expensive modulo operation
    murmur_hash(key, salt) as usize * n >> 32
}

// 2-bit encoded sequence
#[derive(Clone, Debug, Hash, Eq, PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    bit_offset: u8,
    sequence: Vec<u8>,
}

impl Sequence {
    pub fn new() -> Self {
        Self {
            bit_offset: 0,
            sequence: vec![],
        }
    }

    pub fn len(&self) -> usize {
        self.sequence.len() * 4 - ((4 - self.bit_offset) % 4) as usize
    }

    pub fn extend_from_2bit_rep<S: AsRef<[u8]>>(&mut self, s: S) {
        let s = s.as_ref();
        let length = s.len();
        let remainder = self.bit_offset as usize;
        let mut pre = 0;

        // update length
        self.bit_offset = (self.bit_offset + length as u8) % 4;

        if remainder != 0 {
            pre = 4 - remainder;
            // complete current byte
            if length <= pre {
                *self.sequence.last_mut().unwrap() |= s
                    .into_iter()
                    .enumerate()
                    .fold(0, |acc, (i, &c)| acc | (c << (2 * (i + remainder))));
                return;
            }
            *self.sequence.last_mut().unwrap() |= s[..pre]
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * (i + remainder))));
        }

        let (chunks, remainder) = s[pre..].as_chunks::<4>();
        for &[a, b, c, d] in chunks {
            self.sequence.push(a + 4 * b + 16 * c + 64 * d);
        }
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * i))) as u8;
            self.sequence.push(b);
        }
    }

    pub fn check_genome<S: AsRef<[u8]>>(&self, s: S, pos: usize) -> bool {
        let s = s.as_ref();
        let length = s.len();
        if pos + length > self.len() {
            return false;
        }
        let (mut byte_pos, bit_pos) = (pos / 4, pos % 4);

        let mut pre = 0;
        if bit_pos != 0 {
            pre = 4 - bit_pos;
            let cur = self.sequence[byte_pos] >> (bit_pos * 2);
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
            if b != self.sequence[byte_pos] {
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
            return b == self.sequence[byte_pos] & mask;
        }

        true
    }

    pub fn from_genome<S: AsRef<[u8]>>(s: S) -> Self {
        let s = s.as_ref();
        let (chunks, remainder) = s.as_chunks::<4>();
        let mut sequence: Vec<u8> = chunks
            .into_iter()
            .map(|&[a, b, c, d]| {
                (index_in_acgt(a)
                    + 4 * index_in_acgt(b)
                    + 16 * index_in_acgt(c)
                    + 64 * index_in_acgt(d)) as u8
            })
            .collect();
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (index_in_acgt(c) << (2 * i)))
                as u8;
            sequence.push(b);
        }
        Self { bit_offset: remainder.len() as u8, sequence }
    }

    pub fn from_2bc<S: AsRef<[u8]>>(s: S) -> Self {
        let s = s.as_ref();
        let (chunks, remainder) = s.as_chunks::<4>();
        let mut sequence: Vec<u8> = chunks
            .into_iter()
            .map(|&[a, b, c, d]| a + 4 * b + 16 * c + 64 * d)
            .collect();
        if remainder.len() > 0 {
            let b = remainder
                .into_iter()
                .enumerate()
                .fold(0, |acc, (i, &c)| acc | (c << (2 * i))) as u8;
            sequence.push(b);
        }
        Self { bit_offset: remainder.len() as u8, sequence }
    }

    pub fn murmur_hash_n(&self, salt: u32, n: usize) -> usize {
        // we ignore length since we will only hash sequences of the same length
        murmur_hash_n(&self.sequence, salt, n)
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
}
