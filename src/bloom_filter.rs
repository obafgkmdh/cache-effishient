use crate::util::{murmur_hash_64, murmur_hash_step};
use log::debug;
use serde::{Deserialize, Serialize};

const LN_2: f64 = 0.6931471805599453094_f64;

#[derive(Serialize, Deserialize)]
pub struct BloomFilter {
    pub n_bits: usize,
    pub n_hashes_minus_one: u32,
    bv: Vec<u64>,
}

impl BloomFilter {
    pub fn with_fpr(fpr: f64, n_keys: usize) -> Self {
        let n_bits = (-(n_keys as f64) * fpr.log2() / LN_2).ceil() as usize;
        let n_hashes = ((n_bits as f64) / (n_keys as f64) * LN_2).ceil() as u32;

        debug!("bloom filter (fpr = {fpr}, n_keys = {n_keys}): {n_bits} bits, {n_hashes} hashes");

        Self {
            n_bits,
            n_hashes_minus_one: n_hashes - 1,
            bv: vec![0; (n_bits + 63) / 64],
        }
    }

    #[inline(always)]
    fn insert_using_hash(&mut self, h: u32) {
        let loc = h as usize * self.n_bits >> 32;
        // hot path; only check bounds in debug mode
        debug_assert!(loc < self.n_bits);
        let v = unsafe { self.bv.get_unchecked_mut(loc / 64) };
        *v |= 1 << (loc % 64);
    }

    #[inline(always)]
    fn query_using_hash(&self, h: u32) -> u64 {
        let loc = h as usize * self.n_bits >> 32;
        // hot path; only check bounds in debug mode
        debug_assert!(loc < self.n_bits);
        let v = unsafe { self.bv.get_unchecked(loc / 64) };
        *v >> (loc % 64) & 1
    }

    pub fn insert_key<S: AsRef<[u8]>>(&mut self, key: S) {
        let h = murmur_hash_64(key.as_ref());
        let (h_high, mut h_low) = ((h >> 32) as u32, h as u32);

        // check first two hashes
        self.insert_using_hash(h_high);
        self.insert_using_hash(h_low);

        // check remaining hashes
        for i in 1..self.n_hashes_minus_one {
            h_low = murmur_hash_step(h_low, i.wrapping_mul(h_high));
            self.insert_using_hash(h_low);
        }
    }

    pub fn query_key<S: AsRef<[u8]>>(&self, key: S) -> bool {
        let h = murmur_hash_64(key.as_ref());
        let (h_high, mut h_low) = ((h >> 32) as u32, h as u32);

        // check first two hashes
        if self.query_using_hash(h_high) == 0 || self.query_using_hash(h_low) == 0 {
            return false;
        }

        // check remaining hashes
        (1..self.n_hashes_minus_one).all(|i| {
            h_low = murmur_hash_step(h_low, i.wrapping_mul(h_high));
            self.query_using_hash(h_low) != 0
        })
    }

    pub fn num_hashes(&self) -> u32 {
        self.n_hashes_minus_one + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intvec_test() {
        let mut bitvec = BloomFilter::with_fpr(0.05, 10);
        bitvec.insert_key("abc");
        bitvec.insert_key("def");

        assert_eq!(bitvec.query_key("abc"), true);
        assert_eq!(bitvec.query_key("def"), true);
    }
}
