use crate::util::murmur_hash_n;
use serde::{Deserialize, Serialize};

const LN_2: f64 = 0.6931471805599453094_f64;

#[derive(Serialize, Deserialize)]
pub struct BloomFilter {
    pub n_bits: usize,
    pub n_hashes: u32,
    bv: Vec<u8>,
}

impl BloomFilter {
    pub fn with_fpr(fpr: f64, n_keys: usize) -> Self {
        let n_bits = (-(n_keys as f64) * fpr.log2() / LN_2).ceil() as usize;
        let n_hashes = ((n_bits as f64) / (n_keys as f64) * LN_2).ceil() as u32;

        Self {
            n_bits,
            n_hashes,
            bv: vec![0; (n_bits + 7) / 8],
        }
    }

    fn get_position<S: AsRef<[u8]>>(&self, key: S, salt: u32) -> usize {
        murmur_hash_n(key.as_ref(), salt, self.n_bits)
    }

    pub fn insert_key<S: AsRef<[u8]>>(&mut self, key: S) {
        for i in 1..=self.n_hashes {
            let loc = self.get_position(&key, i);
            self.bv[loc / 8] |= 1 << (loc % 8);
        }
    }

    pub fn query_key<S: AsRef<[u8]>>(&self, key: S) -> bool {
        for i in 1..=self.n_hashes {
            let loc = self.get_position(&key, i);
            if ((self.bv[loc / 8] >> (loc % 8)) & 1) == 0 {
                return false;
            }
        }
        true
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
