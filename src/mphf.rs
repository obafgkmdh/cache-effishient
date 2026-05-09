use crate::{
    bitvector::{IntVector, BitVector},
    util::{MapLike, murmur_hash_n},
};
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, default::Default};

#[derive(Debug, Serialize, Deserialize)]
pub struct MPHF {
    size: usize,
    bv: BitVector,
    values: IntVector,
}

impl FromIterator<(Vec<u8>, usize)> for MPHF {
    fn from_iter<T: IntoIterator<Item = (Vec<u8>, usize)>>(iter: T) -> Self {
        let mut kv: HashSet<(Vec<u8>, usize)> = iter.into_iter().collect();

        let max_value: usize = kv.iter().map(|(_, v)| *v).max().unwrap();
        let n_bits = max_value.bit_width();

        let size = kv.len();

        let mut values = IntVector::new(n_bits.try_into().unwrap(), size);

        let mut bv: Vec<u64> = Vec::with_capacity(size.div_ceil(64));

        let mut n_processed = 0;
        let mut bit_offset = 0;
        let mut salt = 0;
        while !kv.is_empty() {
            let n_keys = kv.len(); // number of keys to process this round

            // allocate collision bitvec
            let mut coll: Vec<u64> = vec![0; n_keys.div_ceil(64)];

            // allocate space in bv
            bv.resize_with((bit_offset + n_keys).div_ceil(64), Default::default);

            // find collisions
            for (key, value) in kv.iter() {
                let hash = murmur_hash_n(key, salt, n_keys);
                let (index, shift) = (hash / 64, hash % 64);
                let coll_bit = coll[index] >> shift & 1;
                if coll_bit == 1 {
                    // is a collision
                    continue;
                }

                let (bv_index, bv_shift) = ((bit_offset + hash) / 64, (bit_offset + hash) % 64);
                let bv_bit = bv[bv_index] >> bv_shift & 1;
                if bv_bit == 0 {
                    bv[bv_index] |= 1 << bv_shift;
                    values.set(n_processed + hash, *value as u64);
                } else {
                    // mark collision
                    bv[bv_index] &= !(1 << bv_shift);
                    coll[index] |= 1 << shift;
                }
            }

            // compact values
            let mut value_idx = n_processed;
            for i in 0..n_keys {
                let (bv_index, bv_shift) = (bit_offset / 64, bit_offset % 64);
                if bv[bv_index] >> bv_shift & 1 == 1 {
                    let v = values.get(n_processed + i);
                    values.set(value_idx, v);
                    value_idx += 1;
                }
                bit_offset += 1;
            }
            n_processed = value_idx;

            // retain only the keys that collided
            kv.retain(|(key, _)| {
                let hash = murmur_hash_n(key, salt, n_keys);
                let (index, shift) = (hash / 64, hash % 64);
                let coll_bit = coll[index] >> shift & 1;
                coll_bit == 1
            });

            salt += 1;
        }

        Self {
            size,
            bv: BitVector::new(bit_offset, bv),
            values,
        }
    }
}

impl MapLike for MPHF {
    fn get(&self, key: &Vec<u8>) -> Option<usize> {
        let mut salt = 0;
        let mut bit_offset = 0;
        let mut keys_remaining = self.size;
        loop {
            let hash = murmur_hash_n(key, salt, keys_remaining);
            let bv_bit = self.bv[bit_offset + hash];
            if bv_bit {
                let value_idx = self.bv.rank(bit_offset + hash).unwrap();
                return Some(self.values.get(value_idx) as usize);
            }

            bit_offset += keys_remaining;
            let weight = self.bv.rank(bit_offset).unwrap();
            keys_remaining = self.size - weight;
            salt += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mphf_test() {
        let kv: Vec<(Vec<u8>, usize)> = (0..1000)
            .map(|i| (u16::to_le_bytes(i as u16).to_vec(), i))
            .collect();

        let mphf = MPHF::from_iter(kv.iter().cloned());

        for (k, v) in kv.into_iter() {
            assert_eq!(mphf.get(&k), Some(v));
        }
    }
}
