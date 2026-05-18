use crate::{
    bitvector::{BitVector, IntVector},
    util::{MapLike, Sequence},
};
use serde::{Deserialize, Serialize};
use std::{cmp::max, collections::HashSet, default::Default};

#[derive(Debug, Serialize, Deserialize)]
pub struct MPHF {
    size: usize,
    bv: BitVector,
    values: IntVector,
}

impl MapLike for MPHF {
    fn from_hashset(h: HashSet<(Sequence, usize)>) -> Self {
        let mut max_value = 0;

        let mut kv: HashSet<(Sequence, usize)> = h;

        for (_k, v) in kv.iter() {
            max_value = max(max_value, *v);
        }

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
                let hash = key.murmur_hash_n(salt, n_keys);
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
                let hash = key.murmur_hash_n(salt, n_keys);
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

    fn get(&self, key: &Sequence) -> Option<usize> {
        let mut salt = 0;
        let mut bit_offset = 0;
        let mut keys_remaining = self.size;
        loop {
            let hash = key.murmur_hash_n(salt, keys_remaining);
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
        let kv: HashSet<(Sequence, usize)> = (0..1000)
            .map(|i| {
                (
                    Sequence::from_2bc(
                        (0..8)
                            .map(|j| (i >> (2 * j) & 3) as u8)
                            .collect::<Vec<u8>>(),
                    ),
                    i,
                )
            })
            .collect();

        let mphf = MPHF::from_hashset(kv.clone());

        for (k, v) in kv.into_iter() {
            assert_eq!(mphf.get(&k), Some(v));
        }
    }
}
