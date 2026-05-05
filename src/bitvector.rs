use serde::{Deserialize, Serialize};
use std::{cmp::min, ops::Index};

#[derive(Debug, Serialize, Deserialize)]
struct IntVector {
    n_bits: u8, // must be <= 64
    buf: Vec<u64>,
}

impl IntVector {
    fn new(n_bits: u8, size: usize) -> Self {
        let buf_len = (size * n_bits as usize).div_ceil(64);
        let buf: Vec<u64> = vec![0; buf_len];
        Self { n_bits, buf }
    }

    fn get(&self, idx: usize) -> u64 {
        let mask = (1 << self.n_bits as usize) - 1;
        let bit_idx = idx * self.n_bits as usize;
        let (ele_idx, ele_offset) = (bit_idx / 64, bit_idx % 64);
        self.buf
            .get(ele_idx + 1)
            .map_or(0, |&x| x)
            .funnel_shr(self.buf[ele_idx], ele_offset as u32)
            & mask
    }

    // Accumulate ranks starting from `start_index`, using popcount from `iter`
    fn set_accumulated_ranks<I: IntoIterator<Item = u64>>(
        &mut self,
        start_index: usize,
        iter: I,
    ) -> u64 {
        let mut bit_idx = start_index * self.n_bits as usize;
        let mut acc = 0;
        for i in iter {
            let (ele_idx, ele_offset) = (bit_idx / 64, bit_idx % 64);

            self.buf[ele_idx] |= acc << ele_offset;
            if ele_offset + self.n_bits as usize > 64 {
                self.buf[ele_idx + 1] |= acc.unbounded_shr(64 - ele_offset as u32);
            }

            let next_offset = bit_idx + self.n_bits as usize;
            bit_idx = next_offset;
            acc += i;
        }
        acc
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BitVector {
    n_bits: usize,
    qwords: Vec<u64>,
    chunk_nbits: usize,
    n_subchunks: usize,
    subchunk_nbits: usize,
    cumulative_ranks: IntVector,
    cumulative_subranks: IntVector,
}

impl BitVector {
    pub fn new(n_bits: usize, qwords: Vec<u64>) -> Self {
        let n_bit_width = (n_bits + 1).bit_width() as usize;
        let chunk_nbits = n_bit_width * n_bit_width;
        let subchunk_nbits = n_bit_width.div_ceil(2);
        let n_chunks = n_bits.div_ceil(chunk_nbits);
        let n_subchunks = chunk_nbits.div_ceil(subchunk_nbits);

        let mut cumulative_ranks = IntVector::new(n_bit_width as u8, n_chunks);
        let mut cumulative_subranks = IntVector::new(chunk_nbits as u8, n_chunks * n_subchunks);

        let mut bit_idx = 0usize;

        cumulative_ranks.set_accumulated_ranks(
            0,
            (0..n_chunks).map(|i| {
                let bits_left = n_bits - i * chunk_nbits;
                let subchunks_left = bits_left.div_ceil(subchunk_nbits);
                cumulative_subranks.set_accumulated_ranks(
                    i * n_subchunks,
                    (0..min(subchunks_left, n_subchunks)).map(|j| {
                        // We want to read `subchunk_nbits` bits, unless we hit end of chunk
                        let nbits = min((j + 1) * subchunk_nbits, chunk_nbits) - j * subchunk_nbits;
                        let mask = (1 << nbits) - 1;

                        let (qword_idx, cur_shift) = (bit_idx / 64, bit_idx % 64);

                        let next = qwords.get(qword_idx + 1).map_or(0, |&x| x);
                        let cur = qwords[qword_idx];

                        let value = next.funnel_shr(cur, cur_shift as u32);

                        bit_idx += nbits;

                        // Count number of 1s set
                        (value & mask).count_ones() as u64
                    }),
                )
            }),
        );
        Self {
            n_bits,
            qwords,
            chunk_nbits,
            n_subchunks,
            subchunk_nbits,
            cumulative_ranks,
            cumulative_subranks,
        }
    }

    pub fn len(&self) -> usize {
        self.n_bits
    }

    pub fn access(&self, idx: usize) -> bool {
        let (qword_idx, shift) = (idx / 64, idx % 64);
        ((self.qwords[qword_idx] >> shift) & 1) != 0
    }

    pub fn rank(&self, idx: usize) -> Option<usize> {
        if idx > self.n_bits {
            return None;
        }
        let (chunk_idx, chunk_offset) = (idx / self.chunk_nbits, idx % self.chunk_nbits);
        let (subchunk_idx, nbits) = (
            chunk_offset / self.subchunk_nbits,
            chunk_offset % self.subchunk_nbits,
        );

        let chunk_rank = self.cumulative_ranks.get(chunk_idx);
        let subchunk_rank = self
            .cumulative_subranks
            .get(chunk_idx * self.n_subchunks + subchunk_idx);

        let mask = (1 << nbits) - 1;

        let subchunk_bit_idx = idx - nbits;
        let (qword_idx, shift) = (subchunk_bit_idx / 64, subchunk_bit_idx % 64);

        let next = self.qwords.get(qword_idx + 1).map_or(0, |&x| x);
        let cur = self.qwords[qword_idx];

        let final_rank = (next.funnel_shr(cur, shift as u32) & mask).count_ones() as u64;

        return Some((chunk_rank + subchunk_rank + final_rank) as usize);
    }

    pub fn select(&self, rank: usize) -> Option<usize> {
        let mut lb = 1;
        let mut ub = self.n_bits + 1;
        while lb < ub {
            let mid = lb + (ub - lb) / 2;

            if self.rank(mid).unwrap() > rank {
                ub = mid;
            } else {
                lb = mid + 1;
            }
        }
        if lb <= self.n_bits {
            Some(lb - 1)
        } else {
            None
        }
    }
}

impl Index<usize> for BitVector {
    type Output = bool;

    fn index(&self, index: usize) -> &bool {
        if self.access(index) { &true } else { &false }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intvec_test() {
        let mut intvec = IntVector::new(20, 1000);
        intvec.set_accumulated_ranks(0, 1..=1000);

        for i in 0..1000 {
            assert_eq!(intvec.get(i) as usize, i * (i + 1) / 2);
        }
    }

    #[test]
    fn bitvec_rank_test() {
        let bitvec = BitVector::new(13, vec![0x1456]);

        assert_eq!(bitvec.rank(2), Some(1));
        assert_eq!(bitvec.rank(8), Some(4));
    }

    #[test]
    fn bitvec_select_test() {
        let bitvec = BitVector::new(13, vec![0x1456]);

        assert_eq!(bitvec.select(0), Some(1));
        assert_eq!(bitvec.select(5), Some(12));
    }
}
