use std::collections::HashMap;

pub trait MapLike: FromIterator<(Vec<u8>, usize)> {
    fn get(&self, key: &Vec<u8>) -> Option<usize>;
}

impl MapLike for HashMap<Vec<u8>, usize> {
    fn get(&self, key: &Vec<u8>) -> Option<usize> {
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
}
