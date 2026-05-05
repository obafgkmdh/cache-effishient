pub fn index_in_acgt(c: u8) -> usize {
    // Map A, C, G, T to 0, 1, 2, 3
    // The ilog2 compiles into a single `bsr` instruction, which is pretty neat
    (c - 0x3f).checked_ilog2().unwrap_or(1) as usize - 1
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
