use crate::bloom_filter::BloomFilter;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

#[derive(Serialize, Deserialize)]
pub struct DeBruijnGraph {
    k: usize,
    bloom_filter: BloomFilter,
    critical_false_positives: HashSet<Vec<u8>>,
}

impl DeBruijnGraph {
    pub fn new(k: usize, reference_strings: Vec<String>) -> Self {
        let mut seen_edges: HashSet<Vec<u8>> = HashSet::new();
        // count edges
        for string in reference_strings.iter() {
            let mut last_node: Option<&[u8]> = None;
            for window in string.as_str().as_bytes().windows(k) {
                if let Some(last_node) = last_node {
                    let mut key: Vec<u8> = last_node.to_vec();
                    key.push(*window.last().unwrap());
                    seen_edges.insert(key);
                }
                last_node = Some(window);
            }
        }

        // build bloom filter of edges
        let mut bloom_filter = BloomFilter::with_fpr(0.05, seen_edges.len());
        for edge in seen_edges.iter() {
            bloom_filter.insert_key(edge.clone());
        }

        let mut critical_false_positives: HashSet<Vec<u8>> = HashSet::new();
        for string in reference_strings.iter() {
            for window in string.as_str().as_bytes().windows(k) {
                let mut key: Vec<u8> = window.to_vec();
                key.push(0);
                for &c in b"ACGT".iter() {
                    *key.last_mut().unwrap() = c;
                    if bloom_filter.query_key(&key) && !seen_edges.contains(&key) {
                        critical_false_positives.insert(key.clone());
                    }
                }
            }
        }

        Self {
            k,
            bloom_filter,
            critical_false_positives,
        }
    }

    pub fn query<S: AsRef<[u8]>>(&self, q: S) -> bool {
        for window in q.as_ref().windows(self.k + 1) {
            let key: Vec<u8> = window.to_vec();
            if !self.bloom_filter.query_key(&key) {
                return false;
            }
            if self.critical_false_positives.contains(&key) {
                return false;
            }
        }
        true
    }

    pub fn print_stats(&self) {
        println!("k: {}", self.k);
        println!(
            "bloom filter: {} bits, {} hashes",
            self.bloom_filter.n_bits, self.bloom_filter.n_hashes
        );
        println!(
            "critical false positives: {}",
            self.critical_false_positives.len()
        );
    }
}
