use crate::{bitvector::BitVector, util::index_in_acgt};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap};

pub trait HashMapLike: FromIterator<(Vec<u8>, usize)> {
    fn get(&self, index: &Vec<u8>) -> Option<&usize>;
}

#[derive(Serialize, Deserialize)]
pub struct PufferfishIndex<HM: HashMapLike> {
    k: usize,
    h: HM,
    useq: Vec<u8>,
    bv: BitVector,
}

#[derive(Serialize, Deserialize)]
pub struct HashMapWrapper {
    h: HashMap<Vec<u8>, usize>,
}

impl FromIterator<(Vec<u8>, usize)> for HashMapWrapper {
    fn from_iter<T: IntoIterator<Item=(Vec<u8>, usize)>>(iter: T) -> Self {
        Self {
            h: HashMap::from_iter(iter)
        }
    }
}

impl HashMapLike for HashMapWrapper {
    fn get(&self, index: &Vec<u8>) -> Option<&usize> {
        self.h.get(index)
    }
}

// HM can (and should) later be swapped out with a MPFH
pub type DefaultPufferfishIndex = PufferfishIndex<HashMapWrapper>;

impl<HM: HashMapLike> PufferfishIndex<HM> {
    pub fn new(k: usize, reference_strings: Vec<String>) -> Self {
        // Build De Bruijn graph
        let mut nodes: HashMap<Vec<u8>, u8> = HashMap::new();
        for string in reference_strings.iter() {
            let mut last_node: Option<&[u8]> = None;
            for window in string.as_str().as_bytes().windows(k) {
                if let Some(last_node) = last_node {
                    // insert forward edge
                    let mut key: Vec<u8> = last_node.to_vec();
                    let next = *window.last().unwrap();
                    *nodes.entry(key.clone()).or_default() |= 1 << index_in_acgt(next);
                    // insert back edge
                    key.push(next);
                    *nodes.entry(key[1..].to_vec()).or_default() |= 16 << index_in_acgt(key[0]);
                }
                last_node = Some(window);
            }
        }

        // Find unipaths
        // FIXME: this will not find some loops
        let mut useq: Vec<u8> = Vec::new();
        let mut pos: Vec<(Vec<u8>, usize)> = Vec::new();
        let mut bv: Vec<u64> = Vec::new();
        for (node, edge_info) in &nodes {
            let forward_edges = edge_info & 0xf;
            let back_edges = edge_info >> 4;
            if back_edges.is_power_of_two() {
                // not a junction
                continue;
            }
            // walk the graph
            for (i, &c) in b"ACGT".iter().enumerate() {
                if forward_edges >> i & 1 == 1 {
                    // start a new unipath
                    pos.push((node.clone(), useq.len()));

                    let mut unipath: Vec<u8> = node.clone();
                    unipath.push(c);
                    loop {
                        let key: Vec<u8> = unipath[unipath.len() - k..].to_vec();
                        if key == *node {
                            // we looped
                            break;
                        }

                        pos.push((key.clone(), useq.len() + unipath.len() - k));

                        let forward_edge_info = nodes[&key] & 0xf;
                        if !forward_edge_info.is_power_of_two() {
                            // reached junction
                            break;
                        }
                        unipath.push(b"ACGT"[forward_edge_info.ilog2() as usize]);
                    }
                    useq.extend(unipath);
                    let bit_idx = useq.len() - 1;
                    bv.resize_with(useq.len().div_ceil(64), Default::default);
                    bv[bit_idx / 64] |= 1 << (bit_idx % 64);
                }
            }
        }

        let h: HM = HM::from_iter(pos);
        let bv = BitVector::new(useq.len(), bv);
        Self { k, h, useq, bv }
    }

    pub fn query<S: AsRef<[u8]>>(&self, q: S) -> bool {
        for window in q.as_ref().windows(self.k) {
            let pos = self.h.get(&window.to_vec());
            if let Some(&pos) = pos {
                if self.useq[pos..pos + self.k] != *window {
                    return false;
                }
            } else {
                return false;
            }
        }
        true
    }

    pub fn print_stats(&self) {
        println!("k: {}", self.k);
        println!("useq size: {}", self.useq.len());
    }
}
