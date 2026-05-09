use crate::{
    bitvector::BitVector,
    mphf::MPHF,
    util::{MapLike, index_in_acgt},
};
use serde::{Deserialize, Serialize};
use std::{collections::{HashMap, HashSet}, time::Instant};
use log::trace;

#[derive(Debug, Serialize, Deserialize)]
pub struct PufferfishIndex<HM: MapLike> {
    k: usize,
    h: HM,
    useq: Vec<u8>,
    bv: BitVector,
    utab: Vec<u8>,
}

// HM can (and should) later be swapped out with a MPHF
pub type HashMapPufferfishIndex = PufferfishIndex<HashMap<Vec<u8>, usize>>;
pub type DefaultPufferfishIndex = PufferfishIndex<MPHF>;

impl<HM: MapLike> PufferfishIndex<HM> {
    pub fn new<S: AsRef<[u8]>>(k: usize, reference_strings: Vec<S>) -> Self {
        let n_colors = reference_strings.len();
        let color_bytes = n_colors.div_ceil(8);

        trace!("building pufferfish index (k = {})", k);
        let start = Instant::now();

        // Build De Bruijn graph
        // Hashmap stores forward/back edges, bitvec of colors, and if node is start/end
        let mut nodes: HashMap<&[u8], (u8, Vec<u8>, u8)> = HashMap::new();
        for (color, string) in reference_strings.iter().enumerate() {
            let (color_idx, color_bit) = (color / 8, 1u8 << (color % 8));

            let mut last_node: Option<&[u8]> = None;
            for window in string.as_ref().windows(k) {
                // insert node, add color
                let cur_entry = nodes
                    .entry(window)
                    .or_insert_with(|| (0, vec![0; color_bytes], 0));
                cur_entry.1[color_idx] |= color_bit;

                if let Some(last_node) = last_node {
                    // insert back edge
                    let bit_pattern_prev = 16 << index_in_acgt(last_node[0]);
                    cur_entry.0 |= bit_pattern_prev;

                    // insert forward edge
                    let next = *window.last().unwrap();
                    let bit_pattern_next = 1 << index_in_acgt(next);
                    nodes.get_mut(last_node).unwrap().0 |= bit_pattern_next;
                } else {
                    // mark node as start of sequence
                    cur_entry.2 |= 1;
                }
                last_node = Some(window);
            }

            if let Some(last_node) = last_node {
                // mark last node as end of sequence
                nodes.get_mut(last_node).unwrap().2 |= 2;
            }
        }

        trace!("build graph: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        // Find junctions (all places where a unipath starts)
        let mut junctions: HashSet<&[u8]> = HashSet::new();
        for (node, (edge_info, _, start_end_info)) in &nodes {
            if start_end_info & 1 == 1 {
                // sequence starts are junctions
                junctions.insert(node);
                continue;
            }

            let back_edges = edge_info >> 4;
            let one_back_edge = back_edges.is_power_of_two();

            if !one_back_edge {
                // 0 or multiple back edges is a junction
                junctions.insert(node);
                continue;
            }

            // look at previous node
            let mut prev_node = vec![b"ACGT"[back_edges.ilog2() as usize]];
            prev_node.extend(&node[..k - 1]);

            let (prev_edge_info, _, prev_start_end_info) = nodes[&prev_node[..]];

            if prev_start_end_info & 2 == 2 {
                // nodes following sequence ends are junctions
                junctions.insert(node);
                continue;
            }

            let prev_forward_edges = prev_edge_info & 0xf;
            let one_prev_forward_edge = prev_forward_edges.is_power_of_two();
            if !one_prev_forward_edge {
                // multiple forward edges on prev node is a junction
                junctions.insert(node);
            }
        }

        trace!("find junctions: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let mut useq: Vec<u8> = Vec::new();
        let mut pos: Vec<(Vec<u8>, usize)> = Vec::new();
        let mut bv: Vec<u64> = Vec::new();
        let mut utab: Vec<u8> = Vec::new();

        // Create a unipath at each junction
        for node in junctions.iter() {
            let mut unipath: Vec<u8> = node.to_vec();
            let mut key: Vec<u8> = node.to_vec();
            let mut useq_idx = useq.len();
            pos.push((key.clone(), useq_idx));
            useq_idx += 1;

            let (edge_info, colors, _) = &nodes[&key[..]];
            let mut forward_edges = edge_info & 0xf;

            // create entry in utab
            utab.extend(colors);

            loop {
                if !forward_edges.is_power_of_two() {
                    // multiple paths
                    break;
                }
                let next = b"ACGT"[forward_edges.ilog2() as usize];
                key.remove(0);
                key.push(next);

                if junctions.contains(&key[..]) {
                    // reached a junction
                    break;
                }

                let (edge_info, _, _) = &nodes[&key[..]];
                forward_edges = edge_info & 0xf;

                unipath.push(next);
                pos.push((key.clone(), useq_idx));
                useq_idx += 1;
            }
            useq.extend(unipath);

            // set bit in bv
            let bit_idx = useq.len() - 1;
            bv.resize_with(useq.len().div_ceil(64), Default::default);
            bv[bit_idx / 64] |= 1 << (bit_idx % 64);
        }

        trace!("find unipaths: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let h: HM = HM::from_iter(pos);

        trace!("create map: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let bv = BitVector::new(useq.len(), bv);

        trace!("create bv: {}ms", (Instant::now() - start).as_millis());

        Self {
            k,
            h,
            useq,
            bv,
            utab,
        }
    }

    pub fn query<S: AsRef<[u8]>>(&self, q: S) -> bool {
        // TODO: this should really query consecutive k-mers instead of using windows
        for window in q.as_ref().windows(self.k) {
            let pos = self.h.get(&window.to_vec());
            if let Some(pos) = pos {
                if self.useq[pos..pos + self.k] != *window {
                    // k-mer not in useq
                    return false;
                }
                let rank1 = self.bv.rank(pos);
                let rank2 = self.bv.rank(pos + self.k - 1);
                if rank1 != rank2 {
                    // crossed useq boundary
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_test() {
        let sequences = vec!["CTAAGAT", "CGATGCA", "TAAGAGG"];
        let index = DefaultPufferfishIndex::new(3, sequences);

        assert!(index.query("CTAAGAT"));
        assert!(index.query("CGATGCA"));
        assert!(index.query("TAAGAGG"));
    }
}
