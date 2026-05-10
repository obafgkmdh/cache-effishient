use crate::{
    bitvector::BitVector,
    mphf::MPHF,
    util::{MapLike, Sequence, index_in_acgt},
};
use log::trace;
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    time::Instant,
};

#[derive(Debug, Serialize, Deserialize)]
pub struct PufferfishIndex<HM: MapLike> {
    k: usize,
    h: HM,
    useq: Sequence,
    bv: BitVector,
    utab: Vec<u8>,
}

// HM can (and should) later be swapped out with a MPHF
pub type HashMapPufferfishIndex = PufferfishIndex<HashMap<Sequence, usize>>;
pub type DefaultPufferfishIndex = PufferfishIndex<MPHF>;

impl<HM: MapLike> PufferfishIndex<HM> {
    pub fn new<S: AsRef<[u8]>>(k: usize, reference_strings: Vec<S>) -> Self {
        let n_colors = reference_strings.len();
        let color_bytes = n_colors.div_ceil(8);

        trace!("building pufferfish index (k = {})", k);
        let start = Instant::now();

        let reference_strings: Vec<Vec<u8>> = reference_strings
            .into_iter()
            .map(|s| {
                s.as_ref()
                    .into_iter()
                    .map(|&c| index_in_acgt(c) as u8)
                    .collect()
            })
            .collect();

        trace!(
            "to 2-bit representation: {}ms",
            (Instant::now() - start).as_millis()
        );

        let start = Instant::now();

        // Build De Bruijn graph
        let mut starts: HashSet<&[u8]> = HashSet::new();
        let mut ends: HashSet<&[u8]> = HashSet::new();

        // Hashmap stores forward/back edges
        let mut nodes: HashMap<&[u8], u8> = HashMap::new();
        for string in reference_strings.iter() {
            let mut last_node: Option<&[u8]> = None;
            for window in string.windows(k) {
                // insert node
                let cur_entry = nodes.entry(window).or_insert(0);

                if let Some(last_node) = last_node {
                    // insert back edge
                    let bit_pattern_prev = 16 << last_node[0];
                    *cur_entry |= bit_pattern_prev;

                    // insert forward edge
                    let next = *window.last().unwrap();
                    let bit_pattern_next = 1 << next;
                    *nodes.get_mut(last_node).unwrap() |= bit_pattern_next;
                } else {
                    // mark node as start of sequence
                    starts.insert(window);
                }
                last_node = Some(window);
            }

            if let Some(last_node) = last_node {
                // mark last node as end of sequence
                ends.insert(last_node);
            }
        }

        trace!("build graph: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        // Find junctions (all places where a unipath starts)
        // Hashmap will store color information
        let mut junctions: HashMap<&[u8], Vec<u8>> =
            HashMap::from_iter(starts.iter().map(|&k| (k, vec![0; color_bytes])));
        for (node, edge_info) in &nodes {
            if starts.contains(node) {
                continue;
            }

            let back_edges = edge_info >> 4;
            let one_back_edge = back_edges.is_power_of_two();

            if !one_back_edge {
                // multiple back edges is a junction
                junctions.insert(node, vec![0; color_bytes]);
                continue;
            }

            // look at previous node
            let mut prev_node = vec![back_edges.ilog2() as u8];
            prev_node.extend(&node[..k - 1]);

            if ends.contains(&prev_node[..]) {
                // nodes following sequence ends are junctions
                junctions.insert(node, vec![0; color_bytes]);
                continue;
            }

            let prev_edge_info = nodes[&prev_node[..]];

            let prev_forward_edges = prev_edge_info & 0xf;
            let one_prev_forward_edge = prev_forward_edges.is_power_of_two();
            if !one_prev_forward_edge {
                // multiple forward edges on prev node is a junction
                junctions.insert(node, vec![0; color_bytes]);
            }
        }

        trace!("find junctions: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        for (color, string) in reference_strings.iter().enumerate() {
            let (color_idx, color_bit) = (color / 8, 1u8 << (color % 8));

            for window in string.windows(k) {
                // add color
                junctions.entry(window).and_modify(|info| {
                    info[color_idx] |= color_bit;
                });
            }
        }

        trace!(
            "mark junction colors: {}ms",
            (Instant::now() - start).as_millis()
        );

        let start = Instant::now();

        let mut useq: Sequence = Sequence::new();
        let mut pos: Vec<(Sequence, usize)> = Vec::new();
        let mut bv: Vec<u64> = Vec::new();
        let mut utab: Vec<u8> = Vec::new();

        // Create a unipath at each junction
        for (&node, colors) in junctions.iter() {
            let mut unipath: Vec<u8> = node.to_vec();
            let mut key: Vec<u8> = node.to_vec();
            let mut useq_idx = useq.len();
            pos.push((Sequence::from_2bc(&key), useq_idx));
            useq_idx += 1;

            let edge_info = &nodes[&key[..]];
            let mut forward_edges = edge_info & 0xf;

            // create entry in utab
            utab.extend(colors);

            loop {
                if !forward_edges.is_power_of_two() {
                    // multiple paths
                    break;
                }
                let next = forward_edges.ilog2() as u8;
                key.remove(0);
                key.push(next);

                if junctions.contains_key(&key[..]) {
                    // reached a junction
                    break;
                }

                let edge_info = &nodes[&key[..]];
                forward_edges = edge_info & 0xf;

                unipath.push(next);
                pos.push((Sequence::from_2bc(&key), useq_idx));
                useq_idx += 1;
            }
            useq.extend_from_2bit_rep(unipath);

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
            let pos = self.h.get(&Sequence::from_genome(window));
            if let Some(pos) = pos {
                if !self.useq.check_genome(window, pos) {
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
    fn basic_hashmap_test() {
        let sequences = vec!["CTAAGAT", "CGATGCA", "TAAGAGG"];
        let index = HashMapPufferfishIndex::new(3, sequences);

        assert!(index.query("CTAAGAT"));
        assert!(index.query("CGATGCA"));
        assert!(index.query("TAAGAGG"));
    }

    #[test]
    fn basic_mphf_test() {
        let sequences = vec!["CTAAGAT", "CGATGCA", "TAAGAGG"];
        let index = DefaultPufferfishIndex::new(3, sequences);

        assert!(index.query("CTAAGAT"));
        assert!(index.query("CGATGCA"));
        assert!(index.query("TAAGAGG"));
    }
}
