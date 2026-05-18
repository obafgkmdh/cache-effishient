use rayon::prelude::*;

use crate::{
    bitvector::BitVector,
    mphf::MPHF,
    util::{HyperLogLog, MapLike, Sequence, index_in_acgt, murmur_hash_64},
};
use log::{debug, trace};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    time::Instant,
};
use xorf::{BinaryFuse8, Filter};

#[derive(Debug, Serialize, Deserialize)]
pub struct PufferfishIndex<HM: MapLike> {
    k: usize,
    n_colors: usize,
    h: HM,
    useq: Sequence,
    bv: BitVector,
    utab: Vec<u8>,
    identifiers: Vec<String>,
}

pub type HashMapPufferfishIndex = PufferfishIndex<HashMap<Sequence, usize>>;
pub type DefaultPufferfishIndex = PufferfishIndex<MPHF>;

impl<HM: MapLike> PufferfishIndex<HM> {
    pub fn new<S: Into<String>, T: AsRef<[u8]>>(k: usize, reference_strings: Vec<(S, T)>) -> Self {
        let n_colors = reference_strings.len();
        let color_bytes = n_colors.div_ceil(8);

        // Hashmap will store color information
        let mut junctions: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();

        trace!("building pufferfish index (k = {})", k);
        let start = Instant::now();

        // convert to 2-bit representation
        let mut identifiers: Vec<String> = Vec::with_capacity(n_colors);
        let reference_strings: Vec<Vec<u8>> = reference_strings
            .into_iter()
            .map(|(identifier, sequence)| {
                identifiers.push(identifier.into());
                sequence
                    .as_ref()
                    .into_iter()
                    .map(|&c| index_in_acgt(c) as u8)
                    .collect()
            })
            .collect();

        trace!("to 2-bit repr: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let mut unique_nodes_sketch = HyperLogLog::new(11);
        let mut unique_edges_sketch = HyperLogLog::new(11);

        for string in reference_strings.iter() {
            unique_nodes_sketch.insert(&string[..k]);
            for window in string.windows(k + 1) {
                // insert edge
                unique_edges_sketch.insert(&window);
                unique_nodes_sketch.insert(&window[1..]);
            }
        }

        let unique_edges_bound = unique_edges_sketch.upper_bound(2.0);
        let unique_nodes_bound = unique_nodes_sketch.upper_bound(2.0);

        trace!(
            "approx. count unique edges: {}ms",
            (Instant::now() - start).as_millis()
        );
        debug!(
            "approx. bounds on unique edges: {}, unique nodes: {}",
            unique_edges_bound, unique_nodes_bound
        );

        let start = Instant::now();

        // Build De Bruijn graph
        let mut ends: HashSet<&[u8]> = HashSet::with_capacity(reference_strings.len());

        let mut unique_hashes: HashSet<u64> = HashSet::with_capacity(unique_edges_bound);
        for string in reference_strings.iter() {
            // mark starts as junctions
            junctions.insert(string[..k].to_vec(), vec![0; color_bytes]);

            // mark ends
            ends.insert(&string[string.len() - k..]);

            for window in string.windows(k + 1) {
                // insert edge
                unique_hashes.insert(murmur_hash_64(window));
            }
        }

        let unique_hashes: Vec<u64> = unique_hashes.into_iter().collect();

        let edges_filter = BinaryFuse8::try_from(&unique_hashes).unwrap();

        trace!(
            "build xor filter: {}ms",
            (Instant::now() - start).as_millis()
        );

        // Find junctions (all places where a unipath starts)

        let start = Instant::now();

        let (
            // false positive forward edges
            mut critical_false_positives,
            mut maybe_branch,
            mut maybe_multiple_back_edges,
        ): (HashSet<Vec<u8>>, HashSet<&[u8]>, HashSet<Vec<u8>>) = reference_strings
            .par_iter()
            .map(|string| {
                let mut critical_false_positives: HashSet<Vec<u8>> = HashSet::new();
                let mut maybe_branch: HashSet<&[u8]> = HashSet::new();
                let mut maybe_multiple_back_edges: HashSet<Vec<u8>> = HashSet::new();

                // we don't check backwards edges from start, since we already know starts are
                // junctions
                for window in string.windows(k + 1) {
                    // check other forward edges
                    let mut key: Vec<u8> = window.to_vec();
                    for c in [1, 2, 3] {
                        key[k] = window[k] ^ c;
                        if edges_filter.contains(&murmur_hash_64(&key)) {
                            critical_false_positives.insert(key.clone());
                            maybe_branch.insert(&window[..k]);
                        }
                    }

                    // reset key
                    key[k] = window[k];

                    // check other backwards edges
                    for c in [1, 2, 3] {
                        key[0] = window[0] ^ c;
                        if edges_filter.contains(&murmur_hash_64(&key)) {
                            maybe_multiple_back_edges.insert(key.clone());
                        }
                    }
                }

                // check forward edges from end
                let mut key: Vec<u8> = string[string.len() - k..].to_vec();
                key.push(0);

                for c in [0, 1, 2, 3] {
                    key[k] = c;
                    if edges_filter.contains(&murmur_hash_64(&key)) {
                        critical_false_positives.insert(key.clone());
                    }
                }

                (
                    critical_false_positives,
                    maybe_branch,
                    maybe_multiple_back_edges,
                )
            })
            .reduce(
                || (HashSet::new(), HashSet::new(), HashSet::new()),
                |(mut a1, mut a2, mut a3), (b1, b2, b3)| {
                    a1.extend(b1);
                    a2.extend(b2);
                    a3.extend(b3);
                    (a1, a2, a3)
                },
            );

        trace!(
            "find cfp candidates: {}ms",
            (Instant::now() - start).as_millis()
        );

        debug!(
            "{} cfp candids, {} branch candids, {} multiple back edges candids",
            critical_false_positives.len(),
            maybe_branch.len(),
            maybe_multiple_back_edges.len()
        );

        let start = Instant::now();

        // remove critical true positives
        for string in reference_strings.iter() {
            for window in string.windows(k + 1) {
                critical_false_positives.remove(window);
                if maybe_multiple_back_edges.remove(window) {
                    junctions
                        .entry(window[1..].to_vec())
                        .or_insert_with(|| vec![0; color_bytes]);
                }
            }
        }

        trace!(
            "identify true positives: {}ms",
            (Instant::now() - start).as_millis()
        );
        debug!(
            "{} cfps, {} junctions so far",
            critical_false_positives.len(),
            junctions.len()
        );

        let is_real_edge = |key: &[u8]| -> bool {
            edges_filter.contains(&murmur_hash_64(key)) && !critical_false_positives.contains(key)
        };

        let start = Instant::now();

        maybe_branch.retain(|node| {
            let mut key: Vec<u8> = node.to_vec();
            key.push(0);

            let mut count = 0;
            for c in [0, 1, 2, 3] {
                key[k] = c;
                if is_real_edge(&key) {
                    count += 1;
                    if count > 1 {
                        return true;
                    }
                }
            }

            false
        });

        let branch_or_end: HashSet<&[u8]> = maybe_branch.union(&ends).cloned().collect();

        trace!(
            "identify branches: {}ms",
            (Instant::now() - start).as_millis()
        );
        debug!("{} branches/ends", branch_or_end.len());

        let start = Instant::now();

        // Find the rest of the junctions
        for node in branch_or_end.iter() {
            let mut key: Vec<u8> = node.to_vec();
            key.push(0);
            for c in [0, 1, 2, 3] {
                key[k] = c;
                if is_real_edge(&key) {
                    // nodes following branch/end are junctions
                    junctions.insert(key[1..].to_vec(), vec![0; color_bytes]);
                }
            }
        }

        trace!(
            "find remaining junctions: {}ms",
            (Instant::now() - start).as_millis()
        );
        debug!("number of junctions: {}", junctions.len());

        let start = Instant::now();

        for (color, string) in reference_strings.iter().enumerate() {
            let (color_idx, color_bit) = (color / 8, 1u8 << (color % 8));

            for window in string.windows(k) {
                // add color
                junctions.entry(window.to_vec()).and_modify(|info| {
                    info[color_idx] |= color_bit;
                });
            }
        }

        trace!(
            "mark junction colors: {}ms",
            (Instant::now() - start).as_millis()
        );

        let start = Instant::now();

        let estimated_useq_size = junctions.len() * (k - 1) + unique_nodes_bound;

        debug!("estimated useq size bound: {}", estimated_useq_size);

        let mut useq: Sequence = Sequence::with_capacity(estimated_useq_size);
        let mut pos: Vec<(Sequence, usize)> = Vec::with_capacity(unique_nodes_bound);
        let mut bv: Vec<u64> = Vec::with_capacity(estimated_useq_size.div_ceil(64));
        let mut utab: Vec<u8> = Vec::with_capacity(junctions.len() * color_bytes);

        // Create a unipath at each junction
        for (node, colors) in junctions.iter() {
            let mut unipath: Vec<u8> = node.clone();
            let mut key: Vec<u8> = node.clone();
            let mut useq_idx = useq.len();
            pos.push((Sequence::from_2bc(&key), useq_idx));
            useq_idx += 1;

            // create entry in utab
            utab.extend(colors);

            loop {
                if branch_or_end.contains(&key[..]) {
                    // no or multiple paths
                    break;
                }
                key.push(0);
                for c in [0, 1, 2, 3] {
                    key[k] = c;
                    if is_real_edge(&key) {
                        // found forward edge
                        break;
                    }
                }
                key.remove(0);

                if junctions.contains_key(&key) {
                    // reached a junction
                    break;
                }

                unipath.push(key[k - 1]);
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
        debug!("number of unipaths: {}", utab.len() / color_bytes);
        debug!("useq size: {}", useq.len());

        let start = Instant::now();

        let h: HM = HM::from_iter(pos);

        trace!("create map: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let bv = BitVector::new(useq.len(), bv);

        trace!("create bv: {}ms", (Instant::now() - start).as_millis());

        Self {
            k,
            n_colors,
            h,
            useq,
            bv,
            utab,
            identifiers,
        }
    }

    // query a single k-mer, return which unitig it's in
    pub fn query_kmer<S: AsRef<[u8]>>(
        &self,
        kmer: S,
        hint: Option<(usize, usize)>,
    ) -> Option<(usize, usize)> {
        let kmer = kmer.as_ref();
        debug_assert!(kmer.len() == self.k);

        if let Some((prev_pos, prev_rank)) = hint {
            let pos_start = prev_pos + self.k;
            // does this k-mer appear right after the previous one?
            if self.useq.check_genome(kmer, pos_start) {
                let pos_end = pos_start + self.k - 1;
                let rank2 = self.bv.rank(pos_end).unwrap();
                // is this k-mer in the same unitig as the previous one?
                if rank2 == prev_rank {
                    return Some((pos_start, prev_rank));
                }
            }
        }

        let pos = self.h.get(&Sequence::from_genome(kmer));
        if let Some(pos) = pos {
            if !self.useq.check_genome(kmer, pos) {
                // k-mer not in useq
                return None;
            }
            let rank1 = self.bv.rank(pos).unwrap();
            let rank2 = self.bv.rank(pos + self.k - 1).unwrap();
            if rank1 != rank2 {
                // crossed useq boundary
                return None;
            }
            return Some((pos, rank1));
        }
        None
    }

    pub fn query<S: AsRef<[u8]>>(&self, q: S) -> Vec<&str> {
        let q = q.as_ref();
        let color_bytes = self.n_colors.div_ceil(8);
        let mut colors: Vec<u8> = vec![0xff; color_bytes];
        let mut chunks_iterator = q.chunks_exact(self.k);
        let mut hint: Option<(usize, usize)> = None;
        while let Some(window) = chunks_iterator.next() {
            if let Some((pos, rank)) = self.query_kmer(window, hint) {
                let utab_offset = rank * color_bytes;
                for i in 0..color_bytes {
                    colors[i] &= self.utab[utab_offset + i];
                }
                hint = Some((pos, rank));
            } else {
                return vec![];
            }
        }
        if chunks_iterator.remainder().len() > 0 {
            // check last k-mer
            if let Some((_pos, rank)) = self.query_kmer(&q[q.len() - self.k..], hint) {
                let utab_offset = rank * color_bytes;
                for i in 0..color_bytes {
                    colors[i] &= self.utab[utab_offset + i];
                }
            } else {
                return vec![];
            }
        }
        let mut found_colors = vec![];
        // find colors that describe every k-mer
        for (i, mut c) in colors.into_iter().enumerate() {
            while c > 0 {
                let lowest_one = c.trailing_zeros();
                found_colors.push(&self.identifiers[i * 8 + lowest_one as usize][..]);
                c &= !(1 << lowest_one);
            }
        }
        found_colors
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
        let sequences = vec![("a", "CTAAGAT"), ("b", "CGATGCA"), ("c", "TAAGAGG")];
        let index = HashMapPufferfishIndex::new(3, sequences);

        assert!(index.query("CTAAGAT").contains(&"a"));
        assert!(index.query("CGATGCA").contains(&"b"));
        assert!(index.query("TAAGAGG").contains(&"c"));
    }

    #[test]
    fn basic_mphf_test() {
        let sequences = vec![("a", "CTAAGAT"), ("b", "CGATGCA"), ("c", "TAAGAGG")];
        let index = DefaultPufferfishIndex::new(3, sequences);

        assert!(index.query("CTAAGAT").contains(&"a"));
        assert!(index.query("CGATGCA").contains(&"b"));
        assert!(index.query("TAAGAGG").contains(&"c"));
    }
}
