use crate::{
    bitvector::BitVector,
    bloom_filter::BloomFilter,
    mphf::MPHF,
    util::{HyperLogLog, MapLike, Sequence, index_in_acgt},
};
use log::{debug, trace};
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

        // convert to 2-bit representation
        let reference_strings: Vec<Vec<u8>> = reference_strings
            .into_iter()
            .map(|s| {
                s.as_ref()
                    .into_iter()
                    .map(|&c| index_in_acgt(c) as u8)
                    .collect()
            })
            .collect();

        trace!("to 2-bit repr: {}ms", (Instant::now() - start).as_millis());

        let start = Instant::now();

        let mut unique_edges_sketch = HyperLogLog::new(8);

        for string in reference_strings.iter() {
            for window in string.windows(k) {
                // insert edge
                unique_edges_sketch.insert(&window);
            }
        }

        let estimated_unique_edges = unique_edges_sketch.count();

        debug!(
            "approx. count unique edges: {}ms",
            (Instant::now() - start).as_millis()
        );

        let start = Instant::now();

        // Build De Bruijn graph
        let mut starts: HashSet<&[u8]> = HashSet::new();
        let mut ends: HashSet<&[u8]> = HashSet::new();

        let mut edges_filter = BloomFilter::with_fpr(0.01, estimated_unique_edges);

        for string in reference_strings.iter() {
            // mark starts/ends
            starts.insert(&string[..k]);
            ends.insert(&string[string.len() - k..]);

            for window in string.windows(k + 1) {
                // insert edge
                edges_filter.insert_key(&window);
            }
        }

        trace!(
            "build bloom filter: {}ms",
            (Instant::now() - start).as_millis()
        );

        let mut critical_false_positives: HashSet<Vec<u8>> = HashSet::new();
        let mut maybe_branch: HashSet<&[u8]> = HashSet::new();
        let mut maybe_multiple_back_edges: HashSet<&[u8]> = HashSet::new();

        let start = Instant::now();

        for string in reference_strings.iter() {
            // check backwards edges from start
            let mut key: Vec<u8> = vec![0];
            key.extend(&string[..k]);

            for c in [0, 1, 2, 3] {
                key[0] = c;
                if edges_filter.query_key(&key) {
                    critical_false_positives.insert(key.clone());
                }
            }

            for window in string.windows(k + 1) {
                // check other forward edges
                let mut key: Vec<u8> = window.to_vec();
                for c in [1, 2, 3] {
                    key[k] = window[k] ^ c;
                    if edges_filter.query_key(&key) {
                        critical_false_positives.insert(key.clone());
                        maybe_branch.insert(&window[..k]);
                    }
                }

                // reset key
                key[k] = window[k];

                // check other backwards edges
                for c in [1, 2, 3] {
                    key[0] = window[0] ^ c;
                    if edges_filter.query_key(&key) {
                        critical_false_positives.insert(key.clone());
                        maybe_multiple_back_edges.insert(&window[1..]);
                    }
                }
            }

            // check forward edges from end
            let mut key: Vec<u8> = string[string.len() - k..].to_vec();
            key.push(0);

            for c in [0, 1, 2, 3] {
                key[k] = c;
                if edges_filter.query_key(&key) {
                    critical_false_positives.insert(key.clone());
                }
            }
        }

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
                critical_false_positives.remove(&window.to_vec());
            }
        }

        trace!(
            "identify true positives: {}ms",
            (Instant::now() - start).as_millis()
        );
        debug!("{} cfps", critical_false_positives.len());

        let is_real_edge = |key: &[u8]| -> bool {
            edges_filter.query_key(key.as_ref()) && !critical_false_positives.contains(key.as_ref())
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

        // Find junctions (all places where a unipath starts)
        // Hashmap will store color information
        let mut junctions: HashMap<Vec<u8>, Vec<u8>> =
            HashMap::from_iter(starts.iter().map(|&k| (k.to_vec(), vec![0; color_bytes])));

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

        for node in maybe_multiple_back_edges {
            let mut key: Vec<u8> = vec![0];
            key.extend(node);
            let mut count = 0;
            for c in [0, 1, 2, 3] {
                key[0] = c;
                if is_real_edge(&key) {
                    count += 1;
                    if count > 1 {
                        // multiple back edges is a junction
                        junctions.insert(node.to_vec(), vec![0; color_bytes]);
                        break;
                    }
                }
            }
        }

        trace!("find junctions: {}ms", (Instant::now() - start).as_millis());
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

        let mut useq: Sequence = Sequence::new();
        let mut pos: Vec<(Sequence, usize)> = Vec::new();
        let mut bv: Vec<u64> = Vec::new();
        let mut utab: Vec<u8> = Vec::new();

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

                if junctions.contains_key(&key[..]) {
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
