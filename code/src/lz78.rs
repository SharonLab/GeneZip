//
//  lz78.rs
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;


pub struct LenBases {
    len_bases: Vec<usize>,
}

impl LenBases  {
    pub fn new(max_depth: usize) -> Self {
        let mut len_bases = vec![0_usize; max_depth + 1]; // +1
        // // len_bases is the starting index for each depth i AND the number of indices
        // // required for depth i-1
        // {0, 4, 20, 84, 340, 1364, 5460, 21844, 87380, 349524, 1398100, 5592404, 22369620, 89478484, 357913940, 1431655764, 5726623060, 22906492244};
        for i in 1..=max_depth {
            len_bases[i] = len_bases[i-1] + (4 << (2*(i-1)));
        }

        LenBases {
            len_bases
        }
    }

    pub fn bases(&self) -> &[usize] { &self.len_bases }
}

// Returns the memory size required for keeping all inner nodes up to depth max_depth,
// where max_depth may contain only leaves
fn calc_mem_size(len_bases: &LenBases, max_depth: usize) -> usize {
    (len_bases.bases()[max_depth-1] / 8) + 1
}


////////////////////////////////////////////////////////////////////////////////
// LZ78
////////////////////////////////////////////////////////////////////////////////
pub struct LZ78<'a> {
    mem: Vec<u8>, // Bit array that keeps the inner nodes
    mem_size: usize, // Memory size allocated for mem
    max_depth: usize, // Maximum depth including leaves, specified by the user.
    leaf_count: usize, // Total number of leaves (paths) in the tree. Used for calculating log-loss
    full_depth: usize, // Maximum depth in which all inner nodes are present
    name: String,
    num_nodes_in_depth: Vec<usize>, // Keeps the number of inner nodes in each depth
    len_bases: &'a LenBases,
    buffer_size: usize,
}

impl<'a> LZ78<'a> {
    fn check_bit(&self, idx: usize) -> bool {
        (self.mem[idx >> 3] & (128_u8 >> (idx as u8 & 7))) > 0
    }

    pub fn get_name(&self) -> &str { self.name.as_str() }

    pub fn new(name: &str, max_depth: usize, len_bases: &'a LenBases, file_path: &Path, buffer_size: usize) -> Self {
        let mem_size = calc_mem_size(len_bases, max_depth);

        let mut lz = LZ78 {
            mem: vec![0_u8; mem_size],
            max_depth,
            mem_size,
            name: name.to_string(),
            leaf_count: 4,
            full_depth: 0,
            num_nodes_in_depth: vec![0; max_depth+1],
            len_bases,
            buffer_size,
        };
        lz.num_nodes_in_depth[0] = 1;
        lz.build(file_path);
        lz
    }


    // For a sequence s=s1s2...sn, the index in the bit memory can be calculated
    // using m(s1..si) = (4^0+..+4^(i-1))-1 + 4*m(s1..si-1)
    // There is some bug here
    fn build(&mut self, file_path: &Path) {
        let mut curr_depth = 1; // len is at least 1
        let mut curr_sequence: usize = 0_usize;  // sequence value.

        let mut file_reader = BufReader::new(File::open(file_path)
            .unwrap_or_else(|_| panic!("E: failed to read {}", file_path.display())));

        let mut bytes_buffer = vec![0_u8; self.buffer_size];
        let mut desc_line = false;

        loop {
            match file_reader.read(&mut bytes_buffer) {
                Ok(0) => break,
                Ok(n) => {
                    for p in &bytes_buffer[0..n] {
                        if *p == b'>' {
                            desc_line = true;
                            curr_depth = 1;
                            curr_sequence = 0;
                        } else if desc_line {
                            if *p == b'\n' {
                                desc_line = false;
                            }
                        } else if *p != b'\n' {
                            let p = if *p >= b'a' {
                                p - 32
                            } else {
                                *p
                            };

                            if p == b'N' {
                                curr_depth = 1;
                                curr_sequence = 0;
                                continue;
                            }

                            // For (p >> 1) this is what we get:
                            // A = 0b100000
                            // C = 0b100001
                            // G = 0b100011
                            // T = 0b101010
                            // Order will be ACTG
                            curr_sequence |= ((p >> 1) & 3) as usize;

                            // As far as I understand this can only happen if self.max_depth == 1
                            if curr_depth > self.max_depth-1 {
                                curr_depth = 1;
                                curr_sequence = 0;
                                continue;
                            }

                            //
                            let current_index = self.len_bases.bases()[curr_depth - 1] + curr_sequence;
                            assert!((current_index >> 3) < self.mem_size);
                            if !self.check_bit(current_index) {
                                self.add_node(current_index, curr_depth);
                                curr_depth = 1;
                                curr_sequence = 0;
                                continue;
                            }

                            if curr_depth == self.max_depth - 1 {
                                curr_depth = 1;
                                curr_sequence = 0;
                            } else {
                                curr_sequence <<= 2;
                                curr_depth += 1;
                            }
                        }
                    }
                },
                Err(e) => panic!("E: failed reading file {}, got {:?}", file_path.display(), e),
            };
        }

        // Check what is the maximum level with all nodes
        self.full_depth = 0;
        while self.full_depth + 1 < self.max_depth && self.num_nodes_in_depth[self.full_depth + 1] == (4 << (2 * self.full_depth)) {
            self.full_depth += 1;
        }
    }

    fn add_node(&mut self, node_index: usize, depth: usize) {
        self.mem[node_index >> 3] |= 128 >> (node_index & 7);
        self.num_nodes_in_depth[depth] += 1;

        // One leaf became an inner node, 4 new leaves were added, 3 new leaves added in total
        self.leaf_count += 3;
    }

    pub fn average_log_score(&self, file_path: &Path) -> f64 {
        let mut nchars= 0;
        let mut actual_nchars: usize = 0;
        let mut leaf_count: usize = 0;

        let mut curr_depth: usize = 1;
        let mut curr_sequence: usize = 0;

        let mut file_reader = BufReader::new(File::open(file_path)
            .unwrap_or_else(|_| panic!("E: failed to read {}", file_path.display())));

        let mut bytes_buffer = vec![0_u8; self.buffer_size];
        let mut desc_line = false;

        loop {
            match file_reader.read(&mut bytes_buffer) {
                Ok(0) => break,
                Ok(n) => {
                    for p in &bytes_buffer[0..n] {
                        if *p == b'>' {
                            desc_line = true;
                            curr_depth = 1;
                            curr_sequence = 0;
                        } else if desc_line {
                            if *p == b'\n' {
                                desc_line = false;
                            }
                        } else if *p != b'\n' {
                            // To upper
                            let p = if *p >= b'a' {
                                p - 32
                            } else {
                                *p
                            };
                            if p == b'N' {
                                curr_depth = 1;
                                curr_sequence = 0;
                                continue;
                            }

                            // For (p >> 1) this is what we get:
                            // A = 0b100000
                            // C = 0b100001
                            // G = 0b100011
                            // T = 0b101010
                            // Order will be ACTG
                            let i = (p as usize >> 1) & 3;
                            nchars += 1;
                            curr_sequence |= i;
                            let current_index = self.len_bases.bases()[curr_depth - 1] + curr_sequence;

                            // The first part is an optimization: no need to check if the node
                            // exists for depths with all inner nodes present
                            if curr_depth <= self.full_depth || (curr_depth < self.max_depth && self.check_bit(current_index)) {
                                curr_sequence <<= 2;
                                curr_depth += 1;
                            } else {
                                leaf_count += 1;
                                curr_depth = 1;
                                curr_sequence = 0;
                                actual_nchars = nchars;
                                continue;
                            }
                        }
                    }
                },
                Err(e) => panic!("E: failed reading file {}, got {:?}", file_path.display(), e),
            };
        }

        (self.leaf_count as f64).log2() * leaf_count as f64 / actual_nchars as f64
    }

    fn num_inner_nodes(&self) -> usize {
        (0..self.max_depth)
            .map(|i| self.num_nodes_in_depth[i])
            .sum()
    }

    fn max_complete_depth(&self) -> usize
    {
        self.full_depth
    }

    fn longest_path_root_to_leaf(&self) -> usize {
        let node = self.num_nodes_in_depth.iter()
            .enumerate()
            .find(|(_, &num_nodes_in_depth)| num_nodes_in_depth == 0);

        if let Some((i, _)) = node {
            i
        } else {
            self.max_depth
        }
    }
}

impl<'a> std::fmt::Display for LZ78<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Name:                      {}\n\
                   Node array size:           {}\n\
                   Number of inner nodes:     {}\n\
                   Max complete depth:        {}\n\
                   Longest path (root->leaf): {}\n\
                   Number of inner node in each depth (% of possible nodes):\n\
                   Depth\tNNodes\tNFull\t% of full\n\
                   0\t1\t1\t100.0\n", self.name,
               self.mem_size,
               self.num_inner_nodes(),
               self.max_complete_depth(),
               self.longest_path_root_to_leaf())?;

        for i in 1..self.max_depth {
            let full_n = 4_usize << (2*(i-1));
            writeln!(f, "{}\t{}\t{}\t{:.1}", i,
                                             self.num_nodes_in_depth[i],
                                             full_n,
                                             100.0*(self.num_nodes_in_depth[i] as f64)/(full_n as f64))?;
        }
        writeln!(f, "\nNumber of leaves:\t{}", self.leaf_count)
    }
}

