//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::path::Path;
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;

pub enum FastaPartType {
    ID(Vec<u8>),
    Nuc(u8)
}

pub struct FastaRecordIterator {
    fni: FastaNucltudiesIterator,
    buffer_size: usize,
}

impl FastaRecordIterator {
    pub fn new(fasta: &Path, buffer_size: usize) -> Self {
        Self {
            fni: FastaNucltudiesIterator::new(fasta, buffer_size),
            buffer_size,
        }
    }

    fn next_component(&mut self) -> Option<FastaPartType> {
        match self.fni.next_char() {
            None => None,
            Some(c) => {
                if c == b'>' {
                    let mut id = Vec::with_capacity(self.buffer_size);
                    let mut got_all_id = false;
                    while let Some(c) = self.fni.next_char() {
                        if [b' ', b'\t'].contains(&c) {
                            got_all_id = true;
                        }
                        if [b'\n', b'\r'].contains(&c) {
                            break
                        }
                        if ! got_all_id {
                            id.push(c);
                        }
                    }
                    Some(FastaPartType::ID(id))
                } else if c == b'\n' {
                    self.next_component() // Recursive call
                } else {
                    Some(FastaPartType::Nuc(c & 0b11011111)) // this maps low case letters to capital, and capital to self.
                }
            }
        }
    }
}

impl Iterator for FastaRecordIterator {
    type Item = FastaPartType;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_component()
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::fasta_records_iterator::{FastaPartType, FastaRecordIterator};

    #[test]
    fn use_gc_test() {
        let correct_nucs = [b'A', b'T', b'A', b'T', b'G', b'G', b'C', b'T', b'A', b'T', b'A', b'T', b'G', b'C', b'T', b'A', b'T', b'A', b'G', b'C', b'A', b'T', b'T', b'T', b'T', b'G', b'G', b'T', b'G', b'A', b'C', b'G', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'A', b'T', b'G'];
        let correct_ids = ["1", "2"];

        let mut nucs_iter = correct_nucs.into_iter();
        let mut ids_iter = correct_ids.into_iter();

        let fasta = PathBuf::from("../tests/gc_tests.fna");

        for component in FastaRecordIterator::new(&fasta, 512) {
            match component {
                FastaPartType::Nuc(n) => assert_eq!(n, nucs_iter.next().unwrap()),
                FastaPartType::ID(id) => assert_eq!(std::str::from_utf8(id.as_slice()).unwrap(), ids_iter.next().unwrap()),
            }
        }

        assert_eq!(None, nucs_iter.next());
        assert_eq!(None, ids_iter.next());
    }

    #[test]
    fn presentation_test_add_n() {
        let correct_nucs = [b'C', b'G', b'G', b'T', b'N', b'T'];
        let correct_ids = ["test"];

        let mut nucs_iter = correct_nucs.into_iter();
        let mut ids_iter = correct_ids.into_iter();

        let fasta = PathBuf::from("../tests/presentation_test_add_N.fna");

        for component in FastaRecordIterator::new(&fasta, 512) {
            match component {
                FastaPartType::Nuc(n) => assert_eq!(n, nucs_iter.next().unwrap()),
                FastaPartType::ID(id) => assert_eq!(std::str::from_utf8(id.as_slice()).unwrap(), ids_iter.next().unwrap()),
            }
        }

        assert_eq!(None, nucs_iter.next());
        assert_eq!(None, ids_iter.next());
    }
}