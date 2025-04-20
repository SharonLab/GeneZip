//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::ffi::OsStr;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use flate2::bufread::MultiGzDecoder;

pub struct FastaNucltudiesIterator {
    path: PathBuf,
    stream: Box<dyn BufRead + Sync>,
    buffer: Vec<u8>,
    buffer_index: usize,
    buffer_filled: usize,
}

impl Clone for FastaNucltudiesIterator {
    fn clone(&self) -> Self {
        Self::new(&self.path, self.buffer.len())
    }
}

impl FastaNucltudiesIterator {
    pub fn new(fasta: &Path, buffer_size: usize) -> Self {
        FastaNucltudiesIterator {
            path: fasta.to_path_buf(),
            stream: Self::open_fasta(fasta),
            buffer: vec![0_u8; buffer_size],
            buffer_index: 0,
            buffer_filled: 0,
        }
    }

    fn open_fasta(fasta: &Path) -> Box<dyn BufRead + Sync> {
        let f = BufReader::new(match File::open(fasta) {
            Err(e) => panic!("E: Can not open fasta file at '{}', got '{}'", fasta.display(), e),
            Ok(f) => f,
        });
        if Some(OsStr::new("gz")) == fasta.extension() {
            Box::new(BufReader::new(MultiGzDecoder::new(f)))
        } else {
            Box::new(f)
        }
    }

    #[inline]
    pub fn next_char(&mut self) -> Option<u8> {
        if self.buffer_filled <= self.buffer_index {
            match self.stream.read(&mut self.buffer) {
                Ok(n) => {
                    self.buffer_index = 0;
                    self.buffer_filled = n;
                },
                Err(e) => panic!("ERROR: Failed reading fasta file at '{}', got {:?}", self.path.display(), e),
            }
        }

        if self.buffer_filled > self.buffer_index {
            let byte = Some(self.buffer[self.buffer_index]);
            self.buffer_index += 1;
            byte
        } else {
            None
        }
    }

    #[inline]
    fn next_nuc(&mut self) -> Option<u8> {
        match self.next_char() {
            None => None,
            Some(c) => {
                if c == b'>' {
                    while let Some(c) = self.next_char() {
                        if c == b'\n' {
                            break
                        }
                    }
                    Some(b'N') // This marks the skip between sequences
                } else if c == b'\n' {
                    self.next_nuc() // Recursive call
                } else {
                    Some(c & 0b11011111) // this maps low case letters to capital, and capital to self.
                }
            }
        }
    }

    #[allow(dead_code)]
    pub fn read(&mut self, buffer: &mut Vec<u8>) -> usize {
        let mut results = 0_usize;
        for pos in buffer {
            match self.next_nuc() {
                Some(nuc) => {
                    *pos = nuc;
                    results += 1;
                },
                None => break,
            }
        }
        results
    }
    
    pub fn get_path(&self) -> &Path { self.path.as_path() }
}

impl Iterator for FastaNucltudiesIterator {
    type Item = u8;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.next_nuc()
    }
}

impl Display for FastaNucltudiesIterator {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.path.display())
    }
}