use std::fmt::{Display, Formatter};
use std::path::PathBuf;
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;

#[derive(Clone)]
pub struct CachedFastaNucltudiesIterator {
    buffer: Vec<u8>,
    buffer_index: usize,
    source: PathBuf,
}

impl From<FastaNucltudiesIterator> for CachedFastaNucltudiesIterator {
    fn from(value: FastaNucltudiesIterator) -> Self {
        let source = value.get_path().to_path_buf();
        Self {
            buffer: value.collect::<Vec<u8>>(),
            buffer_index: 0,
            source,
        }
    }
}

impl Iterator for CachedFastaNucltudiesIterator {
    type Item = u8;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer_index < self.buffer.len() {
            self.buffer_index += 1;
            Some(self.buffer[self.buffer_index - 1])
        } else {
            None
        }
    }
}

impl Display for CachedFastaNucltudiesIterator {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.source.display())
    }
} 

