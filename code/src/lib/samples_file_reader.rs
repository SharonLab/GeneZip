//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use crate::reference_sequence::ReferenceSequence;
use crate::taxonomy::Taxonomy;

#[derive(Debug)]
#[allow(dead_code)] // IoError(std::io::Error) will raise a warning, as the error data isn't yet used
pub enum SampleErrorType {
    IoError(std::io::Error),
    None,
}

#[derive(Debug)]
#[allow(dead_code)]
pub struct SampleError {
    context: String,
    error: SampleErrorType,
}

impl SampleError {
    pub fn new(context: &str, error: SampleErrorType) -> Self {
        SampleError {
            context: context.to_string(),
            error,
        }
    }
}

impl Display for SampleError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.context)
    }
}

#[derive(Clone)]
pub struct Sample {
    name: String,
    path: PathBuf,
    taxonomy: Option<Taxonomy>,
    line_number: usize,
}

impl Sample {
    pub fn new(name: &str, path: &str, taxonomy: Option<&str>, line_number: usize) -> Self {
        Sample {
            name: name.to_string(),
            path: PathBuf::from(path),
            taxonomy: taxonomy.map(Taxonomy::from),
            line_number,
        }
    }

    pub fn get_name(&self) -> &str { &self.name }
    pub fn get_path(&self) -> &Path { &self.path }
    pub fn get_taxonomy(&self) -> &Option<Taxonomy> { &self.taxonomy }
    pub fn get_line_number(&self) -> usize { self.line_number }
}

impl Debug for Sample {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sample({}, {}, {:?}, {})", self.get_name(), self.get_path().display(), self.get_taxonomy(), self.get_line_number())
    }
}
impl From<&ReferenceSequence> for Sample {
    fn from(value: &ReferenceSequence) -> Self {
        Self {
            name: value.get_name().to_string(),
            path: value.get_fasta_path().to_path_buf(),
            taxonomy: value.get_kmer_cluster().cloned(),
            line_number: 0,
        }
    }
}

impl Ord for Sample {
    fn cmp(&self, other: &Self) -> Ordering {
        self.get_path().cmp(other.get_path())
    }
}

impl PartialOrd for Sample {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Sample {
    fn eq(&self, other: &Self) -> bool {
        self.get_path().eq(other.get_path())
    }
}
impl Eq for Sample {

}

impl std::hash::Hash for Sample {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.get_path().hash(state)
    }
}

pub struct SampleSource {
    path: PathBuf,
    include_taxonomy: bool,
}

impl SampleSource {
    pub fn new(path: &Path, include_taxonomy: bool) -> Self {
        Self {
            path: path.to_path_buf(),
            include_taxonomy,
        }
    }

    fn get_path(&self) -> &Path { &self.path }
    fn get_include_taxonomy(&self) -> bool { self.include_taxonomy }
}

impl IntoIterator for SampleSource {
    type Item = Result<Sample, SampleError>;
    type IntoIter = SampleIterator;

    fn into_iter(self) -> Self::IntoIter {
        SampleIterator::from(&self)
    }
}

pub struct SampleIterator {
    source: PathBuf,
    reader: Box<dyn BufRead>,
    include_taxonomy: bool,
    line_number: usize,
}

impl From<SampleSource> for SampleIterator {
    fn from(value: SampleSource) -> Self {
        SampleIterator::from(&value)
    }
}
impl From<&SampleSource> for SampleIterator {
    fn from(value: &SampleSource) -> Self {
        SampleIterator {
            source: value.path.to_path_buf(),
            reader: Box::new(BufReader::new(match File::open(value.get_path()) {
                Ok(f) => f,
                Err(e) => panic!("E: Tried to open '{}' to read samples from, but encountered the following error: '{}'", value.get_path().display(), e),
            })),
            include_taxonomy: value.get_include_taxonomy(),
            line_number: 0,
        }
    }
}

impl SampleIterator {
    fn next(&mut self) -> Option<Result<Sample, SampleError>> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Err(e) => Some(Err(SampleError::new(&format!("E: Trying to parse samples from '{}' led to '{}'", self.source.display(), e), SampleErrorType::IoError(e)))),
            Ok(0) => None,
            Ok(_) => {
                let mut split = line.strip_suffix('\n').unwrap().split('\t');
                let name = match split.next() {
                    Some(n) => n,
                    None => return Some(Err(SampleError::new(&format!("E: Invalid line (number {}) in file {}, expected a tab-delimited line.", self.line_number, self.source.display()), SampleErrorType::None))),
                };
                let file_path = match split.next() {
                    Some(fp) => fp,
                    None => return Some(Err(SampleError::new(&format!("E: Invalid line (number {}) in file {}, expected a tab-delimited line with at-least 2 fields.", self.line_number, self.source.display()), SampleErrorType::None)))
                };
                let taxonomy = if self.include_taxonomy {
                    match split.next() {
                        Some(t) => Some(t),
                        None => return Some(Err(SampleError::new(&format!("E: Invalid line (number {}) in file {}, expected a tab-delimited line with at-least 3 fields. This may happen if taxonomy column is missing, you may set --kmer 0 or add the missing column.", self.line_number, self.source.display()), SampleErrorType::None))),
                    }
                } else {
                    None
                };
                self.line_number += 1;
                Some(Ok(Sample::new(name, file_path, taxonomy, self.line_number - 1)))
            }
        }
    }
}

impl Iterator for SampleIterator {
    type Item = Result<Sample, SampleError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next()
    }
}


// impl rayon::iter::ParallelIterator for SampleIterator {
//     type Item = Result<Sample, SampleError>;
//
//     fn drive_unindexed<C>(mut self, consumer: C) -> C::Result
//     where
//         C: UnindexedConsumer<Self::Item>,
//         C::Result: Result<Sample, SampleError>
//     {
//         self.next().unwrap()
//     }
// }