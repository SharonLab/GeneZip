//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use crate::classifier::Classifier;
use serde::{Serialize, Deserialize};
use gzp::{ZWriter, deflate::Bgzf, Compression};
use gzp::par::compress::{ParCompress, ParCompressBuilder};
use gzp::par::decompress::{ParDecompressBuilder, ParDecompress};

#[derive(Serialize, Deserialize)]
pub struct Database {
    classifier: Classifier,
    max_depth: usize,
    kmer_size: Option<usize>,
}

impl Database {
    pub fn new(classifier: Classifier, max_depth: usize, kmer_size: Option<usize>) -> Database {
        Database {
            classifier,
            max_depth,
            kmer_size
        }
    }

    pub fn save(&self, destination: &Path, threads_limit: usize) -> Result<(), Box<dyn std::error::Error>> {
        let mut gzw: ParCompress<Bgzf> = ParCompressBuilder::new()
            .compression_level(Compression::new(9))
            .num_threads(threads_limit)?
            .from_writer(BufWriter::new(File::create(destination)?));
        bincode::serialize_into(&mut gzw, &self)?;
        gzw.finish()?;
        Ok(())
    }

    pub fn load(source: &Path, threads_limit: usize) -> Result<Self, Box<dyn std::error::Error>> {
        let bin_self_reader: ParDecompress<Bgzf> = ParDecompressBuilder::new()
            .num_threads(threads_limit)?
            .from_reader(BufReader::new(File::open(source)?));
        Ok(bincode::deserialize_from(bin_self_reader)?)
    }

    #[allow(dead_code)]
    pub fn get_max_depth(&self) -> usize { self.max_depth }
    pub fn get_kmer_size(&self) -> &Option<usize> { &self.kmer_size }
    pub fn get_classifier(&self) -> &Classifier { &self.classifier }
}