use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use serde::{Serialize, Deserialize};
use flate2::bufread::MultiGzDecoder;
use flate2::{Compression};
use flate2::write::{GzEncoder};
use hashbrown::HashMap;
use crate::print_kmer::calc_kmers;
use crate::samples_file_reader::SampleError;

#[derive(Serialize, Deserialize)]
pub struct KmerDatabase {
    kmer_size: usize,
    kmers: HashMap<String, ndarray::Array1<f64>>
}

fn kmers_vec2hashmap(v: Vec<(String, Result<ndarray::Array1<f64>, String>)>) -> HashMap<String, ndarray::Array1<f64>> {
    let mut results = HashMap::with_capacity(v.len());
    for (name, kmer) in v {
        let k = match kmer {
            Err(e) => panic!("E: Failed to calculate k-mers for '{}', got '{}'", name, e),
            Ok(k) => k,
        };
        results.insert(name, k);
    }

    results
}

impl KmerDatabase {
    pub fn new(name2fasta_path: &Path, buffer_size:usize, kmer_size: usize) -> Result<KmerDatabase, SampleError> {
        Ok(KmerDatabase {
            kmer_size,
            kmers: kmers_vec2hashmap(calc_kmers(name2fasta_path, kmer_size, buffer_size)?),
        })
    }

    pub fn save(&self, destination: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let mut gzw = GzEncoder::new(BufWriter::new(File::create(destination)?), Compression::best());
        bincode::serialize_into(&mut gzw, &self)?;
        gzw.finish()?.flush()?;
        Ok(())
    }

    #[allow(dead_code)]
    pub fn load(source: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        let bin_self_reader = MultiGzDecoder::new(BufReader::new(File::open(source)?));
        Ok(bincode::deserialize_from(bin_self_reader)?)
    }

    #[allow(dead_code)]
    pub fn get_kmer_size(&self) -> usize { self.kmer_size }
    #[allow(dead_code)]
    pub fn get_kmers(&self) -> &HashMap<String, ndarray::Array1<f64>> { &self.kmers }
}