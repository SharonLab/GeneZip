use std::path::{Path, PathBuf};
use crate::get_gc::calc_gc;
use crate::kmer::create_normalized_profile;
use crate::lz78::{LenBases, LZ78};
use crate::taxonomy::Taxonomy;

use serde::{Serialize, Deserialize};
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;

#[derive(Serialize, Deserialize)]
pub struct ReferenceSequence {
    prediction_model: LZ78,
    gc: f64,
    kmer: Option<ndarray::Array1<f64>>,
    name: String,
    // We want to filter by (k)4mer in some cases, it can be done at genus level, however may be done on other parameters so variable name is kept versatile
    kmer_cluster: Option<Taxonomy>,
    self_value: f64,
    fasta_path: PathBuf,
}

impl ReferenceSequence {
    pub fn new(fasta_path: &Path, name: &str, kmer_size: &Option<usize>, buffer_size: usize, lz_lenbases: LenBases, lzmax_depth: usize, kmer_cluster: &Option<Taxonomy>) -> Self {
        let prediction_model = LZ78::new(lzmax_depth, lz_lenbases, FastaNucltudiesIterator::new(fasta_path, buffer_size));
        let self_value = prediction_model.average_log_score(FastaNucltudiesIterator::new(fasta_path, buffer_size));
        ReferenceSequence {
            prediction_model,
            gc: calc_gc(FastaNucltudiesIterator::new(fasta_path, buffer_size)),
            kmer: kmer_size.map(|k| create_normalized_profile(k, FastaNucltudiesIterator::new(fasta_path, buffer_size), &false).1.unwrap_or_else(|_| panic!("ERROR: failed to create kmer for {}, quitting", fasta_path.display()))),
            name: name.to_string(),
            kmer_cluster: kmer_cluster.clone(),
            self_value,
            fasta_path: fasta_path.to_path_buf(),
        }
    }

    pub fn get_prediction_model(&self) -> &LZ78 { &self.prediction_model }
    pub fn get_gc(&self) -> f64 { self.gc }
    pub fn get_kmer(&self) -> &Option<ndarray::Array1<f64>> { &self.kmer }
    #[allow(dead_code)]
    pub fn get_name(&self) -> &str { self.name.as_str() }
    pub fn get_kmer_cluster(&self) -> Option<&Taxonomy> { self.kmer_cluster.as_ref() }
    pub fn get_self_value(&self) -> f64 { self.self_value }
    pub fn get_fasta_path(&self) -> &Path { &self.fasta_path }
}


impl std::fmt::Display for ReferenceSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Name:                      {}\n\
                   {}
                   0\t1\t1\t100.0\n",
               self.name,
               self.prediction_model)
    }
}