use std::io::{BufWriter, Write};
use std::path::{Path};
use hashbrown::HashMap;
use ndarray::stack;
use ndarray_stats::CorrelationExt;
use rayon::iter::{IntoParallelRefIterator, ParallelBridge};
use rayon::iter::ParallelIterator;
use crate::kmer::create_normalized_profile;
use crate::samples_file_reader;
use crate::samples_file_reader::SampleError;

pub struct KmerClassifier {
    models: Box<HashMap<String, ndarray::Array1<f64>>>,
    models_order: Vec<String>,
    kmer_size: usize,
}

fn rate(kmer1: &ndarray::Array1<f64>, kmer2: &ndarray::Array1<f64>) -> f64 {
    let merged_array = stack!(ndarray::Axis(0), kmer1.clone(), kmer2.clone());
    match merged_array.pearson_correlation() {
        Ok(c) => c.mean().unwrap(),
        Err(e) => panic!("E: tried to calculate the correlation between genome to mode, got {:?}, quitting", e),
    }
}

impl KmerClassifier {
    pub fn new(kmer_size: usize) -> Self {
        KmerClassifier {
            models: Box::new(HashMap::new()),
            models_order: Vec::new(),
            kmer_size,
        }
    }

    pub fn add_model(&mut self, model: &mut ndarray::Array1<f64>, name: &str) {
        self.models_order.push(name.to_string());
        self.models.insert(name.to_string(), model.clone());
    }

    pub fn batch_add_model(&mut self, name2file: &Path, buffer_size: usize) -> Result<(), SampleError> {
        let mut models: Vec<_> = samples_file_reader::SampleSource::new(name2file, false)
            .into_iter()
            .collect::<Result<Vec<_>,SampleError>>()?
            .par_iter()
            .map(|sample| (sample.get_name().to_string(), create_normalized_profile(self.kmer_size, sample.get_path(), buffer_size).1.unwrap_or_else(|_| panic!("ERROR: failed to create kmer for {}, quitting", sample.get_path().display()))))
            .collect();

        models.iter_mut()
            .for_each(|(name, model)| self.add_model(model, name.as_str()));
        Ok(())
    }

    pub fn print_prediction<W: Write>(&self, name: &str, fout: &mut BufWriter<W>, prediction: &Vec<(&String, f64)>) -> std::io::Result<()> {
        write!(fout, "{name}")?;

        let mut model2score = HashMap::new();
        let best_model_name = Self::get_best_model_name(prediction);

        for &(model_name, score) in prediction {
            model2score.insert(model_name, score);
        }

        for model_name in &self.models_order {
            let score = model2score[model_name];
            write!(fout, "\t{score:.5}")?
        }

        if let Some(best_model_name) = best_model_name {
            writeln!(fout, "\t{best_model_name}")?;
        }

        Ok(())
    }

    pub fn get_best_model_name<'a>(prediction: &Vec<(&'a String, f64)>) -> Option<&'a String> {
        let mut best_model_name = None;
        let mut best_score = None;

        for &(model_name, score) in prediction {
            if best_score.is_none() || best_score.unwrap() > score {
                best_score = Some(score);
                best_model_name = Some(model_name);
            }
        }

        best_model_name
    }

    pub fn predict(&self, file_path: &Path, buffer_size: usize) -> Vec<(&String, f64)> {
        let genome_kmer = match create_normalized_profile(self.kmer_size, file_path, buffer_size).1 {
            Ok(vector) => vector,
            Err(e) => panic!("E: tried to create k({})-mer for genome {}, but got {:?}, quitting", self.kmer_size, file_path.display(), e),
        };

        self.models_order.iter()
            .par_bridge()
            .map(|model_name| {
                let model = self.models.get(model_name).unwrap();
                (model_name, rate(model, &genome_kmer))
            } )
            .collect::<Vec<_>>()
    }

    pub fn print_header<W: Write>(&self, fout: &mut BufWriter<W>) -> std::io::Result<()> {
        write!(fout, "Genome_name")?;

        for model_name in &self.models_order {
            write!(fout, "\t{model_name}")?;
        }

        writeln!(fout, "\tBest_hit")
    }
}