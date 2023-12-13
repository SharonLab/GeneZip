//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::io::{BufWriter, Write};
use std::path::{Path};
use std::string::String;
use std::collections::HashSet;

use rayon::prelude::*;
use hashbrown::{HashMap};
use ndarray::{Array1, stack};
use ndarray_stats::CorrelationExt;

use crate::get_gc::calc_gc;
use crate::kmer::create_normalized_profile;

use crate::lz78::{LenBases, LZ78};
use crate::reference_sequence::ReferenceSequence;
use crate::taxonomy::{TaxonomicRank, Taxonomy};
use serde::{Serialize, Deserialize};
use crate::samples_file_reader::{Sample, SampleError, SampleSource};

#[derive(Serialize, Deserialize)]
pub struct Classifier {
    models: Box<HashMap<String, ReferenceSequence>>,
    len_bases: LenBases,
    models_order: Vec<String>,
}

impl Classifier {
    pub fn new(len_bases: LenBases) -> Self {
        Classifier {
            models: Box::new(HashMap::new()),
            len_bases,
            models_order: Vec::new(),
        }
    }

    pub fn add_model(&mut self, name: &str, model: ReferenceSequence) {
        self.models.insert(name.to_string(), model);
        self.models_order.push(name.to_string());
        self.models_order.sort();
    }

    pub fn batch_add_model(&mut self, name2file: &Path, max_depth: usize, buffer_size: usize, kmer_size: &Option<usize>) -> Result<(), SampleError> {
        let samples = SampleSource::new(name2file, kmer_size.is_some())
            .into_iter()
            .collect::<Result<Vec<Sample>, SampleError>>()?;
        let mut models: Vec<_> = samples
            .par_iter()
            .map(|sample|  (sample.get_name().to_string(), Some(ReferenceSequence::new(sample.get_path(), sample.get_name(), kmer_size, buffer_size, self.len_bases.clone(), max_depth, sample.get_taxonomy()))))
            .collect();

        models.iter_mut()
            .for_each(|(name, model)| self.add_model(name.as_str(), model.take().unwrap()));
        Ok(())
    }

    pub fn print_prediction<W: Write>(&self, name: &str, fout: &mut BufWriter<W>, prediction: &Vec<(&String, Option<f64>)>) -> std::io::Result<()> {
        write!(fout, "{name}")?;

        let mut model2score = HashMap::new();
        let best_model_name = Classifier::get_best_model_name(prediction);

        for &(model_name, score) in prediction {
            model2score.insert(model_name, score);
        }

        for model_name in &self.models_order {
            let score = model2score.get(model_name).unwrap_or(&None);
            match score {
                Some(score) => write!(fout, "\t{score:.5}")?,
                None => write!(fout, "\tNA")?,
            }
        }

        if let Some(best_model_name) = best_model_name {
            writeln!(fout, "\t{best_model_name}")?;
        }

        Ok(())
    }

    pub fn get_best_model_name<'a>(prediction: &Vec<(&'a String, Option<f64>)>) -> Option<&'a String> {
        let mut best_model_name = None;
        let mut best_score = None;

        for &(model_name, score) in prediction {
            if let Some(score) = score {
                if best_score.is_none() || best_score.unwrap() > score {
                    best_score = Some(score);
                    best_model_name = Some(model_name);
                }
            }
        }

        best_model_name
    }


    fn filter_models_by_gc<'a>(&self, file_path: &Path, gc_limit: f64, buffer_size: usize, models_to_check: &HashSet<&'a String>) -> HashSet<&'a String> {
        let genome_gc = calc_gc(file_path, buffer_size);

        models_to_check.into_par_iter()
            .filter_map(|&model_name| {
                if (self.models.get(model_name).unwrap().get_gc() - genome_gc).abs() < gc_limit {
                    Some(model_name)
                } else {
                    None
                }
            } )
            .collect::<HashSet<_>>()
    }

    fn get_model_kmer_correlation(&self, model_name: &String, genome_kmer: &Array1<f64>, file_path: &Path) -> (&Taxonomy, f64) {
        let model_kmer_vector = match self.models.get(model_name) {
            None => panic!("E: Tried to access non-existing model by its name, this should never happen, quitting"),
            Some(reference_sequence) => match reference_sequence.get_kmer() {
                Some(kmer_vecotr) => kmer_vecotr,
                None => panic!("E: Tried to access kmer vector for model {}, however, this model has no kmer vector, quitting.", model_name),
            }
        };

        if model_kmer_vector.len() != genome_kmer.len() {
            panic!("E: Model {} and genome {} have different k for k-mer filtration, quitting", model_name, file_path.display());
        }

        // let merged_array = arr2(&vec![genome_kmer, model_kmer_vector]);
        let merged_array = stack!(ndarray::Axis(0), genome_kmer.clone(), model_kmer_vector.clone());
        let correlation = match merged_array.pearson_correlation() {
            Ok(c) => c.mean().unwrap(),  // TODO: make sure this is indead the coefficient
            Err(e) => panic!("E: tried to calculate the correlation between genome {} and model {}, got {:?}, quitting", file_path.display(), model_name, e),
        };

        (self.models[model_name].get_kmer_cluster().unwrap_or_else(|| panic!("E: tried to filter by kmer but model {} has no kmer cluster, quitting", model_name)),
         correlation)
    }

    fn filter_models_by_kmer<'b>(&self, file_path: &Path, kmer_cluster_limit: usize, buffer_size: usize, models_to_check: &HashSet<&'b String>) -> HashSet<&'b String> {
        let genome_kmer = match create_normalized_profile(kmer_cluster_limit, file_path, buffer_size).1 {
            Ok(vector) => vector,
            Err(e) => panic!("E: tried to create k({})-mer for genome {}, but got {:?}, quitting", kmer_cluster_limit, file_path.display(), e),
        };

        let best_kmer_cluster = models_to_check.into_par_iter()
            .map(|&model_name | self.get_model_kmer_correlation(model_name, &genome_kmer, file_path))
            .max_by(|pair_a, pair_b| pair_a.1.total_cmp(&pair_b.1) ).unwrap().0;

        models_to_check.iter()
            .filter(|&&model_name| self.models[model_name].get_kmer_cluster().unwrap_or_else(||panic!("E: tried to filter by kmer but model {} has no kmer cluster, quitting", model_name)).equal_to_rank(best_kmer_cluster, &TaxonomicRank::Genus))
            .copied()
            .collect::<HashSet<&String>>()
    }


    pub fn predict(&self, file_path: &Path, gc_limit: Option<f64>, buffer_size: usize, kmer_cluster_limit: &Option<usize>, reflect: bool) -> Vec<(&String, Option<f64>)> {
        let models_to_check = match gc_limit {
            None => HashSet::from_iter(self.models_order.iter()),
            Some(gc_limit) => self.filter_models_by_gc(file_path, gc_limit, buffer_size, &HashSet::from_iter(&self.models_order)),
        };

        let models_to_check = match kmer_cluster_limit {
            None => models_to_check,
            Some(kmer_limit) => self.filter_models_by_kmer(file_path, *kmer_limit, buffer_size, &models_to_check),
        };

        let self_reflection = if reflect {
            let model = LZ78::new(self.len_bases.get_max_depth(), self.len_bases.clone(), file_path, buffer_size);
            let self_value = model.average_log_score(file_path);
            Some((model, self_value))
        } else { None };

        models_to_check.iter()
            .par_bridge()
            .map(|&model_name| {
                let model = self.models.get(model_name).unwrap();
                let gz_model_genome = model.get_prediction_model().average_log_score(file_path);
                if let Some((self_model, self_value)) = &self_reflection {
                    (model_name, Some( (gz_model_genome + self_model.average_log_score(model.get_fasta_path()) ) / ( self_value + model.get_self_value() )  ))
                } else {
                    (model_name, Some(gz_model_genome))
                }

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

    pub fn print_stats<W: Write>(&self, fout: &mut BufWriter<W>) -> std::io::Result<()> {
        write!(fout, "\nNumber of models: {}\n", self.models.len())?;
        write!(fout, "some stats for each model:\n\n")?;

        for model in self.models.values() {
            writeln!(fout, "--------------------------------------------------")?;
            write!(fout, "{model}")?;
        }
        writeln!(fout)
    }

    fn iter_references<'a>(&'a self) -> Box<dyn Iterator<Item=& ReferenceSequence> + 'a> {
        Box::new(self.models.values())
    }
}

impl<'a> IntoIterator for &'a Classifier {
    type Item = Sample;
    type IntoIter = Box<dyn Iterator<Item=Sample> + 'a>;
    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter_references().map(Sample::from))
    }
}