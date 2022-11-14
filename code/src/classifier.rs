//
//  classifier.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use hashbrown::HashMap;

use crate::lz78::{LZ78, LenBases};

pub struct Classifier<'a>  {
    models: Box<HashMap<String, LZ78<'a>>>,
    len_bases: &'a LenBases,
    models_order: Vec<String>,
}

impl<'a> Classifier<'a> {
    pub fn new(len_bases: &'a LenBases) -> Self {
        Classifier {
            models: Box::new(HashMap::new()),
            len_bases,
            models_order: Vec::new(),
        }
    }

    pub fn add_model(&mut self, name: &str, model: LZ78<'a>) {
        self.models.insert(name.to_string(), model);
        self.models_order.push(name.to_string());
        self.models_order.sort();
    }

    pub fn batch_add_model(&mut self, name2file: &Path, max_depth: usize, buffer_size: usize) {
        let file_reader = BufReader::new(File::open(name2file)
            .unwrap_or_else(|_| panic!("E: failed to read {}", name2file.display())));

        let mut models: Vec<_> = file_reader.lines().map(|l| l.unwrap_or_else(|_| panic!("E: Illegal line in file {}", name2file.display())))
            .into_iter()
            .par_bridge()
            .map(|line| {
                let mut split = line.split('\t');
                let name = split.next().unwrap_or_else(|| panic!("E: invalid line in file {}", name2file.display()));
                let file_path = split.next().expect("E: invalid line");

                (name.to_string(), Some(LZ78::new(name, max_depth, self.len_bases, &PathBuf::from(file_path), buffer_size)))
            })
            .collect();

        models.iter_mut()
            .for_each(|(name, model)| self.add_model(name.as_str(), model.take().unwrap()));
    }

    pub fn print_prediction<W: Write>(&self, name: &str, fout: &mut BufWriter<W>, prediction: &Vec<(&String, f64)>) -> std::io::Result<()> {
        write!(fout, "{}", name)?;

        let mut model2score = HashMap::new();
        let best_model_name = self.get_best_model_name(prediction);

        for &(model_name, score) in prediction {
            model2score.insert(model_name, score);
        }

        for model_name in self.models_order.iter() {
            write!(fout, "\t{:.5}", model2score.get(model_name).unwrap_or(&f64::NAN))?;
        }

        if let Some(best_model_name) = best_model_name {
            writeln!(fout, "\t{}", best_model_name)?;
        }

        Ok(())
    }

    pub fn get_best_model_name(&self, prediction: &Vec<(&'a String, f64)>) -> Option<&'a String> {
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

    pub fn predict(&self, file_path: &Path) -> Vec<(&String, f64)> {
        let model_name2score: Vec<_> = self.models.iter()
            .par_bridge()
            .map(|(model_name, model)| {
                let score = model.average_log_score(file_path);
                (model_name, score)
            } )
            .collect();
        model_name2score
    }

    pub fn print_header<W: Write>(&self, fout: &mut BufWriter<W>) -> std::io::Result<()> {
        write!(fout, "Genome_name")?;

        for model_name in self.models_order.iter() {
            write!(fout, "\t{}", model_name)?;
        }

        writeln!(fout, "\tBest_hit")
    }

    pub fn print_stats<W: Write>(&self, fout: &mut BufWriter<W>) -> std::io::Result<()> {
        write!(fout, "\nNumber of models: {}\n", self.models.len())?;
        write!(fout, "some stats for each model:\n\n")?;

        for model in self.models.values() {
            writeln!(fout, "--------------------------------------------------")?;
            write!(fout, "{}", model)?;
        }
        writeln!(fout)
    }
}
