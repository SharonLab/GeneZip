//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::io::{BufWriter, Error, Write};
use std::path::{Path};
use rayon::prelude::*;
use crate::kmer::{create_normalized_profile, create_printed_vector, REV_ALPHABET_VALUES, ALPHABET_MAP};
use crate::samples_file_reader::{SampleError, SampleSource};


pub enum PrintKmersError {
    IoError(std::io::Error),
    SampleError(SampleError)
}

impl From<std::io::Error> for PrintKmersError {
    fn from(value: Error) -> Self {
        PrintKmersError::IoError(value)
    }
}

impl From<SampleError> for PrintKmersError {
    fn from(value: SampleError) -> Self {
        PrintKmersError::SampleError(value)
    }
}

impl Debug for PrintKmersError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            PrintKmersError::IoError(e) => e.fmt(f),
            PrintKmersError::SampleError(e) => e.fmt(f),
        }
    }
}

pub type KMerResults = Result<ndarray::Array1<f64>, String>;
pub fn calc_kmers(input: &Path, k: usize, buffer_size: usize) -> Result<Vec<(String, KMerResults)>, SampleError> {
    Ok(SampleSource::new(input, false)
        .into_iter()
        .collect::<Result<Vec<_>, SampleError>>()?
        .par_iter()
        .map(|sample| (sample.get_name().to_string(), create_normalized_profile(k, sample.get_path(), buffer_size).1))
        .collect())
}

pub fn print_kmers(input: &Path, output: &Path, k: usize, ratio: bool, buffer_size: usize) -> Result<(), PrintKmersError> {
    let mut output_stream = {
        let fout = File::create(output).unwrap_or_else(|_| panic!("E: Cannot create output file '{}'", output.display()));
        BufWriter::new(fout)
    };
    write!(output_stream, "Genome")?;
    for kem in create_printed_vector(k, &ALPHABET_MAP, &REV_ALPHABET_VALUES) {
        write!(output_stream, "\t{kem}")?;
    }
    writeln!(output_stream)?;

    let kmers = calc_kmers(input, k, buffer_size)?;
    for (name, kmer_vec) in kmers {
        match kmer_vec {
            Ok(v) => { write!(output_stream, "{}", name)?;
                if ratio {
                    for freq in v {
                        write!(output_stream, "\t{}", freq)?;
                    }
                } else {
                    for freq in v * 100.0 {
                        write!(output_stream, "\t{}", freq)?;
                    }
                }
                writeln!(output_stream)?;
            },
            Err(e) => eprintln!("ERROR: {}\t{:?}", name, e),
        };
    }

    Ok(())
}