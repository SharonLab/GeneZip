//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::fmt::{Debug, Formatter};
use std::fs::File;
use std::io::{BufWriter, Error, Write};
use std::path::Path;

use rayon::prelude::*;
use tempdir::TempDir;

use crate::contig_naming::sequence_id2str;
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;
use crate::fasta_records_iterator;
use crate::kmer::{ALPHABET_MAP, create_normalized_profile, create_printed_vector, REV_ALPHABET_VALUES};
use crate::samples_file_reader::{SampleError, SampleSource};

pub enum PrintKmersError {
    IoError(Error),
    SampleError(SampleError)
}

impl From<Error> for PrintKmersError {
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
        .map(|sample| (sample.get_name().to_string(), create_normalized_profile(k, FastaNucltudiesIterator::new(sample.get_path(), buffer_size), &false).1))
        .collect())
}

pub fn calc_kmers_meta(fasta: &Path, k: usize, buffer_size: usize) -> Vec<(String, KMerResults)> {
    let work_dir = TempDir::new("genezip").expect("ERROR: failed to create a temporary folder, your tmp may be full, quitting!");
    let temp_fasta_path = work_dir.path().join("temp.fasta");
    let mut prev_record_id: Option<Vec<u8>> = None;
    let mut temp_fasta_stream: Option<BufWriter<File>> = None;
    let mut results = Vec::new();
    for record_part in fasta_records_iterator::FastaRecordIterator::new(fasta, buffer_size) {
        match record_part {
            fasta_records_iterator::FastaPartType::ID(id) => {
                if let Some(prev_id) = prev_record_id.as_ref() {
                    if let Some(mut temp_fasta_stream) = temp_fasta_stream.take() {
                        temp_fasta_stream.flush().expect("E: Failed to flush temporary fasta file");
                        drop(temp_fasta_stream);

                        results.push((sequence_id2str(prev_id.as_slice()).to_string(), create_normalized_profile(k, FastaNucltudiesIterator::new(&temp_fasta_path, buffer_size), &false).1));
                    }
                }

                temp_fasta_stream = {
                    let seq_file_file = File::create(temp_fasta_path.clone()).unwrap_or_else(|_| panic!("ERROR: failed to create {}, a file needed foe meta kmer, quitting", temp_fasta_path.display()));
                    Some(BufWriter::new(seq_file_file))
                };

                let contig = sequence_id2str(id.as_slice());
                writeln!(temp_fasta_stream.as_mut().expect("E: Temporary fasta creation failed"), ">{}", contig).expect("E: Failed to write sequence id into temporary fasta file");
                prev_record_id = Some(id);
            },
            fasta_records_iterator::FastaPartType::Nuc(nuc) => {
                match temp_fasta_stream.as_mut() {
                    Some(temp_fasta_stream) => write!(temp_fasta_stream, "{}", nuc as char).expect("E: Failed to write into temporary file"),
                    None => panic!("E: The fasta file '{}' is malformed! a non-description line shows up before description line", fasta.display()),
                }
            }
        }
    }

    if let Some(prev_id) = prev_record_id.as_ref() {
        if let Some(mut temp_fasta_stream) = temp_fasta_stream.take() {
            temp_fasta_stream.flush().expect("E: Failed to flush temporary fasta file");
            drop(temp_fasta_stream);

            results.push((sequence_id2str(prev_id.as_slice()).to_string(), create_normalized_profile(k, FastaNucltudiesIterator::new(&temp_fasta_path, buffer_size), &false).1));
        }
    }

    results
}

pub fn print_kmers(input: &Path, output: &Path, k: usize, ratio: bool, buffer_size: usize, meta: bool) -> Result<(), PrintKmersError> {
    let mut output_stream = {
        let fout = File::create(output).unwrap_or_else(|_| panic!("E: Cannot create output file '{}'", output.display()));
        BufWriter::new(fout)
    };
    write!(output_stream, "Genome")?;
    for kem in create_printed_vector(k, &ALPHABET_MAP, &REV_ALPHABET_VALUES) {
        write!(output_stream, "\t{kem}")?;
    }
    writeln!(output_stream)?;

    let kmers = if meta {
        calc_kmers_meta(input, k, buffer_size)
    } else {
        calc_kmers(input, k, buffer_size)?
    };
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