//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use hashbrown::HashMap;
use tempdir::TempDir;
use std::process::{Command, Stdio};
use rayon::prelude::*;
use crate::samples_file_reader;
use crate::samples_file_reader::{Sample, SampleError};

fn get_reference2queries(gene_zip_best_hit_table: &Path) -> HashMap<String, Vec<String>> {
    let file = File::open(gene_zip_best_hit_table).expect("ERROR: failed to open GeneZip best hit table for reading, quitting");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Find the index of the Best_hit column. TODO: make the column name and separator into a parameter
    let first_line = match lines.next() {
        Some(Ok(line)) => line,
        Some(Err(e)) => panic!("ERROR: failed to read GeneZip best hit table, got the following error: {e:?}"),
        None => panic!("ERROR: GeneZip best hit table is empty, quitting"),
    };

    let Some((best_hit_column_index, _)) = first_line.split('\t').enumerate().find(|(_, column)| *column == "Best_hit") else {
        panic!("ERROR: GeneZip best table has no 'Best_hit' column, quitting")
    };

    let mut reference2queries = HashMap::new();
    for line in lines {
        let line = line.expect("ERROR: reading GeneZip best hit table failed in the middle, quitting!");
        if !line.is_empty() {
            let mut line_parts = line.split('\t');
            let genome_name = line_parts.next().expect("ERROR: no genome name found in GeneZip best hit table, quitting!");
            if let Some((_, best_hit)) =  line_parts.enumerate().find(|(index, _)| *index + 1 == best_hit_column_index) {
                let entry = reference2queries.entry_ref(&best_hit.to_string()).or_insert(Vec::new());
                entry.push(genome_name.to_string());
            } else {
                panic!("ERROR: GeneZip best hit table is uneven, quitting!")
            };
        }
    }

    reference2queries
}

fn run_fastani<I>(reference_genome: &str, query_genomes: I, threads_limit: usize) -> HashMap<String, f64>
where I: IntoIterator<Item=String> {
    let work_dir = TempDir::new("genezip").expect("ERROR: failed to create a temporary folder, your tmp may be full, quitting!");
    let query_path = work_dir.path().join("query.list");
    let output_path = work_dir.path().join("output.tsv");

    let mut results = HashMap::new();

    {   // Create a file with list of genomes for fastANI input
        let f = File::create(query_path.clone()).unwrap_or_else(|_| panic!("ERROR: failed to create {}, a file needed to run fastANI, quitting", query_path.display()));
        let mut buf_writer = BufWriter::new(f);
        for qg in query_genomes {
            results.insert(qg.to_string(), -1_f64);  // Assume no similarity
            writeln!(buf_writer, "{qg}").unwrap_or_else(|_| panic!("ERROR: failed to write into {}, a file needed to run fastANI, quitting", query_path.display()));
        }
    }

    let _ani_command = Command::new("fastANI")
        .arg("-r")
        .arg(reference_genome)
        .arg("--ql")
        .arg(query_path.display().to_string())
        .arg("-o")
        .arg(output_path.display().to_string())
        .arg("-t")
        .arg(format!("{threads_limit}"))
        .stderr(Stdio::null())
        .stdout(Stdio::null())
        .status()
        .expect("ERROR: failed to run fastANI");

    {
        let f = File::open(output_path.clone()).unwrap_or_else(|_| panic!("ERROR: failed to read {}, fastANIs output file", output_path.display()));
        let buf_reader = BufReader::new(f);
        for line in buf_reader.lines() {
            let line = line.expect("ERROR: failed to read a line from fastANI output, quitting");
            if ! line.is_empty() {
                let mut line_parts = line.split('\t');
                let genome_name = line_parts.next().expect("ERROR: fastANI results file has a line that is not empty but has no genome name");
                line_parts.next();  // Skip
                let ani_value = line_parts
                    .next()
                    .expect("ERROR: fastANI results file has a line that is not empty but has no ani value")
                    .parse::<f64>()
                    .expect("ERROR: fastANI results file has a line with ani value that is not f64");
                    results.insert(genome_name.to_string(), ani_value);
            }
        }
    }

    if let Err(e) = work_dir.close() {
        eprintln!("WARNING: tried to close and remove temporary folder, got an error: {e}");
    }

    results
}

fn transform_path2ani_into_name2ani(path2ani: &HashMap<String, f64>, name2path: &HashMap<String, String>) -> HashMap<String, f64> {
    name2path
        .iter()
        .filter_map(|(name, path)| path2ani.get(path).map(|ani| (name.clone(), *ani)))
        .collect::<HashMap<String, f64>>()
}

pub fn create_fastani_run<I>(genezip_output_table: &Path, output_path: &Path, testing_path: &Path, traning_path: I) -> Result<(), SampleError> where
    I: IntoIterator<Item=Result<Sample, SampleError>> {
    let reference2queries = get_reference2queries(genezip_output_table);
    let training_name2path = traning_path.into_iter()
        .map(|sample| sample.map (|sample| (sample.get_name().to_string(), sample.get_path().display().to_string())))
        .collect::<Result<HashMap<String, String>, SampleError>>()?;
    // let training_name2path = samples_file_reader::SampleSource::new(traning_path, false)
    //     .into_iter()
    //     .map(|sample| sample.map (|sample| (sample.get_name().to_string(), sample.get_path().display().to_string())))
    //     .collect::<Result<HashMap<String, String>, SampleError>>()?;
    let testing_name2path = samples_file_reader::SampleSource::new(testing_path, false)
        .into_iter()
        .map(|sample| sample.map (|sample| (sample.get_name().to_string(), sample.get_path().display().to_string())))
        .collect::<Result<HashMap<String, String>, SampleError>>()?;

    let output_file = File::create(output_path).unwrap_or_else(|_| panic!("ERROR: Failed to create the output file {}, quitting", output_path.display()));
    let mut buf_output_stream = BufWriter::new(output_file);
    writeln!(buf_output_stream, "genome_name\treference\tfastANI").expect("ERROR: failed to write output to file");

    reference2queries.into_iter()
        .par_bridge()
        .map(|(reference, query)| {
            let path2ani = run_fastani(training_name2path.get(&reference).expect("ERROR: A model name found with no path"),
                                       query.iter().map(|name| testing_name2path[name].clone()),
                                       1);

            let name2ani = transform_path2ani_into_name2ani(&path2ani, &testing_name2path);
            (reference, name2ani)
        })
        .collect::<Vec<(String, HashMap<String, f64>)>>()
        .iter()
        .for_each(|(reference, name2ani)| for (name, ani) in name2ani {
            writeln!(buf_output_stream, "{name}\t{reference}\t{ani}").expect("ERROR: failed to write output to file");
        });

    Ok(())
}