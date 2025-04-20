//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31

// ANI module for the main GeneZip tool.

use std::fmt::Formatter;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use hashbrown::HashMap;
use tempdir::TempDir;
use std::process::{Command, ExitStatus};
use rayon::prelude::*;
use crate::ani_calculator_tool::AniCalculatorTool;
use crate::samples_file_reader;
use crate::samples_file_reader::{Sample, SampleError, SampleErrorType};

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

fn run_ani_command(reference_genome: &str, threads_limit: usize, ani_calculator_tool: AniCalculatorTool, query_path: &Path, output_path: &Path) -> std::io::Result<ExitStatus> {
    // TODO: not the best looking solutions, make it better
    match ani_calculator_tool {
        AniCalculatorTool::FastANI => Command::new("fastANI")
            .arg("-r")
            .arg(reference_genome)
            .arg("--ql")
            .arg(query_path.display().to_string())
            .arg("-o")
            .arg(output_path.display().to_string())
            .arg("-t")
            .arg(threads_limit.to_string())
            .spawn()?.wait(),
            // .status(),c
        AniCalculatorTool::Skani => Command::new("skani")
            .arg("dist")
            .arg("-r")
            .arg(reference_genome)
            .arg("--ql")
            .arg(query_path.display().to_string())
            .arg("-o")
            .arg(output_path.display().to_string())
            .arg("-t")
            .arg(threads_limit.to_string())
            .spawn()?.wait()
            // .status(),
    }
    /*
    match ani_calculator_tool {
        AniCalculatorTool::FastANI => &mut Command::new("fastANI"),
        AniCalculatorTool::Skani => Command::new("skani").arg("dist"),
    }.arg("-r")
        .arg(reference_genome)
        .arg("--ql")
        .arg(query_path.display().to_string())
        .arg("-o")
        .arg(output_path.display().to_string())
        .arg("-t")
        .arg(threads_limit.to_string())
        // .stderr(Stdio::null())
        // .stdout(Stdio::null())
        .status()
     */
}

enum ANIRunError {
    IOError(std::io::Error),
    ExitStatusError(String),
    CommandError(String),
    ANIFormat(String)
}

impl From<std::io::Error> for ANIRunError {
    fn from(value: std::io::Error) -> Self {
        ANIRunError::IOError(value)
    }
}

impl std::fmt::Display for ANIRunError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ANIRunError::IOError(e) => e.fmt(f),
            ANIRunError::ExitStatusError(e) => e.fmt(f),
            ANIRunError::CommandError(e) => e.fmt(f),
            ANIRunError::ANIFormat(e) => e.fmt(f),
        }
    }
}


fn get_genome_name<'a>(ani_calculator_tool: AniCalculatorTool, line_parts: &mut dyn Iterator<Item=&'a str>) -> Result<&'a str, ANIRunError> {
    match ani_calculator_tool {
        AniCalculatorTool::FastANI => {
            let temp = match line_parts.next() {
                Some(name) => Ok(name),
                None => return Err(ANIRunError::ANIFormat("ERROR: ANI results file has a line that is not empty but has no genome name".to_string())),
            };
            line_parts.next();  // Skip reference genome
            temp
        },
        AniCalculatorTool::Skani => {
            line_parts.next();  // Skip reference genome
            match line_parts.next() {
                Some(name) => Ok(name),
                None => return Err(ANIRunError::ANIFormat("ERROR: ANI results file has a line that is not empty but has no genome name".to_string())),
            }
        }
    }
}

fn collect_ani_results(ani_calculator_tool: AniCalculatorTool, results: &mut HashMap<String, f64>, output_path: &Path) -> Result<(), ANIRunError> {
    let f = File::open(output_path)?;
    let buf_reader = BufReader::new(f);
    let mut skip_lines = match ani_calculator_tool {
        AniCalculatorTool::FastANI => 0,
        AniCalculatorTool::Skani => 1,
    };

    for line in buf_reader.lines() {
        let line = line?;
        if ! line.is_empty() {
            if skip_lines > 0 { // skani has a header line we need to skip
                skip_lines -= 1;
                continue;
            }

            let mut line_parts = line.split('\t');
            let genome_name = get_genome_name(ani_calculator_tool, &mut line_parts)?;

            let ani_value = match line_parts.next() {
                Some(value) => match value.parse::<f64>() {
                    Ok(value) => value,
                    Err(_) => return Err(ANIRunError::ANIFormat("ERROR: ANI results file has a line with ani value that is not f64".to_string())),
                },
                None => return Err(ANIRunError::ANIFormat("ERROR: ANI results file has a line that is not empty but has no ani value".to_string())),
            };
            results.insert(genome_name.to_string(), ani_value);
        }
    }
    
    Ok(())
}

fn run_ani_collection<I>(reference_genome: &str, query_genomes: I, threads_limit: usize, ani_calculator_tool: AniCalculatorTool) -> Result<HashMap<String, f64>, ANIRunError>
where I: IntoIterator<Item=String> {
    let work_dir = TempDir::new("genezip")?;
    let query_path = work_dir.path().join("query.list");
    let output_path = work_dir.path().join("output.tsv");

    let mut results = HashMap::new();

    {   // Create a file with list of genomes for ANI input
        let f = File::create(query_path.clone())?;
        let mut buf_writer = BufWriter::new(f);
        for qg in query_genomes {
            results.insert(qg.to_string(), -1_f64);  // Assume no similarity
            writeln!(buf_writer, "{qg}")?;
        }
    }

    match run_ani_command(reference_genome, threads_limit, ani_calculator_tool, &query_path, &output_path) {
        Err(e) => return Err(ANIRunError::CommandError(format!("E: failed to run '{}', due to '{}'", ani_calculator_tool, e))),
        Ok(exit_status) => if ! exit_status.success() {
            return Err(ANIRunError::ExitStatusError(format!("E: failed to run '{}', exit status is '{}'", ani_calculator_tool, exit_status)));
        }
    }

    collect_ani_results(ani_calculator_tool, &mut results, &output_path)?;

    work_dir.close()?;
    Ok(results)
}

fn transform_path2ani_into_name2ani(path2ani: &HashMap<String, f64>, name2path: &HashMap<String, String>) -> HashMap<String, f64> {
    name2path
        .iter()
        .filter_map(|(name, path)| path2ani.get(path).map(|ani| (name.clone(), *ani)))
        .collect::<HashMap<String, f64>>()
}

pub fn create_ani_run<I>(genezip_output_table: &Path, output_path: &Path, testing_path: &Path, traning_path: I, ani_calculator_tool: AniCalculatorTool) -> Result<(), SampleError> where
    I: IntoIterator<Item=Result<Sample, SampleError>> {
    let reference2queries = get_reference2queries(genezip_output_table);
    let training_name2path = traning_path.into_iter()
        .map(|sample| sample.map (|sample| (sample.get_name().to_string(), sample.get_path().display().to_string())))
        .collect::<Result<HashMap<String, String>, SampleError>>()?;
    
    let testing_name2path = samples_file_reader::SampleSource::new(testing_path, false)
        .into_iter()
        .map(|sample| sample.map (|sample| (sample.get_name().to_string(), sample.get_path().display().to_string())))
        .collect::<Result<HashMap<String, String>, SampleError>>()?;

    let output_file = match File::create(output_path) {
        Ok(f) => Ok(f),
        Err(e) => Err(SampleError::new(&format!("Failed to create the output file {}, quitting", output_path.display()), SampleErrorType::IoError(e))),
    }?;
    
    let mut buf_output_stream = BufWriter::new(output_file);
    writeln!(buf_output_stream, "genome_name\treference\tANI").expect("ERROR: failed to write output to file");

    let name_ref_2_ani= reference2queries.into_iter()
        .par_bridge()
        .map(|(reference, query)| {
            let path2ani = run_ani_collection(training_name2path.get(&reference).expect("ERROR: A model name found with no path"),
                                              query.iter().map(|name| testing_name2path[name].clone()),
                                              1, ani_calculator_tool);

            let name2ani = match path2ani {
                Ok(path2ani) => Ok(transform_path2ani_into_name2ani(&path2ani, &testing_name2path)),
                Err(e) => Err(e),
            };
            (reference, name2ani)
        })
        .collect::<Vec<(String, Result<HashMap<String, f64>, ANIRunError>)>>();

        for (reference, name2ani) in name_ref_2_ani {
            match name2ani {
                Ok(name2ani) => for (name, ani) in name2ani {
                    writeln!(buf_output_stream, "{name}\t{reference}\t{ani}").expect("ERROR: failed to write output to file");
                },
                Err(e) => {
                    return Err(SampleError::new(&format!("E: running ani for '{}' failed because '{}'", reference, e), SampleErrorType::None))
                },
            }
        }

    Ok(())
}