extern crate core;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use chrono::Utc;

use crate::cli::Usage;
use crate::classifier::Classifier;
use crate::lz78::LenBases;

mod lz78;
mod classifier;
mod cli;

// <const MAX_DEPTH: usize, const BUFFER_SIZE: usize>
fn run_lz_classifier(usage: &Usage) {
    // Open the output stream
    let mut output_stream = {
        let fout = File::create(usage.get_out_file()).expect("E: Cannot create output file %s\n\n");
        BufWriter::new(fout)
    };

    let log_stream =  &mut BufWriter::new(Box::new(std::io::stderr()) as Box<dyn Write>);

    if usage.get_print_statistics() {
        writeln!(log_stream, "GeneZip, {}", usage.get_version()).expect("E: Failed to write log");

        let now = Utc::now();
        writeln!(log_stream, "{}\tStarting", now.to_rfc2822()).expect("E: Failed to write log");
    }

    let len_bases:LenBases = LenBases::new(usage.get_max_depth());
    let mut classifier: Classifier = Classifier::new(&len_bases);


    if usage.get_print_statistics() {
        let now = Utc::now();
        writeln!(log_stream, "{}\tTraining", now.to_rfc2822()).expect("E: Failed to write log");
    }

    classifier.batch_add_model(usage.get_training_name2file_file(), usage.get_max_depth(), usage.get_buffer_size());

    let predication_reader = {
        let f = File::open(usage.get_prediction_name2file_file())
            .unwrap_or_else(|_| panic!("E: failed to read file {}", usage.get_training_name2file_file().display()));
        BufReader::new(f)
    };

    if usage.get_print_statistics() {
        let now = Utc::now();
        writeln!(log_stream, "{}\tPredicting", now.to_rfc2822()).expect("E: Failed to write log");
    }

    classifier.print_header(&mut output_stream).expect("E: failed to write output");
    for line in predication_reader.lines() {
        let line = line.expect("E: failed to read line from prediction file");

        let mut split = line.split('\t');
        let name = split.next().unwrap_or_else(|| panic!("E: invalid line in file {}", usage.get_prediction_name2file_file().display()));
        let file_path = PathBuf::from(split.next().expect("E: invalid line"));
        let model_name2score = classifier.predict(file_path.as_path());
        classifier.print_prediction(name, &mut output_stream, &model_name2score).expect("E: Failed to write prediction");
    }

    if usage.get_print_statistics() {
        classifier.print_stats(log_stream).expect("E: failed to write into log");
        let now = Utc::now();
        writeln!(log_stream, "{}\tDone", now.to_rfc2822()).expect("E: Failed to write log");
    }
}

fn main() {
    if usize::BITS != 64 {
        eprintln!("W: developed and tested for 64 bit systems, you are using a {} system.", usize::BITS);
    }
    let usage = Usage::new();
    if let Some(jobs) = usage.get_jobs() {
        rayon::ThreadPoolBuilder::new().num_threads(jobs).build_global().unwrap();
    }
    run_lz_classifier(&usage);
}