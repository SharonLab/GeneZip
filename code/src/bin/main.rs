//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 2023/08/31



use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use chrono::Utc;

extern crate GeneZipLib;

use GeneZipLib::cli::{Usage, UserTask};
use GeneZipLib::ani_rating::create_fastani_run;
use GeneZipLib::kmer_database::KmerDatabase;
use GeneZipLib::kmer_prediction::KmerClassifier;
use GeneZipLib::print_kmer::print_kmers;
use GeneZipLib::samples_file_reader::{Sample, SampleError, SampleSource};
use GeneZipLib::database;
use GeneZipLib::use_classifier::{create_lz_classifier, meta_predict_using_lz_classifier, predict_using_lz_classifier};



fn run_ani<I>(log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
           ani_path: &Path,
           gz_output_file: &Path,
           prediction_name2file: &Path,
           training_name2file: I) -> Result<(), SampleError> where
    I: IntoIterator<Item=Result<Sample, SampleError>> {
    if let Some(log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tStarting fastANI step", now.to_rfc2822()).expect("E: Failed to write log");
    }

    create_fastani_run(gz_output_file,
                       ani_path,
                       prediction_name2file,
                       training_name2file)
}

fn kmer_prediction(log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                   buffer_size: usize,
                   kmer_size: usize,
                   prediction_name2file: &Path,
                   output_file: &Path,
                   name2file: &Path) -> Result<(), SampleError> {
    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tStarting kmer classifier creation", now.to_rfc2822()).expect("E: Failed to write log");
    }

    let mut classifier = KmerClassifier::new(kmer_size);

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tTraining", now.to_rfc2822()).expect("E: Failed to write log");
    }

    if let Err(e) = classifier.batch_add_model(name2file, buffer_size) {
        eprintln!("{}", e);
    }

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tKmer Prediction is now ready", now.to_rfc2822()).expect("E: Failed to write log");
    }

    //////// Prediction

    // Open the output stream
    let mut output_stream = {
        let fout = File::create(output_file).unwrap_or_else(|_| panic!("E: Cannot create output file '{}'", output_file.display()));
        BufWriter::new(fout)
    };

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tPredicting", now.to_rfc2822()).expect("E: Failed to write log");
    }

    classifier.print_header(&mut output_stream).unwrap_or_else(|_| panic!("E: Failed to write header into output file '{}'", output_file.display()));

    for sample in SampleSource::new(prediction_name2file, false) {
        let sample = sample?;
        let model_name2score = classifier.predict(sample.get_path(), buffer_size);
        classifier.print_prediction(sample.get_name(), &mut output_stream, &model_name2score).unwrap_or_else(|_| panic!("E: Failed to write prediction into '{}", output_file.display()));
    }

    output_stream.flush().unwrap_or_else(|_| panic!("E: Failed to flush output stream into '{}'", output_file.display()));

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tGeneZip prediction is ready", now.to_rfc2822()).expect("E: Failed to write log");
    }

    Ok(())
}

fn is_file_missing(file_path: &Path) -> bool {
    if file_path.exists() {
        false
    } else {
        eprintln!("E: The file '{}' is missing", file_path.display());
        true
    }
}

fn compute_task(usage: &Usage) {
    let mut log_stream =  if usage.get_print_statistics() {
        Some(BufWriter::new(Box::new(std::io::stderr()) as Box<dyn Write>))
    } else {
        None
    };

    if let Some(ref mut log_stream) = log_stream {
        writeln!(log_stream, "GeneZip, {}", usage.get_version()).expect("E: Failed to write log");

        let now = Utc::now();
        writeln!(log_stream, "{}\tStarting", now.to_rfc2822()).expect("E: Failed to write log");
    }

    match usage.get_task() {
        UserTask::BuildDB => {
            let training_name2file = usage.get_training_name2file_file().expect("E: Trying to read path to training file provided by the user, however, the user did not provide that. This should never happen");

            if ! is_file_missing(training_name2file) {
                let md = usage.get_max_depth().expect("E: Trying to use user-provided max depth, however, the user did not provide max depth. This should never happen.");
                let database_path = usage.get_database_path().expect("E: Trying to use user-provided database path, however, the user did not provde that information. This should never happen.");

                let classifier = create_lz_classifier(log_stream.as_mut(),
                                                      md,
                                                      training_name2file,
                                                      usage.get_buffer_size(),
                                                      &usage.get_kmer_size());

                let database = database::Database::new(classifier, md, usage.get_kmer_size());
                if let Err(e) = database.save(database_path, usage.get_jobs().unwrap_or(0)) {
                    panic!("E: Failed to save GeneZip database to '{}', encountered the following error: '{}'", database_path.display(), e)
                }
            }
        },
        UserTask::DBPredict => {
            let database_path = usage.get_database_path().expect("E: Trying to use user-provided database path, however, the user did not provde that information. This should never happen.");
            let prediction_name2file = usage.get_prediction_name2file_file().expect("E: Trying to get prediction input path, however, the user was not asked to provide that. This should never happen.");
            if ! is_file_missing(database_path) && ! is_file_missing(prediction_name2file) {
                let database = match database::Database::load(database_path, usage.get_jobs().unwrap_or(0)) {
                    Ok(db) => db,
                    Err(e) => panic!("E: Failed to read GeneZip database from '{}', encountered the following error: '{}'", database_path.display(), e),
                };
                let gz_output_path = usage.get_out_file().expect("E: Trying to open the output file, but no path was provided by user. This should never happen.");
                if let Err(e) = predict_using_lz_classifier(log_stream.as_mut(),
                                            usage.get_buffer_size(),
                                            database.get_kmer_size(),
                                            usage.get_gc_limit(),
                                            database.get_classifier(),
                                            prediction_name2file,
                                                            gz_output_path,
                                            usage.get_reflect()) {
                    eprintln!("{}", e);
                } else if let Some(ani_path) = usage.get_ani_out_file() {
                    let samples_iterator: Box<dyn Iterator<Item=Result<Sample, SampleError>>> = if let Some(training_name2file) = usage.get_training_name2file_file() {
                        Box::new(SampleSource::new(training_name2file, false).into_iter())
                    } else {
                        Box::new(database.get_classifier().into_iter().map(Ok))
                    };
                    if let Err(e) = run_ani(log_stream.as_mut(),
                                            ani_path,
                                            gz_output_path,
                                            prediction_name2file,
                                            samples_iterator) {
                        eprintln!("{}", e);
                    }
                }
            }
        },
        UserTask::Predict => {
            let md = usage.get_max_depth().expect("E: Trying to use user-provided max depth, however, the user did not provide max depth. This should never happen.");
            let gz_output_path = usage.get_out_file().expect("E: Trying to get the output file, but no path was provided by user. This should never happen.");
            let prediction_name2file = usage.get_prediction_name2file_file().expect("E: Trying to get prediction input path, however, the user was not asked to provide that. This should never happen.");
            let training_name2file = usage.get_training_name2file_file().expect("E: Trying to read path to training file provided by the user, however, the user did not provide that. This should never happen");

            if ! is_file_missing(prediction_name2file) && ! is_file_missing(training_name2file) {
                let classifier = create_lz_classifier(log_stream.as_mut(),
                                                      md,
                                                      training_name2file,
                                                      usage.get_buffer_size(),
                                                      &usage.get_kmer_size());
                if let Err(e) = predict_using_lz_classifier(log_stream.as_mut(),
                                            usage.get_buffer_size(),
                                            &usage.get_kmer_size(),
                                            usage.get_gc_limit(),
                                            &classifier,
                                            prediction_name2file,
                                            gz_output_path,
                                            usage.get_reflect()) {
                    eprintln!("{}", e);
                }
                if let Some(ani_path) = usage.get_ani_out_file() {
                    if let Err(e) = run_ani(log_stream.as_mut(),
                            ani_path,
                            gz_output_path,
                            prediction_name2file,
                            SampleSource::new(training_name2file, false)) {
                        eprintln!("{}", e);
                    }
                }
            }
        },
        UserTask::PrintKmer => {
            let training_name2file_file = usage.get_training_name2file_file().expect("E: Trying to read path to sequences file provided by the user, however, the user did not provide that. This should never happen");
            if ! is_file_missing(training_name2file_file) {
                match print_kmers(training_name2file_file,
                                  usage.get_out_file().expect("E: Trying to get the output file, but no path was provided by user. This should never happen."),
                                  usage.get_kmer_size().expect("E: Trying to get k chosen by the user, however, no such parameter was taken. This should never happen."),
                                  usage.get_ratio(),
                                  usage.get_buffer_size(),
                                  usage.get_meta()) {
                    Ok(()) => (),
                    Err(e) => eprintln!("E: Failed to print kmers, got the following error: {:?}", e),
                }
            }
        },
        UserTask::BuildKmer => {
            let training_name2_file = usage.get_training_name2file_file().expect("E: Trying to read path to sequences file provided by the user, however, the user did not provide that. This should never happen");
            if ! is_file_missing(training_name2_file) {
                match KmerDatabase::new(training_name2_file,
                                        usage.get_buffer_size(),
                                        usage.get_kmer_size().expect("E: Trying to get k chosen by the user, however, no such parameter was taken. This should never happen.")) {
                    Ok(kmers) => if let Err(e) = kmers.save(usage.get_out_file().expect("E: Trying to get the output file, but no path was provided by user. This should never happen.")) {
                        eprintln!("E: Failed to save k-mers database due to the following error: '{}'", e)
                    },
                    Err(e) => eprintln!("E: {:}", e),
                }
            }
        },
        UserTask::KMerPredict => {
            let training_name2file = usage.get_training_name2file_file().expect("E: Trying to read path to sequences file provided by the user, however, the user did not provide that. This should never happen");
            let prediction_name2file = usage.get_prediction_name2file_file().expect("E: Trying to get prediction input path, however, the user was not asked to provide that. This should never happen.");
            if ! is_file_missing(training_name2file) && ! is_file_missing(prediction_name2file) {
                let kmer = usage.get_kmer_size().expect("E: Trying to get k chosen by the user, however, no such parameter was taken. This should never happen.");
                if kmer < 1 {
                    eprintln!("E: Can not create kmer frequencies prediction with k=0");
                } else if let Err(e) = kmer_prediction(log_stream.as_mut(),
                                    usage.get_buffer_size(),
                                    kmer,
                                    training_name2file,
                                    usage.get_out_file().expect("E: Trying to get the output file, but no path was provided by user. This should never happen."),
                                    prediction_name2file) {
                        eprintln!("{}", e);
                    }
            }
        },
        UserTask::MetaPredict => {
            let training_name2file = usage.get_training_name2file_file().expect("E: Trying to read path to sequences file provided by the user, however, the user did not provide that. This should never happen");
            let prediction_name2file = usage.get_prediction_name2file_file().expect("E: Trying to get prediction input path, however, the user was not asked to provide that. This should never happen.");
            let gz_output_path = usage.get_out_file().expect("E: Trying to get the output file, but no path was provided by user. This should never happen.");
            if ! is_file_missing(training_name2file) && ! is_file_missing(prediction_name2file) {
                let md = usage.get_max_depth().expect("E: Trying to use user-provided max depth, however, the user did not provide max depth. This should never happen.");
                let classifier = create_lz_classifier(log_stream.as_mut(),
                                                      md,
                                                      training_name2file,
                                                      usage.get_buffer_size(),
                                                      &usage.get_kmer_size());
                if let Err(e) = meta_predict_using_lz_classifier(log_stream.as_mut(),
                                                                            usage.get_buffer_size(),
                                                                            &classifier,
                                                                            prediction_name2file,
                                                                            gz_output_path,
                                                                            usage.get_genes(),
                                                                            usage.get_min_genes(),
                                                                            usage.get_gc_limit()) {
                    eprintln!("{}", e);
                }
            }
        },
    }

    if let Some(mut log_stream) = log_stream {
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

    compute_task(&usage)
}