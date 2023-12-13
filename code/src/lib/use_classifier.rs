use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use chrono::Utc;
use crate::classifier::Classifier;
use crate::lz78::LenBases;
use crate::{fasta_records_iterator, samples_file_reader};
use crate::contig_naming::{are_genes_of_same_contig, get_contig_name, sequence_id2str};
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;
use crate::fasta_record::FastaRecord;
use crate::samples_file_reader::SampleError;

pub fn create_lz_classifier(log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                            max_depth: usize,
                            name2file: &Path,
                            buffer_size: usize,
                            kmer_size: &Option<usize>) -> Classifier {
    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tStarting classifier creation", now.to_rfc2822()).expect("E: Failed to write log");
    }
    let len_bases: LenBases = LenBases::new(max_depth);
    let mut classifier: Classifier = Classifier::new(len_bases);

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tTraining", now.to_rfc2822()).expect("E: Failed to write log");
    }

    if let Err(e) = classifier.batch_add_model(name2file, max_depth, buffer_size, kmer_size) {
        eprintln!("{}", e);
    }

    if let Some(log_stream) = log_stream {
        classifier.print_stats(log_stream).expect("E: Failed to write into log");
        let now = Utc::now();
        writeln!(log_stream, "{}\tGeneZip stats are ready", now.to_rfc2822()).expect("E: Failed to write log");
    }

    classifier
}

pub fn predict_using_lz_classifier(log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                               buffer_size: usize,
                               kmer_size: &Option<usize>,
                               gc_limit: Option<f64>,
                               classifier: &Classifier,
                               prediction_name2file: &Path,
                               output_file: &Path,
                               reflect: bool) -> Result<(), SampleError> {
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
    for sample in samples_file_reader::SampleSource::new(prediction_name2file, false).into_iter() {
        let sample = sample?;
        let model_name2score = classifier.predict(FastaNucltudiesIterator::new(sample.get_path(), buffer_size), gc_limit, kmer_size, reflect);
        classifier.print_prediction(sample.get_name(), &mut output_stream, &model_name2score).unwrap_or_else(|_| panic!("E: Failed to write prediction into '{}", output_file.display()));
    }

    output_stream.flush().unwrap_or_else(|_| panic!("E: Failed to flush output stream into '{}'", output_file.display()));

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tGeneZip prediction is ready", now.to_rfc2822()).expect("E: Failed to write log");
    }

    Ok(())
}

pub fn meta_predict_using_lz_classifier(log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                                    buffer_size: usize,
                                    classifier: &Classifier,
                                    fasta: &Path,
                                    output_file: &Path,
                                    genes: bool,
                                    min_genes: usize) -> Result<(), SampleError> {
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

    let mut temp_fasta_stream = None;
    let mut prev_record_id = None;
    let mut found_genes = 0_usize;
    for record_part in fasta_records_iterator::FastaRecordIterator::new(fasta, buffer_size) {
        match record_part {
            fasta_records_iterator::FastaPartType::ID(id) => {
                match prev_record_id.as_ref() {
                    None => { // First sequence
                        let contig = sequence_id2str(id.as_slice());
                        prev_record_id = Some(id.clone());
                        found_genes += 1;
                        temp_fasta_stream = Some(FastaRecord::new(contig, "", buffer_size));
                    },
                    Some(prev_id) => {  // Not the first sequence
                        // replace is true if we need to print and replace the sequence we follow
                        let replace = if genes {
                            if are_genes_of_same_contig(prev_id, &id) {
                                // We write it here to avoid duplicating "id" as the compiler fails to prove it will not be moved to a "else" block.
                                found_genes += 1;
                                if let Some(ref mut f) = temp_fasta_stream.as_mut() {
                                    f.push(b'N'); // This will create a separation between genes of the same genome
                                }
                                None // We continue with this gene
                            } else {
                                Some(id)
                            }
                        } else {
                            Some(id)
                        };
                        if let Some(replace) = replace {
                            if let Some(temp_fasta_stream) = temp_fasta_stream.take() {
                                if !genes || min_genes == 0 || found_genes >= min_genes {
                                    let contig = sequence_id2str(if genes {
                                        get_contig_name(prev_id.as_slice())
                                    } else {
                                        prev_id.as_slice()
                                    });
                                    let model_name2score = classifier.predict(temp_fasta_stream, None, &None, false);
                                    classifier.print_prediction(contig,
                                                                &mut output_stream, &model_name2score).unwrap_or_else(|_| panic!("E: Failed to write prediction into '{}", output_file.display()));
                                }
                            }

                            {
                                found_genes = 1;
                                let contig = sequence_id2str(replace.as_slice());
                                temp_fasta_stream = Some(FastaRecord::new(contig, "", buffer_size));
                            }

                            let _ = prev_record_id.insert(replace.clone());
                        }
                    }
                }
            },
            fasta_records_iterator::FastaPartType::Nuc(nuc) => {
                match temp_fasta_stream.as_mut() {
                    Some(temp_fasta_stream) => temp_fasta_stream.push(nuc),
                    None => panic!("E: The fasta file '{}' is malformed! a non-description line shows up before description line", fasta.display()),
                }
            }
        }
    }

    if let Some(prev_id) = prev_record_id {
        if let Some(temp_fasta_stream) = temp_fasta_stream.take() {

            if !genes || min_genes == 0 || found_genes >= min_genes {
                let model_name2score = classifier.predict(temp_fasta_stream, None, &None, false);
                let contig = sequence_id2str(if genes {
                    get_contig_name(prev_id.as_slice())
                } else {
                    prev_id.as_slice()
                });
                classifier.print_prediction(contig,
                                            &mut output_stream, &model_name2score).unwrap_or_else(|_| panic!("E: Failed to write prediction into '{}", output_file.display()));
            }
        }
    }

    output_stream.flush().unwrap_or_else(|_| panic!("E: Failed to flush output stream into '{}'", output_file.display()));

    if let Some(&mut ref mut log_stream) = log_stream {
        let now = Utc::now();
        writeln!(log_stream, "{}\tGeneZip prediction is ready", now.to_rfc2822()).expect("E: Failed to write log");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::BufRead;
    use std::path::PathBuf;
    use crate::use_classifier::{create_lz_classifier, meta_predict_using_lz_classifier, predict_using_lz_classifier};

    #[test]
    fn test_meta_genes() {
        let output_path = PathBuf::from("../tests/meta_small_sample_predication.tsv");
        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        meta_predict_using_lz_classifier(None,
                                         512,
                                         &classifier,
                                         &PathBuf::from("../data/small_sample_as_meta.fna"),
                                         &output_path,
                                         true,
                                         0).unwrap();
        let lines = ["Genome_name\t18\t20\t24\t4\t8\tBest_hit",
            "4_1\t1.98068\t1.97345\t1.98103\t1.96102\t1.99144\t4",
            "4_2\t1.97270\t1.96790\t1.97642\t1.95492\t1.98538\t4",
            "4_3\t1.96976\t1.96571\t1.97475\t1.95316\t1.98320\t4",
            "8_1\t1.95204\t1.97514\t1.97050\t1.96392\t1.93566\t8",
            "8_2\t1.95018\t1.97452\t1.97047\t1.96315\t1.93396\t8",
            "8_3\t1.95041\t1.97463\t1.96904\t1.96244\t1.93397\t8"];
        for (cl, kl) in std::fs::read(&output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        std::fs::remove_file(&output_path).unwrap();
    }

    #[test]
    fn test_small_example() {
        let output_path = PathBuf::from("../tests/small_sample_predication.tsv");
        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        predict_using_lz_classifier(None,
                                    512,
                                    &None,
                                    None,
                                    &classifier,
                                    &PathBuf::from("../tests/small_example_testing.txt"),
                                    &output_path,
                                    false).unwrap();
        let lines = ["Genome_name\t18\t20\t24\t4\t8\tBest_hit",
            "4\t1.98068\t1.97345\t1.98103\t1.96102\t1.99144\t4",
            "4\t1.97270\t1.96790\t1.97642\t1.95492\t1.98538\t4",
            "4\t1.96976\t1.96571\t1.97475\t1.95316\t1.98320\t4",
            "8\t1.95204\t1.97514\t1.97050\t1.96392\t1.93566\t8",
            "8\t1.95018\t1.97452\t1.97047\t1.96315\t1.93396\t8",
            "8\t1.95041\t1.97463\t1.96904\t1.96244\t1.93397\t8"];
        for (cl, kl) in std::fs::read(&output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        std::fs::remove_file(&output_path).unwrap();
    }
}