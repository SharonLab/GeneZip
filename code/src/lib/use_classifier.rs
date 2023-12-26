use std::io::{BufWriter, Write};
use std::path::Path;

use crate::classifier::Classifier;
use crate::lz78::LenBases;
use crate::{fasta_records_iterator, samples_file_reader};
use crate::contig_naming::{are_genes_of_same_contig, get_contig_name, sequence_id2str};
use crate::fasta_nucleutide_iterator::FastaNucltudiesIterator;
use crate::fasta_record::FastaRecord;
use crate::output_streams::OutputStreams;
use crate::samples_file_reader::SampleError;
use crate::logger::log_event;

pub fn create_lz_classifier(mut log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                            max_depth: usize,
                            name2file: &Path,
                            buffer_size: usize,
                            kmer_size: &Option<usize>) -> Classifier {
    log_event(&mut log_stream, "Starting classifier creation");
    let len_bases: LenBases = LenBases::new(max_depth);
    let mut classifier: Classifier = Classifier::new(len_bases);

    log_event(&mut log_stream, "Training");
    if let Err(e) = classifier.batch_add_model(name2file, max_depth, buffer_size, kmer_size) {
        eprintln!("{}", e);
    }


    if let Some(&mut ref mut log_stream) = log_stream {
        classifier.print_stats(log_stream).expect("E: Failed to write into log");
    }
    log_event(&mut log_stream, "GeneZip stats are ready");

    classifier
}


fn write_classifier_prediction(classifier: &Classifier, sample_name: &str,
                               output_streams: &mut OutputStreams,
                               model_name2score: &(Vec<(&String, Option<f64>)>, usize)) -> Result<(), std::io::Error> {
    classifier.print_prediction(sample_name, output_streams, model_name2score)
        .map_err(|_| std::io::Error::new(std::io::ErrorKind::InvalidData,
                                         format!("ERROR: Failed to write prediction into {}", output_streams)))
}
pub fn predict_using_lz_classifier(mut log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                               buffer_size: usize,
                               kmer_size: &Option<usize>,
                               gc_limit: Option<f64>,
                               classifier: &Classifier,
                               prediction_name2file: &Path,
                               output_streams: &mut OutputStreams,
                               reflect: bool) -> Result<(), SampleError> {
    // Open the output stream
    log_event(&mut log_stream, "Predicting");

    classifier.print_header(output_streams).unwrap_or_else(|_| panic!("E: Failed to write header into output file '{}'", output_streams));
    for sample in samples_file_reader::SampleSource::new(prediction_name2file, false).into_iter() {
        let sample = sample?;
        let model_name2score = classifier.predict(FastaNucltudiesIterator::new(sample.get_path(), buffer_size), gc_limit, kmer_size, reflect);
        if let Err(e) = write_classifier_prediction(classifier, sample.get_name(), output_streams, &model_name2score) {
            panic!("{}", e);
        }
    }

    output_streams.flush().unwrap_or_else(|_| panic!("E: Failed to flush output stream into '{}'", output_streams));

    log_event(&mut log_stream, "GeneZip prediction is ready");

    Ok(())
}


fn meta_predict_using_lz_classifier_helper(temp_fasta_stream: &mut Option<FastaRecord>,
                                           genes: bool,
                                           min_genes: usize,
                                           found_genes: usize,
                                           classifier: &Classifier,
                                           gc_limit: Option<f64>,
                                           prev_id: &Vec<u8>,
                                           output_streams: &mut OutputStreams) {
    if let Some(temp_fasta_stream) = temp_fasta_stream.take() {
        if !genes || min_genes == 0 || found_genes >= min_genes {
            let model_name2score = classifier.predict(temp_fasta_stream, gc_limit, &None, false);
            let contig = sequence_id2str(if genes {
                get_contig_name(prev_id.as_slice())
            } else {
                prev_id.as_slice()
            });
            if let Err(e) = write_classifier_prediction(classifier, contig, output_streams, &model_name2score) {
                panic!("{}", e);
            }
        }
    }
}
pub fn meta_predict_using_lz_classifier(mut log_stream: Option<&mut BufWriter<Box<dyn Write>>>,
                                        buffer_size: usize,
                                        classifier: &Classifier,
                                        fasta: &Path,
                                        output_streams: &mut OutputStreams,
                                        genes: bool,
                                        min_genes: usize,
                                        gc_limit: Option<f64>) -> Result<(), SampleError> {
    log_event(&mut log_stream, "Predicting");

    classifier.print_header(output_streams).unwrap_or_else(|_| panic!("E: Failed to write header into output file '{}'", output_streams));

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
                            meta_predict_using_lz_classifier_helper(&mut temp_fasta_stream, genes, min_genes, found_genes, classifier, gc_limit, prev_id, output_streams);

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
        meta_predict_using_lz_classifier_helper(&mut temp_fasta_stream, genes, min_genes, found_genes, classifier, gc_limit, &prev_id, output_streams);
    }

    output_streams.flush().unwrap_or_else(|_| panic!("E: Failed to flush output stream into '{}'", output_streams));

    log_event(&mut log_stream, "GeneZip prediction is ready");

    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::BufRead;
    use std::path::PathBuf;
    use crate::output_streams::{OutputFileType, OutputStreams};
    use crate::use_classifier::{create_lz_classifier, meta_predict_using_lz_classifier, predict_using_lz_classifier};

    #[test]
    fn test_meta_genes_no_gz_output() {
        let basic_output_path = PathBuf::from("../tests/meta_small_sample_predication_basic_nogz.tsv");
        let mut output_streams = match OutputStreams::new(&[(OutputFileType::BaseGz, basic_output_path.as_path())].into_iter().collect()) {
            Err(e) => panic!("ERROR: failed to create output_streams, got {}", e),
            Ok(os) => os,
        };

        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        meta_predict_using_lz_classifier(None,
                                         512,
                                         &classifier,
                                         &PathBuf::from("../data/small_sample_as_meta.fna"),
                                         &mut output_streams,
                                         true,
                                         0,
                                         None).unwrap();
        // basic_output_path
        let lines = ["Genome_name\tLength\tBest_hit",
            "4_1\t2728891\t4",
            "4_2\t2569875\t4",
            "4_3\t2661614\t4",
            "8_1\t2213305\t8",
            "8_2\t2185249\t8",
            "8_3\t2222558\t8"];

        for (cl, kl) in std::fs::read(&basic_output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        std::fs::remove_file(&basic_output_path).unwrap();
    }
    #[test]
    fn test_meta_genes() {
        let basic_output_path = PathBuf::from("../tests/meta_small_sample_predication_basic.tsv");
        let lz_matrix_path = PathBuf::from("../tests/meta_small_sample_predication_gz.tsv");
        let mut output_streams = match OutputStreams::new(&[(OutputFileType::BaseGz, basic_output_path.as_path()),
                                                                         (OutputFileType::LzValues, lz_matrix_path.as_path())].into_iter().collect()) {
            Err(e) => panic!("ERROR: failed to create output_streams, got {}", e),
            Ok(os) => os,
        };

        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        meta_predict_using_lz_classifier(None,
                                         512,
                                         &classifier,
                                         &PathBuf::from("../data/small_sample_as_meta.fna"),
                                         &mut output_streams,
                                         true,
                                         0,
                                         None).unwrap();
        // basic_output_path
        let lines = ["Genome_name\tLength\tBest_hit",
            "4_1\t2728891\t4",
            "4_2\t2569875\t4",
            "4_3\t2661614\t4",
            "8_1\t2213305\t8",
            "8_2\t2185249\t8",
            "8_3\t2222558\t8"];

        for (cl, kl) in std::fs::read(&basic_output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        // gz_output_path
        let lines = ["Genome_name\t18\t20\t24\t4\t8",
            "4_1\t1.98068\t1.97345\t1.98103\t1.96102\t1.99144",
            "4_2\t1.97270\t1.96790\t1.97642\t1.95492\t1.98538",
            "4_3\t1.96976\t1.96571\t1.97475\t1.95316\t1.98320",
            "8_1\t1.95204\t1.97514\t1.97050\t1.96392\t1.93566",
            "8_2\t1.95018\t1.97452\t1.97047\t1.96315\t1.93396",
            "8_3\t1.95041\t1.97463\t1.96904\t1.96244\t1.93397"];
        for (cl, kl) in std::fs::read(&lz_matrix_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        std::fs::remove_file(&basic_output_path).unwrap();
        std::fs::remove_file(&lz_matrix_path).unwrap();
    }

    #[test]
    fn test_meta_genes_gc() {
        let basic_output_path = PathBuf::from("../tests/meta_small_sample_predication_gc_basic.tsv");
        let lz_matrix_path = PathBuf::from("../tests/meta_small_sample_predication_gc_gz.tsv");
        let mut output_streams = match OutputStreams::new(&[(OutputFileType::BaseGz, basic_output_path.as_path()),
                                                                         (OutputFileType::LzValues, lz_matrix_path.as_path())].into_iter().collect()) {
            Err(e) => panic!("ERROR: failed to create output_streams, got {}", e),
            Ok(os) => os,
        };

        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        meta_predict_using_lz_classifier(None,
                                         512,
                                         &classifier,
                                         &PathBuf::from("../data/small_sample_as_meta.fna"),
                                         &mut output_streams,
                                         true,
                                         0,
                                         Some(2.0)).unwrap();

        // Basic
        let lines = ["Genome_name\tLength\tBest_hit",
            "4_1\t2728891\t4",
            "4_2\t2569875\t4",
            "4_3\t2661614\t4",
            "8_1\t2213305\t8",
            "8_2\t2185249\t8",
            "8_3\t2222558\t8"];
        for (cl, kl) in std::fs::read(&basic_output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        // GZ
        let lines = ["Genome_name\t18\t20\t24\t4\t8",
            "4_1\tNA\t1.97345\t1.98103\t1.96102\tNA",
            "4_2\tNA\t1.96790\t1.97642\t1.95492\tNA",
            "4_3\tNA\t1.96571\t1.97475\t1.95316\t1.98320",
            "8_1\t1.95204\tNA\tNA\tNA\t1.93566",
            "8_2\t1.95018\tNA\tNA\tNA\t1.93396",
            "8_3\t1.95041\tNA\tNA\tNA\t1.93397"];
        for (cl, kl) in std::fs::read(&lz_matrix_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        std::fs::remove_file(&basic_output_path).unwrap();
        std::fs::remove_file(&lz_matrix_path).unwrap();
    }

    #[test]
    fn test_small_example() {
        let basic_output_path = PathBuf::from("../tests/small_sample_predication_basic.tsv");
        let lz_matrix_path = PathBuf::from("../tests/small_sample_predication_gz.tsv");

        let mut output_streams = match OutputStreams::new(&[(OutputFileType::BaseGz, basic_output_path.as_path()),
                                                                         (OutputFileType::LzValues, lz_matrix_path.as_path())].into_iter().collect()) {
            Err(e) => panic!("ERROR: failed to create output_streams, got {}", e),
            Ok(os) => os,
        };

        // max_depth is 12, for consistency with the small sample.
        let classifier = create_lz_classifier(None, 12, &PathBuf::from("../tests/small_example_training.txt"), 512, &None);
        predict_using_lz_classifier(None,
                                    512,
                                    &None,
                                    None,
                                    &classifier,
                                    &PathBuf::from("../tests/small_example_testing.txt"),
                                    &mut output_streams,
                                    false).unwrap();
        // Basic
        let lines = ["Genome_name\tLength\tBest_hit",
            "4\t2728891\t4",
            "4\t2569875\t4",
            "4\t2661614\t4",
            "8\t2213305\t8",
            "8\t2185249\t8",
            "8\t2222558\t8"];
        for (cl, kl) in std::fs::read(&basic_output_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }

        // GZ
        let lines = ["Genome_name\t18\t20\t24\t4\t8",
            "4\t1.98068\t1.97345\t1.98103\t1.96102\t1.99144",
            "4\t1.97270\t1.96790\t1.97642\t1.95492\t1.98538",
            "4\t1.96976\t1.96571\t1.97475\t1.95316\t1.98320",
            "8\t1.95204\t1.97514\t1.97050\t1.96392\t1.93566",
            "8\t1.95018\t1.97452\t1.97047\t1.96315\t1.93396",
            "8\t1.95041\t1.97463\t1.96904\t1.96244\t1.93397"];
        for (cl, kl) in std::fs::read(&lz_matrix_path).unwrap().lines().zip(lines) {
            assert_eq!(cl.unwrap(), kl);
        }


        std::fs::remove_file(&basic_output_path).unwrap();
        std::fs::remove_file(&lz_matrix_path).unwrap();
    }
}