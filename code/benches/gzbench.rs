
use std::path::PathBuf;
use std::time::Duration;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use GeneZipLib::fasta_nucleutide_iterator::FastaNucltudiesIterator;

use GeneZipLib::lz78::{LenBases, LZ78};
use GeneZipLib::use_classifier::{create_lz_classifier, meta_predict_using_lz_classifier, predict_using_lz_classifier};

fn bench_presentation_example_two_seqs(max_depth: usize) {
    let train_path = PathBuf::from("../tests/presentation_train.fna");
    let test_path = PathBuf::from("../tests/presentation_test_two_seq.fna");

    let len_bases = LenBases::new(max_depth);
    let model = LZ78::new(max_depth, len_bases, FastaNucltudiesIterator::new(train_path.as_path(), 1024));
    let prediction = model.average_log_score(FastaNucltudiesIterator::new(test_path.as_path(), 1024));
    // - 1 / 8 * log2(1/22^4)
    assert!(prediction < 2.229715809318649);
    assert!(prediction > 2.229715809318647);
}

fn bench_small_example(max_depth: usize, buffer_size: usize) {
    let output_path = PathBuf::from("../tests/small_sample_predication.tsv");
    // max_depth is 12, for consistency with the small sample.
    let classifier = create_lz_classifier(None, max_depth, &PathBuf::from("../tests/small_example_training.txt"), buffer_size, &None);
    predict_using_lz_classifier(None,
                                buffer_size,
                                &None,
                                None,
                                &classifier,
                                &PathBuf::from("../tests/small_example_testing.txt"),
                                &output_path,
                                false).unwrap();
    std::fs::remove_file(&output_path).unwrap();
}
fn bench_taxonomy_example(max_depth: usize, buffer_size: usize) {
    // run with depth 13, buffer_size 512
    let train_path = PathBuf::from("../tests/taxonomy_test_training.txt");
    let test_path = PathBuf::from("../tests/taxonomy_test_testing.txt");
    let output_path = PathBuf::from("../tests/taxonomy_example_predication.tsv");

    let classifier = create_lz_classifier(None, max_depth, &train_path, buffer_size, &Some(4));
    predict_using_lz_classifier(None,
                                buffer_size,
                                &Some(4),
                                Some(11.0),
                                &classifier,
                                &test_path,
                                &output_path,
                                false).unwrap();
    std::fs::remove_file(&output_path).unwrap();
}

fn bench_tiny_example(max_depth: usize, buffer_size: usize) {
    // run with depth 12, buffer_size 512
    let train_path = PathBuf::from("../tests/tiny_training.txt");
    let test_path = PathBuf::from("../tests/tiny_testing.txt");
    let output_path = PathBuf::from("../tests/tiny_sample_predication.tsv");

    let classifier = create_lz_classifier(None, max_depth, &train_path, buffer_size, &None);
    predict_using_lz_classifier(None,
                                buffer_size,
                                &None,
                                None,
                                &classifier,
                                &test_path,
                                &output_path,
                                false).unwrap();
    std::fs::remove_file(&output_path).unwrap();
}

fn bench_meta_genes(max_depth: usize, buffer_size: usize) {
    let output_path = PathBuf::from("../tests/meta_small_sample_predication.tsv");
    // max_depth is 12, for consistency with the small sample.
    let classifier = create_lz_classifier(None, max_depth, &PathBuf::from("../tests/small_example_training.txt"), buffer_size, &None);
    meta_predict_using_lz_classifier(None,
                                     buffer_size,
                                     &classifier,
                                     &PathBuf::from("../data/small_sample_as_meta.fna"),
                                     &output_path,
                                     true,
                                     0).unwrap();
    std::fs::remove_file(&output_path).unwrap();
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("bench_presentation_example_two_seqs 13", |b| b.iter(|| bench_presentation_example_two_seqs(black_box(13))));
    c.bench_function("bench_small_example 12 512", |b| b.iter(|| bench_small_example(black_box(12), black_box(512))));
    c.bench_function("bench_tiny_example 12 512", |b| b.iter(|| bench_tiny_example(black_box(12), black_box(512))));
}

fn long_criterion_benchmark(c: &mut Criterion) {
    c.bench_function("bench_taxonomy_example 13 512", |b| b.iter(|| bench_taxonomy_example(black_box(13), black_box(512))));
    c.bench_function("bench_meta_genes 12 512", |b| b.iter(|| bench_meta_genes(black_box(12), black_box(512))));
}

criterion_group!{
  name = benches;
  config = Criterion::default().measurement_time(Duration::from_secs(60));
  targets = criterion_benchmark
}

criterion_group!{
  name = longbenches;
  config = Criterion::default().measurement_time(Duration::from_secs(500));
  targets = long_criterion_benchmark
}

criterion_main!(benches, longbenches);