#!/usr/bin/env bash

../../code/target/release/GeneZip --verbose -j 8 train-predict --gc 100 --kmer 0 -i training.txt -t testing.txt -o output/prediction_output --lzvalues output/lz_matrix -d 12 2> output/log
