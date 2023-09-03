#!/usr/bin/env bash

../../code/target/release/GeneZip --verbose train-predict --gc 100 --kmer 0 -i training.txt -t testing.txt -o output/prediction_output -d 12 2> output/log
