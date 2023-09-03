#!/usr/bin/env bash

../../code/target/release/GeneZip --verbose -j 3 train-predict --gc 11 -i training.txt -t testing.txt -o output/prediction_output -d 13 2> output/log
