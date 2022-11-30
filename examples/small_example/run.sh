#!/usr/bin/env bash

../../code/target/release/GeneZip -v -i training.txt -t testing.txt -o output/prediction_output -d 12 -j 8 2> output/log
