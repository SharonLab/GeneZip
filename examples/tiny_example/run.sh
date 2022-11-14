#!/usr/bin/env bash

../../code/target/release/GeneZip -i training.txt -t testing.txt -o output/prediction_output -d 12 2> output/log
