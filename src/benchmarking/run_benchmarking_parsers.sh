#!/bin/bash
METADATA = "benchmarking_metadata"

function run_directory_parser() {
    truth_dir = "./tests/"
    result_dir = "./results/"
    echo "Running directory parser."
    python3 directory_parser.py -t ${truth_dir} -r ${result_dir} -o ${METADATA}
}


function run_result_parsers() {
    prefix = "$(date +'%m-%d-%Y_%H%M')"
    echo $"Running result file parsers."
    python3 benchmark.py -f ${METADATA} -o ${prefix}
}

function main() {
    conda activate enrichseq;
    run_directory_parser;
    run_result_parsers;
}
