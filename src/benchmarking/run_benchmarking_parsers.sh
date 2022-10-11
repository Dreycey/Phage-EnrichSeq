#!/bin/bash
METADATA="benchmarking_metadata";

function run_directory_parser() {
    truth_dir="./tests/";
    result_dir="./results/";
    echo "Running directory parser.";
    echo ${truth_dir};
    if [[ ! -d ${truth_dir} ]] || [[ ! -d ${result_dir} ]]; then
	echo "${truth_dir} or ${result_dir} does not exist.";
        exit 404;
    else
        python3 tools/Phage-EnrichSeq/src/benchmarking/directory_parser.py -t ${truth_dir} -r ${result_dir} -o ${METADATA};
    fi
}


function run_result_parsers() {
    prefix="$(date +'%m-%d-%Y_%H%M')";
    echo "Running result file parsers.";
    python3 tools/Phage-EnrichSeq/src/benchmarking/benchmark.py -f ${METADATA} -o ${prefix};
}

function main() {  
    #conda activate enrichseq;
    run_directory_parser;
    #run_result_parsers;
}

main
