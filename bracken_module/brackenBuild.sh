#!/bin/bash
#
# Build bracken

function buildBracken() {
  echo "Running buildBracken()";
  
  # arguments
  local dbDir=$1;
  local threads=$2;
  local kmer=$3;
  local rlength=$4;
  local krakenPath=$5;


  bracken-build -d ${dbDir} -t ${threads} -k ${kmer} -l ${rlength} -x ${krakenPath}
}

function main() {
  source bracken_module.config;

  #local krakenPath=$(which kraken);

  #if [[ -z "${krakenPath}" ]]; then 
    buildBracken ${dbDir} ${threads} ${kmer} ${rlength} ${krakenPath}
  #else echo "kraken tool not found"
  #fi
}

echo "Running the bracken build script";
main;
