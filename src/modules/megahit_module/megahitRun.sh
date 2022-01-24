#!/bin/bash

usage() {
    echo; echo "Usage: bash $0 --read=single --input=inputfasta/simgenomes.fa --threads=4 --out=megahit_20210307"
    echo "  --read      Read type of fasta files [single, paired, long]"
    echo "  --input1     Input fasta file 1 (just this is single end)"
    echo "  --input2     Input fasta file 2"
    echo "  --threads   Number of threads"
    echo "  --out       Output directory for metaphlan results"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

# check that all the required arguments are used
if [ $# -gt 5 ] || [ $# -lt 4 ]
then
    usage
fi

# parse the commands
while true
do
    case $1 in
    --help|-h)
        usage
        exit;;
    --read=?*)
        readType=${1#*=};;
    --read|read=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --input1=?*)
        echo ${1#*=};
        inFasta=${1#*=};;
    --input2=?*)
        echo ${1#*=};
        inFasta2=${1#*=};;
    --input1|input1=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --input2|input2=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --threads=?*)
        threads=${1#*=};;
    --threads|threads=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --out=?*)
        outdir=${1#*=};;
    --out|out=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --)
        shift
        break;;
    -?*)
        echo "$0: invalid option: $1"
        usage
        exit 1;;
    *)
        break
    esac
    shift
done


function runMegahit() {
  echo "Running runMegahit()";
  # arguments
  local readType=$1;
  local inFasta=$2;
  local inFasta2=$3;
  local threads=$4;
  local outdir=$5;

  if [[ ${readType} == "single" ]]; then
    megahit -r ${inFasta} -t ${threads} -m 1e9 -o ${outdir} --out-prefix megahit_out
  elif [[ ${readType} == "paired" ]]; then
    megahit -1 ${inFasta} -2 ${inFasta2} -t ${threads} -m 1e9 -o ${outdir} --out-prefix megahit_out
  fi

}

function main() {
  echo ${readType};
  echo ${inFasta};
  echo ${inFasta2};
  echo ${threads};
  echo ${outdir};

  runMegahit ${readType} ${inFasta} ${inFasta2} ${threads} ${outdir};
}

echo "Running the megahit run script";
main;
