#!/bin/bash

usage() {
    echo; echo "Usage: bash $0 --read=single --input=inputfasta/simgenomes.fa --threads=4 --out=megahit_20210307"
    echo "  --read      Read type of fasta files [single, paired, long]"
    echo "  --input     Input fasta file"
    echo "  --threads   Number of threads"
    echo "  --out       Output directory for metaphlan results"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

# check that all the required arguments are used
if [ $# -gt 4 ] || [ $# -lt 4 ]
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
    --input=?*)
        echo ${1#*=};
        inFasta=${1#*=};;
    --input|input=)
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
  local threads=$3;
  local outdir=$4;

  if [[ ${readType} == "single" ]]; then
    megahit -r ${inFasta} -t ${threads} -m 1e9 -o ${outdir} --out-prefix megahit_out
  fi

}

function main() {
  echo ${readType};
  echo ${inFasta};
  echo ${threads};
  echo ${outdir};

  runMegahit ${readType} ${inFasta} ${threads} ${outdir};
}

echo "Running the megahit run script";
main;
