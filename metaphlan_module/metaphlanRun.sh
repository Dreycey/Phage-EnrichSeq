#!/bin/bash

usage() {
    echo; echo "Usage: bash $0 --input=assembly/megahitout.fa --type=fasta --out=profiled_metagenome.txt --nproc 2"
    echo "  --input     Path to the kraken report file generated from the kraken module"
    echo "  --type   input file type [fasta, fastq, sam, bowtie2out]"
    echo "  --out       Output file for metaphlan results"
    echo "  --nproc      number of processors to use"
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
    --input=?*)
        inFasta=${1#*=};;
    --input|input=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --type=?*)
        echo ${1#*=};
        inputType=${1#*=};;
    --type|type=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --out=?*)
        outfile=${1#*=};;
    --out|out=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --nproc=?*)
        nproc=${1#*=};;
    --nproc|nproc=)
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


function runMetaphlan() {
  echo "Running runMetaphlan()";
  # arguments
  local inFasta=$1;
  local inputType=$2;
  local outfile=$3;
  local nproc=$4;

  # TODO: change this to incorporate custom DB
  metaphlan ${inFasta} --input_type ${inputType} --bowtie2out metagenome.bowtie2.bz2 --nproc ${nproc} -o ${outfile} \
  --ignore_eukaryotes --ignore_archaea --add_viruses;
}

function main() {
  source metaphlan_module.config; # to obtain metaphlan db file
  echo ${inFasta};
  echo ${inputType};
  echo ${outfile};
  echo ${nproc};

  if [[ -f metagenome.bowtie2.bz2 ]]; then
    rm metagenome.bowtie2.bz2
  fi
  runMetaphlan ${inFasta} ${inputType} ${outfile} ${nproc};
}

echo "Running the metaphlan run script";
main;
