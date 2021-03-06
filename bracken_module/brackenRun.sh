#!/bin/bash

usage() {
    echo; echo "Usage: bash $0 --krakendb=krakenDB/ --input=krakenout.kraken --out=bracken_run1.bracken --read 100"
    echo "  --krakendb   Path to the kraken2 database created after running the kraken module"
    echo "  --input     Path to the .kraken file generated from the kraken module"
    echo "  --out       Output prefix for the .bracken file name"
    echo "  --read      Read length"
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
    --krakendb=?*)
        dbDir=${1#*=};;
    --krakendb|krakendb=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --input=?*)
        echo ${1#*=};
        inputKraken=${1#*=};;
    --input|input=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --out=?*)
        outPrefix=${1#*=};;
    --out|out=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --read=?*)
        rlength=${1#*=};;
    --read|read=)
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


function runBracken() {
  echo "Running runBracken()";
  # arguments
  local dbDir=$1;
  local inputKraken=$2;
  local outPrefix=$3;
  local rlength=$4;
  
  bracken -d ${dbDir} -i ${inputKraken} -o ${outPrefix}.bracken -r ${rlength} -l 'S' -t 10;
}

function main() {
  echo ${dbDir};
  echo ${inputKraken};
  echo ${outPrefix};
  echo ${rlength}; 

  runBracken ${dbDir} ${inputKraken} ${outPrefix} ${rlength};
}

echo "Running the bracken run script";
main;
