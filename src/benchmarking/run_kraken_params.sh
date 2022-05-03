usage() {
    echo; echo "Usage: bash $0 --input=simulatedgenomes_illumina.fa --refdir=database/ref_genomes/ --kmer_length=35 --min_length=31 --seed=7 --assembly=y"
    echo "  --input           Path to the in fasta file"
    echo "  --refdir          Path to the kraken2 database created using the config and kraken2Build.sh"
    echo "  --kmer_length     Size of kmers to create (for database build)"
    echo "  --min_length      Length of minimizer (for database build)"
    echo "  --seed            Spaced seeds size (for database build)"
    echo "  --assembly        y or n (with assembly or without)"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}


function buildDatabase() {
     local dbDir=$1;
     # Flag options are the parameters
     echo "Building the database using kraken2-build"
     # Add-to-library must already be complete. How to check this?

     if ! [[ -x "$(command -v kraken2-build)" ]]; then
       throw_fatal_error "kraken2-build is missing, make sure to activate conda environment for enrichseq. Read the README.md";
     else
       kraken2-build --build --db ${dbDir};
     fi
}

function runMegahit() {
     #arguments

}

function runKraken2() {
    # arguments
   local dbDir=$1;
   local inFasta=$2;
   local kmerLength=$3;
   local minimizerLength=$4;
   local seed=$5;

  if [[ -d ./kraken_output ]]; then
       mkdir "./kraken_output"
  fi
  local krakenOutFile='.kraken_out/K${kmerLength}L${minimizerLength}S${seed}.kraken'

  if [[ -f ${dbDir}/taxo.k2d ]]; then
    kraken2 --use-names --threads 4 --db ${dbDir} --report kraken.report ${inFasta} > ${krakenOutFile}
  else
    echo "Kraken database must be built first"
  fi
}

function deleteHash() {
     # TODO: Delete hash.k2d and seqid2taxid.map files
}

function main() {
     echo $inFasta;
     echo $refDir;
     echo $kmerLength;
     echo $minimizerLength;
     echo $seed;
     echo $assembly;
     local dbDir='krakenParamDB';

     buildDatabase; ## Which params to pass?
     if [[ ${assembly} == 'y' ]]; then
          runMegahit ${};
     runKraken2 ${dbDir} ${inFasta} ${kmerLength} ${minimizerLength} ${seed};
}

echo "Running kraken2 parameter sweep";
main;
