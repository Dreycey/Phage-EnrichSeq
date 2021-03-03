

function runKraken2() {
  echo "Running runKraken2()";
  # arguments
  local dbDir=$1;
  local inFasta=$2;
  local reportFile=$3;
  local krakenOutFile=$4;

  kraken2 --use-names --threads 4 --db ${dbDir} --report ${reportFile} ${inFasta} > ${krakenOutFile}
}


function main() {
  # input arguments
  local dbDir="krakenDB";
  local inFasta="simgenomes1_illumina.fa";
  local reportFile="krakenOut.txt";
  local krakenOutFile="krakenOut.kraken";

  runKraken2 ${dbDir} ${inFasta} ${reportFile} ${krakenOutFile};
}


echo "Running the kraken2 run script";
main;
