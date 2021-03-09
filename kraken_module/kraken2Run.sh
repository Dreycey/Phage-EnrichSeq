usage() {
    echo; echo "Usage: bash $0 --krakendb=krakenDB/ --queryfasta=inputfasta/simulatedgenomes_illumina.fa --report=kraken.report.txt --out=krakenout.kraken"
    echo "  --krakendb   Path to the kraken2 database created using the config and kraken2Build.sh"
    echo "  --queryfasta     Path to the in fasta file"
    echo "  --report    Report file for kraken results"
    echo "  --out       File name for .kraken output"
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
    --queryfasta=?*)
        echo ${1#*=};
        inFasta=${1#*=};;
    --queryfasta|queryfasta=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --report=?*)
        echo ${1#*=};
        reportFile=${1#*=};;
    --report|report=)
        echo "$0: missing argument for '$1' option"
        usage
        exit 1;;
    --out=?*)
        krakenOutFile=${1#*=};;
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

  echo $dbDir;
  echo $inFasta;
  echo $reportFile;
  echo $krakenOutFile;

  # input arguments
  # local dbDir="krakenDB";
  # local inFasta="simgenomes1_illumina.fa";
  # local reportFile="krakenOut.txt";
  # local krakenOutFile="krakenOut.kraken";

  runKraken2 ${dbDir} ${inFasta} ${reportFile} ${krakenOutFile};
}


echo "Running the kraken2 run script";
main;
