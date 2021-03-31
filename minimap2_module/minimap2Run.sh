#!/bin/bash
#
# Run blastN

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                         

usage() {
    echo; echo "Usage: bash $0 --blastdb=blastdb/outputMulti3.fa --queryfasta=blast_testfiles/simulatedgenomes_illumina.fa --out=blastout.txt"
    echo "  --blastdb   Path to the blastn databased created using the config and blastBuild.sh"
    echo "  --queryfasta     Path to the in fasta file"
    echo "  --out       Path to blastN output"
    echo "  -h, --help  Print this help message out"; echo;
    exit 1;
}

# check that all the required arguments are used
if [ $# -gt 3 ] || [ $# -lt 3 ]
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
    --blastdb=?*)
        blastDB=${1#*=};;
    --blastdb|blastdb=)
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
    --out=?*)
        blastOut=${1#*=};;
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


#######################################  
# runBlastN()                                     
# THis function runs the blastN script, thus using blastN to map sequences
# Globals:                                                                      
#   None                                                                        
# Arguments:   
#   blastdb name                                                                 
#   input fasta (query)
#   outputfile      
# Outputs:                                                                      
#   This function returns a blastn output file    
#######################################  
function  runBlastN(){
  echo "running runBlastN()";
  # arguments
  local blastDB=$1;
  local inputFasta=$2;
  local blastOutput=$3;
  # running blastN
  if [[ ! -f ${blastOutput} ]]; then
    blastn -db ${blastDB} -query ${inputFasta} -out ${blastOutput} \
    -outfmt "6 qseqid sseqid sblastnames score evalue pident length mismatch gap open gaps sseq";
  fi
}

#######################################
# Main function runs all of the other methods for the build
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Calls on functions for creating running blastN
#######################################
function main(){
    # input arguments
   # local blastDB="blastdb/outputMulti3.fa";
   # local inFasta="blast_testfiles/simulatedgenomes_illumina.fa";
   # local blastOut="blastOut.txt";

    echo $blastDB;
    echo $inFasta;
    echo $blastOut;

    # set blastdb path if not already set    
    if [[ -z ${BLASTDB} ]]; then
      export BLASTDB=$BLASTDB:"blastdb/"
    fi

    # running underlying methods
    runBlastN ${blastDB} ${inFasta} ${blastOut};
}

echo "Running the BLASTDB run script";
main;
