


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
    local blastDB="outputMulti3.fa";
    local inFasta="simulatedgenomes_illumina.fa";
    local blastOut="blastOut.txt";

    # running underlying methods
    runBlastN ${blastDB} ${inFasta} ${blastOut};
}

echo "Running the BLASTDB run script";
main;
