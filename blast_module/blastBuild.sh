#!/bin/bash
#
# Perform the build for the Blast database


#######################################
# createMultiFasta() 
# creates a multifasta from a directory of fasta files.
# Globals:
#   None
# Arguments:
#   fasta file directory
#   output multifasta file
# Outputs:
#   Creates multifasta file with the second input argument
#######################################
function createMultiFasta() {
  echo "Running createMultiFasta()";
  # arguments
  local genomeDir=$1;
  local outputMultiFasta=$2;
  # iteratively add genomes to the multifasta
  if [[ ! -f ${outputMultiFasta} ]]; then
    for file in ${genomeDir}/*; do 
      #echo ${file}; 
      cat ${file} >> ${outputMultiFasta}; 
    done; 
  fi
}


#######################################                                         
# fixMultifastaNames()                                                     
# The fixMultifastaNames() function makes sure no spaces are present in names                      
# Globals:                                                                      
#   None                                                                        
# Arguments:                                                                    
#   multifasta file path                                                   
# Outputs:                                                                      
#   Creates a new multifasta with corrected names
#######################################
function fixMultifastaNames() {
  echo "Running fixMultifastaNames()";
  # arguments
  local inputMultiFasta=$1;
  local outputMultiFasta=$2;
  # run the script
  if [[ ! -f ${outputMultiFasta} ]]; then 
    python fixMultiFastaNames.py ${inputMultiFasta} ${outputMultiFasta}; 
  fi
}

#######################################                                         
# addActinoToMultiFasta()                                                         
# this function adds the actino database to the multiFasta
# Globals:                                                                      
#   None                                                                        
# Arguments:                                                                    
#   multifasta file path
#   actinoDBfile 
# Outputs:                                                                      
#   Creates athe finalized multi fasta file
#######################################
function addActinoToMultiFasta() {                                                 
  echo "Running addActinoToMultiFasta()";
  # arguments
  local inputMultiFasta=$1;
  local inputActinoDB=$2;
  local outputMultiFasta=$3;
  # run the script
  if [[ ! -f ${outputMultiFasta} ]]; then
    cat ${inputMultiFasta} >> ${outputMultiFasta};
    cat ${inputActinoDB} >> ${outputMultiFasta}; 
  fi
} 


#######################################                                         
# createBlastDB()                                                      
# creates a blast database from the multifasta file.                         
# Globals:                                                                      
#   None                                                                        
# Arguments:                                                                    
#   multifasta file path                            
#   output blastDB name                                                    
# Outputs:                                                                      
#   Creates a balstDB using the input multifasta.                  
#######################################
function createBlastDB() { 
  echo "Running createBlastDB()";
  # iniput arguments
  local inputMultiFasta=$1;
  # make the blast database
  if [[ ! -f ${inputMultiFasta}.nhr || \
        ! -f ${inputMultiFasta}.nin || \
        ! -f ${inputMultiFasta}.nsq ]]; then
    makeblastdb -in ${inputMultiFasta} -dbtype nucl 
  fi
}

#######################################
# Main function runs all of the other methods for the build
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Calls on functions for creating the bash DB
#######################################
function main(){
    # input arguments
    local outputMulti_1="outputMulti1.fa";
    local outputMulti_2="outputMulti2.fa";
    local outputMulti_3="outputMulti3.fa";
    local actinoFilePath="/Users/dreyceyalbin/Desktop/Phage-EnrichSeq_old/downloadPhageGenomes_module/Actinobacteriophages-All.fasta";
    local phagedbpath="/users/dreyceyalbin/desktop/phage-enrichseq_old/downloadphagegenomes_module/phage_genomes/";
    
    # running underlying methods
    createMultiFasta ${phagedbpath} ${outputMulti_1};
    addActinoToMultiFasta ${outputMulti_1} ${actinoFilePath} ${outputMulti_2};
    fixMultifastaNames ${outputMulti_2} ${outputMulti_3};
    createBlastDB ${outputMulti_3};
}

echo "Running the BLASTDB build script";
main;
