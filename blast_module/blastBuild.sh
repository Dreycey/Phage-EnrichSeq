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
function fixMultifastaNames() {

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
function addActinoToMultiFasta() {                                                 

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
    local outputMulti="outputMulti.fa";
    local phagedbpath="/users/dreyceyalbin/desktop/phage-enrichseq_old/downloadphagegenomes_module/phagedbgenomes/";
    # running underlying methods
    createMultiFasta ${phagedbpath} ${outputMulti};
    createBlastDB ${outputMulti};
}

echo "Running the BLASTDB build script";
main;
