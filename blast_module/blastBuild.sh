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
  # arguments
  local genomeDir=$1;
  local outputMultiFasta=$2;
  # iteratively add genomes to the multifasta
  if [[ ! -f blastdb/${outputMultiFasta} ]]; then
    echo "Running createMultiFasta()";
    for file in ${genomeDir}/*; do 
      #echo ${file}; 
      cat ${file} >> ${outputMultiFasta}; 
    done; 
  else
    echo "NOT running createMultiFasta()- already done!";
    echo "    Delete the blastdb/ folder if you'd like to rerun";
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
  # arguments
  local inputMultiFasta=$1;
  local outputMultiFasta=$2;
  # run the script
  if [[ ! -f blastdb/${outputMultiFasta} ]]; then 
    echo "Running fixMultifastaNames()";
    python fixMultiFastaNames.py ${inputMultiFasta} ${outputMultiFasta}; 
  else
    echo "NOT running fixMultifastaNames()- already done!";
    echo "    Delete the blastdb/ folder if you'd like to rerun";
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
  # arguments
  local inputMultiFasta=$1;
  local inputActinoDB=$2;
  local outputMultiFasta=$3;
  # run the script
  if [[ ! -f blastdb/${outputMultiFasta} ]]; then
    echo "Running addActinoToMultiFasta()";
    cat ${inputMultiFasta} >> ${outputMultiFasta};
    cat ${inputActinoDB} >> ${outputMultiFasta}; 
  else
    echo "NOT running addActinoToMultiFasta()- already done!";
    echo "    Delete the blastdb/ folder if you'd like to rerun";
  fi
} 


#######################################                                         
# addBlastTaxdb()                                           
# downloads and adds the taxd to the directory so BLASTN can run properly
# Globals:                                                                      
#   None                                                                        
# Arguments:                                                                           
#   URL for the NCBI taxdb
# Outputs:                                                                      
#   Retrieves the blast taxdb, untars, and adds to blastdb/
#######################################   
function addBlastTaxdb() {
  #  arguments
  local blastTaxURL=$1;
  # get the database, add to directory.
  if [[ ! -f "" || ! -f "" ]]; then
    echo "Running addBlastTaxdb()";
    # wget, untar, and rm tar.gz
    wget ${blastTaxURL};
    tar -zxf taxdb.tar.gz;
    rm taxdb.tar.gz;
    # add to blastdb directory, if exists
    if [[ ! -d "blastdb/" ]]; then                                              
      mkdir blastdb/;                                                           
      mv taxdb.bti blastdb/;                                          
    else                                                                        
      mv taxdb.btd blastdb/;
    fi  
  else
    echo "NOT running addBlastTaxdb()- already done!";                          
    echo "    Delete the blastdb/ folder if you'd like to rerun";  
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
  # input arguments
  local inputMultiFasta=$1;
  # make the blast database
  if [[ ! -f blastdb/${inputMultiFasta}.nhr || \
        ! -f blastdb/${inputMultiFasta}.nin || \
        ! -f blastdb/${inputMultiFasta}.nsq ]]; then
    echo "Running createBlastDB()";
    # make the blast database.
    makeblastdb -in ${inputMultiFasta} -dbtype nucl;
    # move the files to the correct location.
    if [[ ! -d "blastdb/" ]]; then
      mkdir blastdb/;
      mv ${inputMultiFasta}* blastdb/;
    else
      mv ${inputMultiFasta}* blastdb/;
    fi 
  else
    echo "NOT running createBlastDB()- already done!";
    echo "    Delete the blastdb/ folder if you'd like to rerun";
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
    source blast_module.config;
    
    # running underlying methods
    createMultiFasta ${phagedbpath} ${outputMulti_1};
    addActinoToMultiFasta ${outputMulti_1} ${actinoFilePath} ${outputMulti_2};
    fixMultifastaNames ${outputMulti_2} ${outputMulti_3};
    addBlastTaxdb ${blasttaxdbURL};
    createBlastDB ${outputMulti_3};
}

echo "Running the BLASTDB build script";
main;
