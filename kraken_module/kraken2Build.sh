#!/bin/bash

#######
# DESCRIPTION:
#     This bash script installs all of the required files needed for
#     running the EnrichSeq pipeline. This includes phage genomes, 
#     actinodatabase genomes, and also build all of the underlying
#     databases in the pipeline.
#######




# change to the correct directory.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$( dirname "${BASH_SOURCE[0]}" )" 

####
# Description:
#      This function uses ncbi-genome-download to download all 
#      phage genomess into an output genome direrctory.
# TODO:
#      1. add the correct error and warnings.
####
function NCBIPhageGenomeDownload() {
  local genomeDir=$1;
  # download NCBI phage genomes
  while [[ ! -d refseq/viral ]]; do #TODO: ensure this is the right way to make while loop.
    if [[ ! -d refseq/viral ]]; then
      echo "Downloading NCBI phage genomes...";
      #TODO: print out the attempt #
      #TODO: If tries above 3, then exiit script and give advise for nternet connection
      ncbi-genome-download  --formats fasta --assembly-level complete \
                            --genera phage --fuzzy-genus viral --parallel 4 \
                            --flat-output -o ${genomeDir}; 
      # unzip the downloaded genomes
      find ref_genomes/ -name '*.gz' -exec gzip -d {} +;
    else 
      # TODO break while loop, print "Success and output message"
      break
    fi
  done
}

####
# Description:
#      This function calls wget_file() to obtain all ftp
#      DB required files.
####
function DownloadFTPfiles() {
  echo "Running DownloadFTPFiles()";
  local dbDir=$1;
  # Download phagesDB all phages database
  actino_ftp='https://phagesdb.org/media/Actinobacteriophages-All.fasta'
  wget_file ${actino_ftp} Actinobacteriophages-All.fasta Actinobacteriophages-All.fasta
  # Download NCBI taxonomy
  taxdump_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
  wget_file ${taxdump_ftp} new_taxdump.tar.gz ${dbDir}/taxonomy/names.dmp
  # download accession to taxid for WGS
  nucl_wgs_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz'
  wget_file ${nucl_wgs_ftp} nucl_wgs.accession2taxid.gz ${dbDir}/taxonomy/nucl_wgs.accession2taxid 
}

####
# Description:
#      This function uses wget to obtain a file,
#      and performs error handeling if the file
#      is not corrrectly downloaded.
# TODO:
#      1. add the correct error and warnings.
####
function wget_file() {
  local ftp_path=$1;
  local downloaded_file=$2;
  local final_file=$3;
  while [[ ! -f ${downloaded_file} && ! -f ${final_file} ]]; do
    if [[ ! -f ${downloaded_file} && ! -f ${final_file} ]]; then
      # print out attempt number
      echo; echo "      Downloading ${ftp_path}"; echo;
      wget ${ftp_path};
    else;
      # TODO break while loop, print "Success and output message"
      break
    fi
  done
}

####
# TODO:
#      1. split function up into sub functions w/ error handeling.
####
function reorganizeFiles() {
  echo "Running reorganizeFiles()";
  
  # arguments
  local genomeDir=$1;
  local dbDir=$2;
  local actinoOutFile=$3;

  # organize taxdump files into the appropriate location
  if [[ -f new_taxdump.tar.gz ]]; then
    tar -xf new_taxdump.tar.gz --directory ${dbDir}/taxonomy;
    rm new_taxdump.tar.gz;
  fi

  # downloading PhageDB genomes
  if [[ -f nucl_gb.accession2taxid.gz ]]; then
    gzip -d nucl_gb.accession2taxid.gz;
    moveFile nucl_gb.accession2taxid ${dbDir}/taxonomy;
  fi

  # Reformat ActinoDB file to be compatible with kraken
  python ${DIR}/kraken_db_format.py ${dbDir}/taxonomy/names.dmp \
         Actinobacteriophages-All.fasta ${actinoOutFile};
  moveFile ${actinoOutFile} ${genomeDir};
  rm Actinobacteriophages-All.fasta;
}

####
# Description:
#      Moves a file into a directory
# Errors:
#      1. If file empty
#      2. If directory does not exist
# Warnings:
#      1. If directory already contains that file
# TODO:
#      1. add the correct error and warnings as above.
####
function moveFile() {
  local file_to_move=$1;
  local new_directory=$2;
  echo "    moving the file ${file_to_move} to ${new_directory}";
  if [[ -f ${file_to_move} ]]; then
    mv ${file_to_move} ${new_directory}
  fi
}

####
# Description:
#      This method adds genomes to the to-be-build krakenDB
# Errors:
#      1. If python kraken2-build missing
#      2. If ${genomeDir}/ empty
# Warnings:
#      1. if file without find extensions are in ${genomeDir}/
# TODO:
#      1. add the correct error and warnings as above.
####
function addGenomesToDb() {
  echo "Running addGenomesToDb()";
  # arguments
  local genomeDir=$1;
  local dbDir=$2;
  # add all fasta files in the genomes directory to the kraken library
  find ${genomeDir}/ \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
       -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db ${dbDir};
}

####
# Description:
#      This method calls on kraken2-build to create a directory.
# Errors:
#      1. If python kraken2-build missing
####
function buildKrakenDb() {
  echo "Running buildKrakenDb()";
  # arguments
  local dbDir=$1;
  kraken2-build --build --db ${dbDir}; # TODO: throw error if no kraken2-build
}

####
# Description:
#      This method makes a directory if it does not exist.
####
function makeDirectory() {
  local directory_to_make=$1;
  echo "making directory for ${directory_to_make}"
  if [[ ! -d ${directory_to_make} ]]; then
    mkdir ${directory_to_make};
  fi
}

####
# Description:
#      This method uses multifasta2single.py to create sub fasta files.
# Errors:
#      1. If python script missing
#      2. If multifasta doesn't exist
####
function multifasta2fasta(){
  local multi_fasta_file=$1;
  local directory_for_fasta_files=$2;
  local multifasta2fasta_script=multifasta2single.py
  echo "Turning ${multi_fasta_file} into multiple fasta files, then deleting";
  if [[ -f ${multifasta2fasta_script} ]]; then
    if [[ -f ${multi_fasta_file} ]]; then
      python ${DIR}/${multifasta2fasta_script} ${multi_fasta_file} ${directory_for_fasta_files};
      rm ${multi_fasta_file}; #TODO: Should the multifasta be deleted?
    else;
      throw_fatal_error "${multifasta2fasta_script} is missing from ${DIR}";
  else; # ensure this is the correct use of else
    throw_fatal_error "${multi_fasta_file} is missing or not found!!!!";
  fi
}

####
# Description:
#      This method throws a warning then continues.
####
function throw_warning() {
  local warning_msg=$1
  echo; echo "**WARNING:";
  echo "${multifasta2fasta_script} is missing from ${DIR}";
  echo;
}

####
# Description:
#      This method throws an error then exits.
####
function throw_fatal_error() {
  local error_msg=$1
  echo; echo "**ERROR:";
  echo ${error_msg};
  echo "Exiting";
  exit 1; 
}


####
####
# Description: The main function controls the flow of the script.
####
####
function main() {                                                               
  # input arguments                                                             
  if [[ -f ${DIR}/kraken_module.config ]]; then
    source ${DIR}/kraken_module.config;
  else
    echo "The config file (kraken_module.config) is missing!";
    exit 1;
  fi

  # update arguments with full path
  local genomeDir=${DIR}/${genomeDir};
  local dbDir=${DIR}/${dbDir};
  local movedActinoOutFile=${genomeDir}/${actinoOutFile}
  local actinoOutFile=${DIR}/${actinoOutFile};
  echo $actinoOutFile;

  # run the script                                           
  if [[ ! -f ${dbDir}/taxo.k2d ]]; 
  then   
    makeDirectory ${genomeDir};
    makeDirectory ${dbDir};
    makeDirectory ${dbDir}/taxonomy;
    NCBIPhageGenomeDownload ${genomeDir};
    DownloadFTPFiles ${dbDir};                                                
    reorganizeFiles ${genomeDir} ${dbDir} ${actinoOutFile};                   
    addGenomesToDb ${genomeDir} ${dbDir};                                       
    buildKrakenDb ${dbDir};
    multifasta2fasta ${movedActinoOutFile} ${genomeDir};
  else
    echo "The database was already built. If you'd like to rebuild, delete '${dbDir}/taxo.k2d'";
  fi                                                                            
}

echo "Running the kraken2 DB build script";
main;
