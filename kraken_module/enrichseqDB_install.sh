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
#      This function calls wget_file() to obtain all ftp DB required files.
####
function downloadFTPFiles() {
  add_to_log "";
  add_to_log "RUNNING: downloadFTPFiles()";
  local dbDir=$1;
  # Download phagesDB all phages database
  actino_ftp='https://phagesdb.org/media/Actinobacteriophages-All.fasta';
  wget_file ${actino_ftp} Actinobacteriophages-All.fasta Actinobacteriophages-All.fasta;
  # Download NCBI taxonomy
  taxdump_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz';
  wget_file ${taxdump_ftp} new_taxdump.tar.gz ${dbDir}/taxonomy/names.dmp;
  # download accession to taxid for WGS
  nucl_wgs_ftp='https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz';
  wget_file ${nucl_wgs_ftp} nucl_gb.accession2taxid.gz ${dbDir}/taxonomy/nucl_gb.accession2taxid;
}

####
# Description:
#      This function uses ncbi-genome-download to download all phage genomes.
# Notes:
#      This function must be the first to add to ${genomeDir}
# Errors:
#      1. If ${genomeDir} not empty
####
function NCBIPhageGenomeDownload() {
  add_to_log "";
  add_to_log "RUNNING: NCBIPhageGenomeDownload()";
  add_to_log "NOTE: This step takes a while ... (typically about 30 minutes)";
  local genomeDir=$1;
  COUNTER=1
  MAX_ATTEMPTS=3
  # ensure directory is empty
  if [[ "$(ls -A ${genomeDir})" ]]; then
    throw_warning "NCBIPhageGenomeDownload() ${genomeDir} SHOULD BE deleted before rerunning the script! Delete this then try again.";
  fi
  # download NCBI phage genomes
  while [[ ! "$(ls -A ${genomeDir})" ]]; do
    if [[ $COUNTER -gt $MAX_ATTEMPTS ]]; then break; fi;
    if [[ ! "$(ls -A ${genomeDir})" ]]; then
      # DOWNLOAD GENOMES
      add_to_log "Downloading NCBI phage genomes...";
      add_to_log "Attempt # ${COUNTER}";
      ncbi-genome-download  --formats fasta --assembly-level complete \
                            --genera phage --fuzzy-genus viral --parallel 4 \
                            --flat-output -o ${genomeDir}; 
      add_to_log "SUCCESS: genomes downloaded";
      # UNZIP GENOMES
      add_to_log "Unzipping the downloaded NCBI phage genomes...";
      find ref_genomes/ -name '*.gz' -exec gzip -d {} +;
      add_to_log "SUCCESS: unzipping genomes worked";
    else 
      break
    fi
    let COUNTER++
  done;
}

####
# Description:
#      This function uses wget to obtain a file.
# Errors:
#      1. If file not downloaded after preset 25 times.
####
function wget_file() {
  add_to_log "";
  add_to_log "RUNNING: wget_file()";
  local ftp_path=$1;
  local downloaded_file=$2;
  local final_file=$3;
  if [[ ! -f ${downloaded_file} && ! -f ${final_file} ]]; then
    add_to_log "      Downloading ${ftp_path}";
    wget --tries=25 ${ftp_path};
    if [[ -f ${downloaded_file} ]]; then
      add_to_log "SUCCESS: wget_file() downloaded ${downloaded_file}";
    else
      throw_fatal_error "wget_file() could not downlod ${downloaded_file}, try again (with internet!)";
    fi
  else
    add_to_log "wget_file() already downloaded ${downloaded_file}";
  fi
}

####
# TODO:
#      1. split function up into sub functions w/ error handling.
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

  # Reformat the NCBI files to compatible with kraken (and PathOrganizer)
  python ${DIR}/ncbi2krakenHeader.py ${dbDir}/taxonomy/nucl_gb.accession2taxid ${genomeDir}
  
  # Reformat ActinoDB file to be compatible with kraken
  python ${DIR}/kraken_db_format.py ${dbDir}/taxonomy/names.dmp \
         Actinobacteriophages-All.fasta ${actinoOutFile};
  moveFile ${actinoOutFile} ${genomeDir};
  rm Actinobacteriophages-All.fasta;
}

####
# Description:
#      This method adds genomes to the to-be-build krakenDB
# Errors:
#      1. If python kraken2-build missing
#      2. If ${genomeDir}/ empty
# Warnings:
#      1. if file without find extensions are in ${genomeDir}/
####
function addGenomesToDb() {
  add_to_log "";
  add_to_log "RUNNING: addGenomesToDb()";
  # arguments
  local genomeDir=$1;
  local dbDir=$2;
  if ! [[ -x "$(command -v kraken2-build)" ]]; then 
    throw_fatal_error "kraken2-build is missing, make sure environment for enrichseq is used. Read the README.md";
  elif [[ ! "$(ls -A ${genomeDir})" ]]; then
    throw_fatal_error "The genome directory is empty: ${genomeDir}";
  else
    # add all fasta files in the genomes directory to the kraken library
    find ${genomeDir}/ \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
        -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db ${dbDir};
    add_to_log "SUCCESS: addGenomesToDb() added DB files from ${genomeDir}";
  fi
}

####
# Description:
#      This method calls on kraken2-build to create a directory.
# Errors:
#      1. If python kraken2-build missing
####
function buildKrakenDb() {
  add_to_log "";
  add_to_log "RUNNING: buildKrakenDb()";
  # arguments
  local dbDir=$1;
  if ! [[ -x "$(command -v kraken2-build)" ]]; then 
    throw_fatal_error "kraken2-build is missing, make sure environment for enrichseq is used. Read the README.md";
  else 
    kraken2-build --build --db ${dbDir};
    add_to_log "SUCCESS: buildKrakenDb() added files from ${dbDir}";
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
  add_to_log "";
  add_to_log "RUNNING: multifasta2fasta()";
  add_to_log "Turning ${multi_fasta_file} into multiple fasta files, then deleting";
  if [[ -f ${multifasta2fasta_script} ]]; then
    if [[ -f ${multi_fasta_file} && -s ${multi_fasta_file} ]]; then
      python ${DIR}/${multifasta2fasta_script} ${multi_fasta_file} ${directory_for_fasta_files};
      add_to_log "SUCCESS: multifasta2fasta() finished creating fasta files from ${multi_fasta_file}";
      rm ${multi_fasta_file}; #TODO: Should the multifasta be deleted?
      add_to_log "SUCCESS: multifasta2fasta() deleted ${multi_fasta_file}";
    else
      throw_fatal_error "${multifasta2fasta_script} is missing (or empty) from ${DIR}";
    fi
  else
    throw_fatal_error "${multi_fasta_file} is missing or not found in ${DIR}";
  fi
}

####
# Description: 
#      Deletes duplicate genomes in the genome directory
# Errors:
#      1. If db_cleaner.py not found
####
function delete_duplicates() {
  add_to_log "";
  add_to_log "RUNNING: delete_duplicates()";
  local db_cleaner_py='db_cleaner.py'
  local genome_dir=$1
  if [[ ! -f ${DIR}/${db_cleaner_py} ]]; then
    throw_fatal_error "${db_cleaner_py} is not found in ${DIR}! This must be in ${DIR} for the script to work!";
  else
    python ${db_cleaner_py} ${genome_dir}
    add_to_log "SUCCESS: delete_duplicates() deleted duplicate genomes in ${genome_dir}";
  fi
}

####
# Description: (Helper Function)
#      Moves a file into a directory
# Errors:
#      1. If file empty
#      2. If directory does not exist
# Warnings:
#      1. If directory already contains that file
####
function moveFile() {
  add_to_log "";
  add_to_log "RUNNING: moveFile()";
  add_to_log "    moving the file ${file_to_move} to ${new_directory}";
  local file_to_move=$1;
  local new_directory=$2;
  if [[ ! -f ${file_to_move} ]]; then
    throw_fatal_error "moveFile() the following FILE does not exist: ${file_to_move}";
  elif [[ ! -d ${new_directory} ]]; then
    throw_fatal_error "moveFile() the following DIR does not exist: ${new_directory}";
  elif [[ -f ${new_directory}/${file_to_move} ]]; then
    throw_warning "NOTE: moveFile() the file is already in the directory, adding new one...";
    mv ${file_to_move} ${new_directory};
    add_to_log "SUCCESS: moveFile() ${file_to_move} moved to ${new_directory}";
  else
    mv ${file_to_move} ${new_directory};
    add_to_log "SUCCESS: moveFile() ${file_to_move} moved to ${new_directory}";
  fi
}

####
# Description: (Helper Function)
#      This method makes a directory if it does not exist.
####
function makeDirectory() {
  local directory_to_make=$1;
  add_to_log "";
  add_to_log "RUNNING: makeDirectory()";
  add_to_log "making directory for ${directory_to_make}"
  if [[ ! -d ${directory_to_make} ]]; then
    mkdir ${directory_to_make};
    add_to_log "SUCCESS: makeDirectory() made ${directory_to_make}";
  else
    throw_warning "NOTE: makeDirectory() already made ${directory_to_make}";
  fi
}

####
# Description: (Helper Function)
#      This method throws a warning then continues.
####
function throw_warning() {
  local warning_msg=$1
  echo; echo "**WARNING**:";
  echo ${warning_msg};
  add_to_log "**WARNING**: ${warning_msg}";
  echo;
}

####
# Description: (Helper Function)
#      This method throws an error then exits.
####
function throw_fatal_error() {
  local error_msg=$1
  echo; echo "**ERROR**:";
  echo ${error_msg};
  add_to_log "**ERROR**: ${error_msg}";
  echo "Exiting";
  exit 1; 
}

####
# Description: (Helper Function)
#      This function adds a message and time stamp to a log file.
####
function add_to_log() {
  local msg=$1
  DATE=`date "+%Y%m%d"`
  DATE_WITH_TIME=`date "+%Y%m%d-%H%M%S"`;
  echo ${DATE_WITH_TIME} ${msg}
  if [ ! -z ${LOG_FILE} ]; then
    echo ${DATE_WITH_TIME} ${msg} >> ${LOG_FILE};
  else
    echo ${DATE_WITH_TIME} ${msg} >> DB_BUILD_${DATE}.log;
  fi
}

####
# run_PIPELINE()
# Description: The pipeline for building the database.
####
function run_PIPELINE() {
    makeDirectory ${genomeDir};
    makeDirectory ${dbDir};
    makeDirectory ${dbDir}/taxonomy;
    downloadFTPFiles ${dbDir};
    NCBIPhageGenomeDownload ${genomeDir};                                                
    reorganizeFiles ${genomeDir} ${dbDir} ${actinoOutFile};                   
    multifasta2fasta ${movedActinoOutFile} ${genomeDir};
    delete_duplicates ${genomeDir};
    addGenomesToDb ${genomeDir} ${dbDir};                                       
    buildKrakenDb ${dbDir};
}

####
# MAIN()
# Description: The main function controls the flow of the script.
####
function main() {                                                               
  # input arguments                                                             
  if [[ -f ${DIR}/kraken_module.config ]]; then
    source ${DIR}/kraken_module.config;
  else
    throw_fatal_error "The config file (kraken_module.config) is missing!";
  fi

  # update arguments with full path
  local genomeDir=${DIR}/${genomeDir};
  local dbDir=${DIR}/${dbDir};
  local movedActinoOutFile=${genomeDir}/${actinoOutFile}
  local actinoOutFile=${DIR}/${actinoOutFile};

  # start log file
  if [[ -f  ${LOG_FILE} ]]; then
    add_to_log "The LOG file (${LOG_FILE}) already exists. If you'd like to rebuild, enter 'yes'";
    read varname
    if [[ $varname == 'yes' || $varname == 'y' ]]; then
      rm ${LOG_FILE}
      add_to_log "SUCCESS: deleted ${LOG_FILE}";
    fi
  fi

  # Build genome directory     
  if [[ "$(ls -A ${genomeDir})" ]]; then
    add_to_log "The genome directory (${genomeDir}) already exists. If you'd like to rebuild, enter 'yes'";
    read varname
    if [[ $varname == 'yes' || $varname == 'y' ]]; then
      rm -rf ${genomeDir}/
      add_to_log "SUCCESS: deleted ${genomeDir}/";
    fi
  fi

  # Build DB     
  if [[ "$(ls -A ${dbDir})" ]]; then
    add_to_log "The DB (${dbDir}) already exists. If you'd like to rebuild, enter 'yes'";
    read varname
    if [[ $varname == 'yes' || $varname == 'y' ]]; then
      rm -rf ${dbDir}/
      add_to_log "SUCCESS: deleted ${dbDir}/";
    fi
  fi

  # run the DB building pipeline                                     
  if [[ ! -f ${dbDir}/taxo.k2d ]]; then
    run_PIPELINE;
  else
    add_to_log "The database was already built. If you'd like to rebuild, enter 'yes'";
    read varname
    if [[ $varname == 'yes' || $varname == 'y' ]]; then
      add_to_log "you said ${varname} - Rebuilding DB and rerunning the script!";
      rm -rf ${dbDir}/
      add_to_log "SUCCESS: deleted ${dbDir}/";
      run_PIPELINE;
    else 
      add_to_log "you said ${varname} - Exiting since DB rebuild not desired";
    fi
  fi                                                                            
}

echo; echo "Running the EnrichSeq database build script"; echo;
main;
