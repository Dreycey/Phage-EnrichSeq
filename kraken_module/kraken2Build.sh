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

####
# Description:
#      This function uses ncbi-genome-download to download all 
#      phage genomess into an output genome direrctory.
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
####
function wget_file() {
  local ftp_path=$1;
  local downloaded_file=$2;
  local final_file=$3;
  while [[ ! -f ${downloaded_file} && ! -f ${final_file} ]]; do
    if [[ ! -f ${downloaded_file} && ! -f ${final_file} ]]; then
      # print out attempt number
      echo "Downloading ${ftp_path}";
      wget ${ftp_path};
    else;
      # TODO break while loop, print "Success and output message"
      break
    fi
  done
}


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
  if [[ -f nucl_wgs.accession2taxid.gz ]]; then
    gzip -d nucl_wgs.accession2taxid.gz;
    moveFile nucl_wgs.accession2taxid ${dbDir}/taxonomy;
  fi

  # Reformat Action file to be compatible with kraken
  python ${DIR}/kraken_db_format.py ${dbDir}/taxonomy/names.dmp \
         Actinobacteriophages-All.fasta ${actinoOutFile};
  moveFile ${actinoOutFile} ${genomeDir};
  rm Actinobacteriophages-All.fasta;
}

function moveFile() {
  local file_to_move=$1;
  local new_directory=$2;
  echo "moving the file ${file_to_move} to ${new_directory}";
  if [[ -f ${file_to_move} ]]; then
    mv ${file_to_move} ${new_directory}
  fi
}

function addGenomesToDb() {
  echo "Running addGenomesToDb()";
  # arguments
  local genomeDir=$1;
  local dbDir=$2;
  # add all fasta files in the genomes directory to the kraken library
  find ${genomeDir}/ \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) \
       -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db ${dbDir};
}

function buildKrakenDb() {
  echo "Running buildKrakenDb()";
  # arguments
  local dbDir=$1;
  kraken2-build --build --db ${dbDir};
}

function makeDirectory() {
  local directory_to_make=$1;
  echo "making directory for ${directory_to_make}"
  if [[ ! -d ${directory_to_make} ]]; then
    mkdir ${directory_to_make};
  fi
}

function multifasta2fasta(){
  local multi_fasta_file=$1;
  local directory_for_fasta_files=$2;
  echo "Turning ${multi_fasta_file} into multiple fasta files";
  python ${DIR}/multifasta2single.py ${multi_fasta_file} ${directory_for_fasta_files};
  rm ${multi_fasta_file};
}

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
