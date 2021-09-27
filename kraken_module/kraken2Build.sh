#!/bin/bash
# Build a kraken2 database

#######
# DESCRIPTION:
#     This bash script installs all of the required files needed for
#     running the EnrichSeq pipeline. This includes phage genomes, 
#     actinodatabase genomes, and also build all of the underlying
#     databases in the pipeline.
#######


# change to the correct directory.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function downloadRequiredFiles() {
  echo "Running downloadRequiredFiles()";
 
  #arguments
  local dbDir=$1;

  # download NCBI phage genomes
  if [[ ! -d refseq/viral ]]; then
    echo "Downloading NCBI phage genomes...";
    ncbi-genome-download  --formats fasta --assembly-level complete \
                          --genera phage --fuzzy-genus viral --parallel 4 \
                          --flat-output -o ${genomeDir}; 
    # unzip the downloaded genomes
    find ref_genomes/ -name '*.gz' -exec gzip -d {} +;
  fi

  # Download phagesDB all phages database
  if [[ ! -f Actinobacteriophages-All.fasta ]]; then
    echo "Downloading Actinobacteriophage genomes...";
    wget https://phagesdb.org/media/Actinobacteriophages-All.fasta;
  fi
  
  # Download NCBI taxonomy
  if [[ ! -f ${dbDir}/taxonomy/names.dmp && ! -f new_taxdump.tar.gz ]]; then 
    echo "Downloading NCBI taxonomy...";
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz;
  fi

  if [[ ! -f ${dbDir}/taxonomy/nucl_wgs.accession2taxid && ! -f nucl_wgs.accession2taxid.gz ]]; then
    echo "Downloading accession to taxon ID conversion...";
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz;
  fi
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
  fi
  moveFile nucl_wgs.accession2taxid ${dbDir}/taxonomy;
  moveFile nucl_wgs.accession2taxid ${dbDir}/taxonomy;

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
  if [[ -f ${DIR}/kraken_module.config ]];
  then
    source ${DIR}/kraken_module.config;
  else
    echo "The config file (kraken_module.config) is missing!";
    exit 1;
  fi
  # update arguments with full path
  genomeDir=${DIR}/${genomeDir};
  dbDir=${DIR}/${dbDir};
  movedActinoOutFile=${genomeDir}/${actinoOutFile}
  actinoOutFile=${DIR}/${actinoOutFile};
  echo $actinoOutFile;
  # run the script                                           
  if [[ ! -f ${DIR}/${dbDir}/taxo.k2d ]]; 
  then   
    makeDirectory ${genomeDir};
    makeDirectory ${dbDir};
    makeDirectory ${dbDir}/taxonomy;                                      
   downloadRequiredFiles;                                                      
   reorganizeFiles ${genomeDir} ${dbDir} ${actinoOutFile};
    multifasta2fasta ${movedActinoOutFile} ${genomeDir};                   
    addGenomesToDb ${genomeDir} ${dbDir};                                       
    buildKrakenDb ${dbDir};
  else
    echo "The database was already built. If you'd like to rebuild, delete '${dbDir}/taxo.k2d'";
  fi                                                                            
}

echo "Running the kraken2 DB build script";
main;
