#!/bin/bash
# Build a kraken2 database

function downloadRequiredFiles() {
  echo "Running downloadRequiredFiles()";
 
  #arguments
  local dbDir=$1;

  # download NCBI phage genomes
  if [[ ! -d refseq/viral ]]; then
    echo "Downloading NCBI phage genomes..."
    ncbi-genome-download  --formats fasta --assembly-level complete --genera phage --fuzzy-genus viral --parallel 4
  fi

  # Download phagesDB all phages database
  if [[ ! -f Actinobacteriophages-All.fasta ]]; then
    echo "Downloading Actinobacteriophage genomes..."
    wget https://phagesdb.org/media/Actinobacteriophages-All.fasta
  fi
  
  # Download NCBI taxonomy
  if [[ ! -f ${dbDir}/taxonomy/names.dmp && ! -f new_taxdump.tar.gz ]]; then 
    echo "Downloading NCBI taxonomy..."
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
  fi

  if [[ ! -f ${dbDir}/taxonomy/nucl_wgs.accession2taxid && ! -f nucl_wgs.accession2taxid.gz ]]; then
    echo "Downloading accession to taxon ID conversion..."
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
  fi
}


function reorganizeFiles() {
  echo "Running reorganizeFiles()"
  
  # arguments
  local genomeDir=$1;
  local dbDir=$2;
  local actinoOutFile=$3;
  
  mkdir ${genomeDir};
  mkdir ${dbDir};
  mkdir ${dbDir}/taxonomy;

  # Move all .fna.gz files into the container genomes folder then unzip all
  #for file in refseq/viral/*/*fna*; do
   # mv ${file} ${genomeDir};
  #done;
  mv refseq/viral/*/*fna* ${genomeDir};

  gzip -d ${genomeDir}/*;
  rm -rf refseq

  # organize taxdump files into the appropriate location
  if [[ -f new_taxdump.tar.gz ]]; then
    tar -xf new_taxdump.tar.gz --directory ${dbDir}/taxonomy;
    rm new_taxdump.tar.gz
  fi

  if [[ -f nucl_wgs.accession2taxid.gz ]]; then
    gzip -d nucl_wgs.accession2taxid.gz;
    mv nucl_wgs.accession2taxid ${dbDir}/taxonomy;
  fi

  if [[ -f nucl_wgs.accession2taxid ]]; then 
    mv nucl_wgs.accession2taxid ${dbDir}/taxonomy
  fi

  # Reformat Action file to be compatible with kraken
  python kraken_db_format.py ${dbDir}/taxonomy/names.dmp Actinobacteriophages-All.fasta ${actinoOutFile};
  mv ${actinoOutFile} ${genomeDir};
  rm Actinobacteriophages-All.fasta;
}

function addGenomesToDb() {
  echo "Running addGenomesToDb()"
  
  # arguments
  local genomeDir=$1;
  local dbDir=$2;

  # add all fasta files in the genomes directory to the kraken library
  find ${genomeDir}/ \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db ${dbDir};
}

function buildKrakenDb() {
  echo "Running buildKrakenDb()"

  # arguments
  local dbDir=$1;
  
  kraken2-build --build --db ${dbDir};
}

function main() {
  # input arguments
  source kraken_module.config;

  downloadRequiredFiles ${dbDir};
  reorganizeFiles ${genomeDir} ${dbDir} ${actinoOutFile};
  addGenomesToDb ${genomeDir} ${dbDir};
  buildKrakenDb ${dbDir};
}


echo "Running the kraken2 DB build script";
main;

