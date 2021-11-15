![EnrichSeq](figures/EnrichSeq-LOGO.png)

# Phage-EnrichSeq
Finding phages for phage therapy using Enrichment followed by Differential Genome Amplification

## Getting Started.
### Activate enrichseq conda environment with pre-installed dependencies
* To setup the ![conda](https://docs.conda.io/en/latest/miniconda.html) environment, use the following command once conda is installed locally. 
```
conda env create -f environment.yml;
```
* To activate, use the following **always** before running EnrichSeq:
```
conda activate enrichseq;
```

### Download test data
* Currently, the test genomic data is on ![git lfs](https://git-lfs.github.com/), so ensure this is working locally before running the following:
```
cd simulated_test_files/; git lfs pull; cd ../;
```

### Downloading Databases.
To ensure each module runs correctly, the following installation script needs to be ran.

* From the outermost EnrichSeq directory, run the following:
```
cd kraken_module;
bash enrichseqDB_install.sh;
cd ../;
```

*NOTE*: While the above is in the **kraken_module**, the information is needed for other modules.

### Downloading Databases.
To ensure each module runs correctly, the following installation script needs to be ran.

* From the outermost EnrichSeq directory, run the following:
```
cd kraken_module;
bash enrichseqDB_install.sh;
```

*NOTE*: While the above is in the **kraken_module**, the information is needed for other modules.


## Running EnrichSeq.
* USAGE
```
nextflow Nextflow/enrichseq.nf --read <"single" or "paired" or "long"> --fasta <path to fasta> --workdir <path to output directory> --toolpath <local global path to EnrichSeq repository> --dbdir <Path to the kracken DB> --threads <# of threads to use> --genomedir <path to the directory containing genomes> 
```

* EXAMPLE (this is an example test)
```
nextflow Nextflow/enrichseq.nf --read single --fasta simulated_test_files/simulate_genomes_2_1000000reads_illumina.fa --workdir ./outdirectory --toolpath ${PWD} --dbdir kraken_module/krakenDB --threads  4 --genomedir kraken_module/ref_genomes/
``` 