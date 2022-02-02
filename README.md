![example workflow](https://github.com/Dreycey/Phage-EnrichSeq/actions/workflows/github_actions.yml/badge.svg)
![EnrichSeq](figures/EnrichSeq-LOGO.png)

# Phage-EnrichSeq
Finding phages for phage therapy using Enrichment followed by Differential Genome Amplification.

## Getting Started.
### Activate enrichseq conda environment with pre-installed dependencies
To setup the [conda](https://docs.conda.io/en/latest/miniconda.html) environment, use the following command once conda is installed locally.
* for linux
```
conda env create -f environment_linux.yml;
```
* for mac
```
conda env create -f environment_mac.yml;
```
* To activate, use the following **always** before running EnrichSeq:
```
conda activate enrichseq;
```

## Download DB
```
python3 EnrichSeq.py db_build
```

## Running EnrichSeq.
Below are examples using single end and paired end reads. For further usage, please refer to the [EnrichSeq Wiki](https://github.com/Dreycey/Phage-EnrichSeq/wiki).

* USAGE (single end reads)
```
python3 EnrichSeq.py -1 examples/single_end_reads/simulated_test_reads_illumina.fa -o single_out_example
```

* USAGE (paired end reads)
```
python3 EnrichSeq.py -1 examples/paired_end_reads/paired_illumina_1.fa -2 examples/paired_end_reads/paired_illumina_2.fa -o paired_out
```
