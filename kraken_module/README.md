## Run kraken2 for phage classification

# Usage

## Build the kraken2 database

* update the paths in the bash script (TODO: allow for input)

```
bash kraken2Build.sh;
```

## Run kraken2

* update the paths in the bash script (TODO: allow for input)

```
bash kraken2Run.sh;
```

Currently this runs a test on a simulated reads file with the following 
parameters:

```
# <genome path> <relative abundance>
../genomes/Blessica.fasta, 0.25   
../genomes/D29_genome.fa, 0.25
../genomes/perseus_genome.fa, 0.5
```

