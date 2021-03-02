##Creating BlastDB

# Usage

## Build the blast database

* update the paths in the bash script (TODO: allow for input)

```
bash blastBuild.sh;
```

## Query using blastN

* update the paths in the bash script (TODO: allow for input)

```
bash blastRun.sh;
```

Currently this runs a test on a simulated reads file with the following 
parameters:

```
# <genome path> <relative abundance>
../genomes/Blessica.fasta, 0.25   
../genomes/D29_genome.fa, 0.25
../genomes/perseus_genome.fa, 0.5
```

