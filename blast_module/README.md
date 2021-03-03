# BlastN module

## Usage

### Build the blast database

* update the paths in the `blast_module.config` file before running the following:

```
bash blastBuild.sh;
```

### Query using blastN

* pass the correct, requiired, args to the run script
```                                                                             
 bash blastRun.sh --blastdb=blastdb/<db file name> --queryfasta=<fasta input> --out=<outfile name>
```

## Example Run
* for example:
```
 bash blastRun.sh --blastdb=blastdb/outputMulti3.fa --queryfasta=blast_testfiles/simulatedgenomes_illumina.fa --out=blastout.txt
```

Currently this runs a test on a simulated reads file with the following 
parameters:

```
# <genome path> <relative abundance>
../genomes/Blessica.fasta, 0.25   
../genomes/D29_genome.fa, 0.25
../genomes/perseus_genome.fa, 0.5
```

