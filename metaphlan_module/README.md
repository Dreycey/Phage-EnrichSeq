
# Usage

## Build the Metaphlan database

* insert required parameters in bracken_module.config

```
bash metaphlanBuild.sh;
```

## Run metaphlan for classification and relative abundance

```                                                                             
 bash metaphlanRun.sh --input=<input sequence> --type=<input file type> --out=<output file> --nproc <number of cores>
```

## Example Run
* for example:
```
 bash metaphlanRun.sh --input=assembly/megahitout.fa --type=fasta --out=profiled_metagenome.txt --nproc 4
```

Currently this runs a test on a simulated reads file with the following
parameters:

```
# <genome path> <relative abundance>
../genomes/Blessica.fasta, 0.25   
../genomes/D29_genome.fa, 0.25
../genomes/perseus_genome.fa, 0.5
```
