
# Usage

## Build the bracken database

* insert required parameters in bracken_module.config 

```
bash brackenBuild.sh;
```

## Run bracken for relative abundance

```                                                                             
 bash brackenRun.sh --krakendb=kraken_module/<db directory> --input=<.kraken outfile> --out=<outfile prefix> --read=<read length>
```

## Example Run
* for example:
```
 bash brackenRun.sh --krakendb=kraken_module/krakenDB --input=kraken_module/kraken_out_20210303_1615.kraken --out=bracken_out_01 --read=100
```

Currently this runs a test on a simulated reads file with the following 
parameters:

```
# <genome path> <relative abundance>
../genomes/Blessica.fasta, 0.25   
../genomes/D29_genome.fa, 0.25
../genomes/perseus_genome.fa, 0.5
```

