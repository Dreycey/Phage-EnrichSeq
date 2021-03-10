## Run kraken2 for phage classification

# Usage

## Update kraken_module.config file

local genomeDir="ref_genomes"; (where do you want to keep the downloaded reference genomes?)
local dbDir="krakenDB"; (name of the database directory you want created)
local actinoOutFile="actinoReformatted.fa"; (new Actinobacteriophage)

## Build the kraken2 database

```
bash kraken2Build.sh;
```

## Run kraken2

Syntax
```
bash kraken2Run.sh --krakendb=<path to kraken DB> --queryfasta=< input fasta file> --report=< kraken report file output > --out=< kraken output log >;
```

Example
```
bash kraken2Run.sh --krakendb=krakenDB/ --queryfasta=inputfasta/simulatedgenomes_illumina.fa --report=kraken.report.txt --out=krakenout.kraken;
```


