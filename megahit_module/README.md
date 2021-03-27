## Run MegaHIT for genome assembly


* TODO: Change flags to --read for single and --pair1 and --pair2 for paired reads

# Usage

## Run megahit
### Syntax

```
bash megahitRun.sh --read=<single/paired/long> --input=<input fasta file> --threads=<number of threads> --out=<output directory>
```

### Example
```
bash megahitRun.sh --read=single --input=inputfasta/simgenomes.fa --threads=4 --out=megahit_20210307
```
