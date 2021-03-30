Required input: sequencing file (.fa/.fna/.fasta), readtype (single/paired)
[done] Process 1: genome assembly using megaHIT
Process 2: genome assembly using MetaSPAdes
Process 3: assembly comparison/consolidation
*Process 4: classification using Kraken2
Process 5: classification using MetaPhlAn2
Process 6: binning comparison/consolidation
*Process 7: relative abundance using Bracken
*Process 8: make pretty output --> make into several processes later depending on type of output requested

Before running, needs to know:
- location of the simulated fasta files or those done in lab -can be passed as a parameter for nextflow command
- location of databases --> should there be a process to build databases?

## TODO
1. Add database location as option in krakenBuild script
2. Add paired end options for megahitRun script
3. Add paired end options for Nextflow 
4. Conda environment
