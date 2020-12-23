# metagenomics relative abundance module
The goal of this module is run common place software used
for metagenomic binning and relative abundance estimation.
The output from this module can be used to gather information
on how the relative abundance of specific phages changes 
before and after enrichment. 

# USAGE

```

```

## OUTLINE

This pipeline works as the following:
```
          input fasta (simulated or real)
                     |
                     |------------------------ *Assembly module* 
                     |                                 |             
                     |                                / \  
                     |                      metaSPADES   MegaHIT
                     |                               \   /
                     |                           assembleCompare
                     |<_______________________________ |
                     |
                     |
              *Abundance Module*
                    / \
               ____/   \_____                    
               |             |    
               |             |
           kracken2       MetaPhlAn2
                \           /
              abundanceCompare
                     |
                     |
                OUTPUT FILES

```

OUTPUT FILES:
1. figs (abundance changes, kronos)
2. A latex report and pdf based on information from the analysis.
3. An HTML output report.

ADD EVENTUALLY:
1. Read mapping (where  in the genomes are the reads and assemblies from?)
2. SV detection/analysis
3. abundanceCompare: Math/stats to choose confidence values for the relative abundance estimations.
4. abundanceCompare: Math/stats to find fold change throughout the enrichment process

