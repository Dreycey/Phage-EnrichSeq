# Phage-EnrichSeq
Finding phages for phage therapy using Enrichment followed by Differential Genome Amplification


## Description
This repository contains information and code regarding pipeline for finding phages capable of infecting bacteria of interest.


## Major Goals

### M.G. 1 - create a RNA-Seq pipeline. (another option would be creating tools from scratch)
We want to download RNA seq data and create a pipeline (pipeline = RNA-seq tools combined to perform analysis in a specific order, so like DAG that uses the output of one tool as input for another). The pipeline will be written in the language Snakemake (https://snakemake.readthedocs.io/en/stable/). A good example of such a pipeline is here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3433-x.

### M.G. 2 - Create theoretical data yeilded from experiments.
We will simulate DNA sequencing of phage genomes before and after differential genome amplification. What does this mean? This means we will download different phage genomes. We will then simulate reads for each genome using an illumina read simulator (ART, InSilicoSeq, or SimuScop).

InSilicoSeq: https://academic.oup.com/bioinformatics/article/35/3/521/5055123
SimuScop: https://link.springer.com/article/10.1186/s12859-020-03665-5
ART: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm

### M.G. 3 - Convert pipeline from MG 1 to handle theoretical input from MG 2
So what is this? Basically we created a software tool (not our own software, but strung them together) in MG 1. In MG 2 we created data that we would expect from the sequencing experiments (if they were to happen, and were work as we expect), and we will add noise to some the simulated data sets (again, MG 2). 

NOW we want to transform our pipeline to handle the simulated data. And because simulated, we know the ground truth. This allows to test the accuracy of the pipeline for the given task at hand. What's the given task? We want to find specific genomes that have been amplified from enrichment using sequencing data. 







 


