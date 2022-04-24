#! usr/bin/python3
"""
DESCRIPTION:
    This script is used for comparing two different methods of detecting similarity
    between genomes: (1) Jaccard similarity index using k-mer sets, and (2) Mummer's dnadiff.
    It is also used to test optimum k-mer sizes for jaccard index. 

USAGE:
    python compare_methods.py <Phage Name line-delimited file> <Multi Fasta file> <out file>

EXAMPLE:

TODO: 
    1. Make sure the genomic mapping can be parrallized. 
        1.A - see if threads can be used for individual mapping.
    2. Look for dipps in coverage along the genome.
        2.A - perhaps plot this as well! 
"""
import os
import sys
import re
import argparse
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_path}/../../modules/genomeCompare_module")
#from PathOrganizer import PathOrganizer, PathErrors, DuplicateGenomeError
from dna import DNA

# Changing to directory of script.
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)

# GLOBALS.
PATH = os.path.dirname(os.path.abspath(__file__))


def extract_genome_paths(genomes_directory: Path) -> list:
    print(f"Extracting genome paths from: {genomes_directory}")
    genomeList = []
    for genome_fasta in os.listdir(genomes_directory):
        genomeList.append(Path(genomes_directory) / Path(genome_fasta))

    return genomeList

def run_jaccard(genomeList: list, kmerLength: int) -> float:
    '''
    Creates kmers out of the provided genomes to calculate the jaccard index between
    the sets of kmers, then outputs an adjacency matrix.
    '''
    #print(f"Calculating jaccard similarity index.")
    # TODO: change to random 2 genomes
    genome1_kmers = __create_kmers(__fasta_to_genome(genomeList[0]), kmerLength)
    genome2_kmers = __create_kmers(__fasta_to_genome(genomeList[2]), kmerLength)

    # print(os.path.basename(genomeList[0]))
    # print(os.path.basename(genomeList[2]))

    a = set(genome1_kmers)
    b = set(genome2_kmers)

    intersection = len(a.intersection(b))
    union = len(a.union(b))

    return round(intersection / union, 8)


def run_dnadiff(genomes_directory: Path):
    return None


def parse_dnadiff():
    return None


def plot_method_comparison():
    '''
    dnadiff vs jaccard scatter plot
    '''
    return None


def plot_kmer_effect(genomePair: list, maxKmerSize: int, increment: int):
    '''
    Plots the jaccard similarity between two genomes only when
    the kmer size is changed. 
    '''
    plotting_dict = {'K-mer Size': [], 'Jaccard Index': []}
    for k in range(0, maxKmerSize, increment):
        if k != 0:
            plotting_dict['K-mer Size'].append(k)
            plotting_dict['Jaccard Index'].append(run_jaccard(genomePair, k))
    
    kmer_df = pd.DataFrame.from_dict(plotting_dict)
    sns.lineplot(data=kmer_df, x="K-mer Size", y="Jaccard Index", markers=True)
    filename = 'jaccard_vs_kmer-size.png'
    plt.savefig(filename)
    print(f'Line plot stored in {filename}')
    


# PRIVATE METHODS
def __create_kmers(genome, kmerLength = 20) -> list:
    ''' Generates k-mers given a dna sequence and specified k-mer length'''
    kmers = []
    genomeLength = len(genome)
    # assigns 1 kmer if the length passed is bigger than genome length
    if genomeLength < kmerLength or kmerLength <= 0:
        numKmers = 1
        kmerLength = genomeLength
    else:
        numKmers = genomeLength - kmerLength + 1

    for i in range(numKmers):
        kmer = genome[i:i + kmerLength]
        kmers.append(kmer)

    return kmers

def __fasta_to_genome(fasta_path) -> str:
    ''' 
        DESCRIPTION:
            Extracts just the genome from the genome member variable, which should be a file path.
            Currently assumes only single FASTA file is passed (no multi-fasta). 
        
        INPUT:
            FASTA/FASTQ file (.fa, .fna, .fasta, .fastq)
        
        OUTPUT:
            Genome in string format
    '''
    sequence = ""
    with open(fasta_path) as fasta_file:
        fasta_lines = fasta_file.readlines()
        sequence_i = ""
        # If name is found (i.e. this is the correct file), get genome
        for counter, line in enumerate(fasta_lines):
            if line[0] == ">": #and re.search(str(self.taxid), line, re.IGNORECASE):
                if (counter != 0):
                    sequence += (sequence_i.rstrip())
                    sequence_i = ""
            elif line[0] != ">":
                sequence += line.rstrip()
    return sequence

def __create_kmer_list(kmerMax: int, increments: int) -> list:
    kmerSizeList = []

    return kmerSizeList


def parseArgs(argv=None) -> argparse.Namespace:
    '''
    DESCRIPTION:
        This method takes in the arguments from the command and performs
        parsing.

    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-d", "--genome_directory", help="path to the directory of genomes", required=True)
    #parser.add_argument("-o", "--output_dir", help="the output directory", required=False)
    parser.add_argument("-k1","--max_kmer_size", help="size (length) of kmers to create, or the max size", required=False)
    parser.add_argument("-k2", "--kmer_increments", help="amount to increment kmers by for kmer benchmarking", required=False)
    return parser.parse_args(argv)

def main():
    arguments = parseArgs(argv=sys.argv[1:])
    genomeList = extract_genome_paths(arguments.genome_directory)
    plot_kmer_effect(genomeList, int(arguments.max_kmer_size), int(arguments.kmer_increments))


if __name__ == "__main__":
    main()