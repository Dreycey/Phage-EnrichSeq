#! usr/bin/python3
"""
DESCRIPTION:
    This script is used for comparing two different methods of detecting similarity
    between genomes: (1) Jaccard similarity index using k-mer sets, and (2) Mummer's dnadiff.
    It is also used to test optimum k-mer sizes for jaccard index. 

USAGE:
    python compare_methods.py -d <path to reference genomes> -k1 maximum k-mer size -k2 increments
"""
import os
import sys
import re
import argparse
import subprocess
from pathlib import Path
import pandas as pd
import seaborn as sns
import altair as alt
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(f"{current_path}/../../modules/genomeCompare_module")

# Changing to directory of script.
#abspath = os.path.abspath(sys.argv[0])
#dname = os.path.dirname(abspath)
#os.chdir(dname)

# GLOBALS.
#PATH = os.path.dirname(os.path.abspath(__file__))


# def extract_genome_paths(genomes_directory: Path) -> list:
#     print(f"Extracting genome paths from: {genomes_directory}")
#     genomeList = []
#     for genome_fasta in os.listdir(genomes_directory):
#         genomeList.append(Path(genomes_directory) / Path(genome_fasta))

#     return genomeList


def extract_genome_path(genomes_directory: Path, filename: str) -> list:
    print(f"Extracting genome paths from: {genomes_directory}")

    return Path(genomes_directory) / Path(filename)


def populate_genome_dict(genomes_directory: Path) -> dict:
    '''
    Dictionary: [ fullpath : genome str]
    '''
    print(f'Populating genome dictionary from: {genomes_directory}')
    genomeDict = {}
    for genome_filename in os.listdir(genomes_directory):
        fullpath = Path(genomes_directory) / Path(genome_filename)
        genomeDict[fullpath] = __fasta_to_genome(fullpath)

    return genomeDict


def run_jaccard(genomePair: list, kmerLength=35) -> float:
    '''
    Creates kmers out of the provided genomes to calculate the jaccard index between
    the sets of kmers

    INPUT:
        genomePair: list of (2) genomes
    '''
    kmerLength=int(kmerLength)
    genome1_kmers = __create_kmers(__fasta_to_genome(genomePair[0]), kmerLength)
    genome2_kmers = __create_kmers(__fasta_to_genome(genomePair[1]), kmerLength)
    #genome1_kmers = __create_kmers(genomePair[0], kmerLength)
    #genome2_kmers = __create_kmers(genomePair[1], kmerLength)
    
    a = set(genome1_kmers)
    b = set(genome2_kmers)

    intersection = len(a.intersection(b))
    union = len(a.union(b))
    
    return round((intersection / union)*100, 4)


def run_dnadiff(genomePair: list, outputdir: str):
    print('Running dnadiff')
    #subprocess.run(['python', commandPath, '-q', genomeList[0], '-r', genomeList[2]]) 
    subprocess.run(['dnadiff', '-p', Path(outputdir)/Path('dnadiff_out'), genomePair[0], genomePair[1]])
    return Path(outputdir + '/dnadiff_out.report')


def parse_dnadiff(dnadiff_outdir: str) -> float:
    similarity=-1.0
    with open(dnadiff_outdir) as dnadiff_file:
        for line in dnadiff_file:
            if 'AvgIdentity' in line:
                #re.search(str(self.taxid), line, re.IGNORECASE).group(1)
                #val1 = line.rstrip().split()[1]
                #val2 = line.rstrip().split()[2]
                #num1 = float(val1[val1.find("(")+1:val1.find(")")].strip('%'))
                #num2 = float(val2[val2.find("(")+1:val2.find(")")].strip('%'))
                #similarity = (num1 + num2) / 2
                similarity = line.rstrip().split()[1]
                break; # only first occurrence
    return similarity



def plot_method_comparison(genomeList: list):
    '''
    dnadiff vs jaccard scatter plot
    '''
    pass


def plot_simulated_percentages(genomes_directory: Path, original_filename: str, kmer: int, increment: int, outputDir: str):
    '''
    Plots jaccard index vs. simulated percentages.
    
    INPUT:
        genomeDict: dictionary of genomes to conduct pairwise comparisons on, including the original
            [<filename> : <genome str>]
        original_filename: the filename (i.e. dictionary key) of the original, or reference genome to compare with
        kmer: kmer length for jaccard sets
        outputDir: location to store plot
    '''
    print('PLOTTING JACCARD VS. SIMULATED')
    plotting_dict = {'Simulated Percentage': [], 'Jaccard Index': [], 'K-mer Length': []}
    genomeDict = populate_genome_dict(genomes_directory)
   
    genomePair = [ extract_genome_path(genomes_directory, original_filename), None]

    for k in range(0,int(kmer),int(increment)):
        for key in genomeDict.keys():
            percent_str = re.sub('\D', '', str(key))
            genomePair[1] = key
            plotting_dict['Simulated Percentage'].append(float(percent_str))
            jaccard = run_jaccard(genomePair, k)
            plotting_dict['Jaccard Index'].append(jaccard)
            plotting_dict['K-mer Length'].append(k)

    simPercent_df = pd.DataFrame.from_dict(plotting_dict)
    print(simPercent_df)
           
    #ALTAIR PLOT
    chart = alt.Chart(simPercent_df).mark_line(point=True).encode(
        x='Simulated Percentage:Q',
        y='Jaccard Index:Q',
        color='K-mer Length:N'
    )       
    filename = outputDir + 'jaccard_vs_simulated_' + str(kmer) + '-mer.png'
    chart.save(filename)
    print(f'Line plot stored in {filename}')   



def plot_kmer_effect(genomePair: list, maxKmerSize=30, increment=1):
    '''
    Plots the jaccard similarity between two genomes only when
    the kmer size is changed. 
    '''
    plotting_dict = {'K-mer Size': [], 'Jaccard Index': []}
    for k in range(0, maxKmerSize, increment):
        #if k != 0:
            plotting_dict['K-mer Size'].append(k)
            plotting_dict['Jaccard Index'].append(run_jaccard(genomePair, k))
    
    kmer_df = pd.DataFrame.from_dict(plotting_dict)
    
    ## SEABORN PLOT
    # fig, ax = plt.subplots()
    # sns.lineplot(data=kmer_df, x="K-mer Size", y="Jaccard Index", marker='o')
    # ax.set_xlim(1, maxKmerSize)
    # ax.set_ylim(0.0, 1.1)
    # filename = 'jaccard_vs_kmer-size.png'
    # plt.savefig(filename)
    # print(f'Line plot stored in {filename}')

    ## ALTAIR PLOT
    chart = alt.Chart(kmer_df).mark_line(
        point=True
    ).encode(
        x='K-mer Size:Q',
        #alt.X('K-Mer Size:Q',scale=alt.Scale(domain=(0, maxKmerSize))),
        y='Jaccard Index:Q'
        #alt.Y('Jaccard Index:Q', scale=alt.Scale(domain=(0, 1.1),clamp=True))
    ).properties(width=500)
    filename = 'jaccard_vs_kmer-size.png'
    chart.save(filename)
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
    print(f'File: {fasta_path}')
    with open(fasta_path) as fasta_file:
        fasta_lines = fasta_file.readlines()
        sequence_i = ""
        for counter, line in enumerate(fasta_lines):
            if line[0] == ">": 
                if (counter != 0):
                    sequence += (sequence_i.rstrip())
                    sequence_i = ""
            elif line[0] != ">":
                sequence += line.rstrip()
    print(sequence[:21])
    return sequence


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
    parser.add_argument("-d", "--genome_directory", help="path to the directory of genomes", required=True)
    parser.add_argument("-o", "--output_dir", help="the output directory", required=True)
    parser.add_argument("-k1","--max_kmer_size", help="size (length) of kmers to create, or the max size", required=False)
    parser.add_argument("-k2", "--kmer_increments", help="amount to increment kmers by for kmer benchmarking", required=False)
    return parser.parse_args(argv)

def main():
    arguments = parseArgs(argv=sys.argv[1:])
    #genomePair = extract_genome_paths(arguments.genome_directory)
    #plot_kmer_effect(genomeList, int(arguments.max_kmer_size), int(arguments.kmer_increments))
    #print(run_dnadiff(genomeList, arguments.output_dir))
    #dnadiff_val = parse_dnadiff(run_dnadiff(genomeList, arguments.output_dir))
    #print(f'dnadiff = {dnadiff_val}')
    #print(f'jaccard = {run_jaccard(genomeList, int(arguments.max_kmer_size))}')
    #genomeDict = populate_genome_dict(arguments.genome_directory)
    plot_simulated_percentages(arguments.genome_directory, 'genome_100.fa', arguments.max_kmer_size,arguments.kmer_increments, arguments.output_dir)

if __name__ == "__main__":
    main()
