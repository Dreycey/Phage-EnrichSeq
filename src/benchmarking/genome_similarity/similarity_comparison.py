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
import time
import random
from pathlib import Path
from itertools import combinations
from venv import create
import pandas as pd
import seaborn as sns
import altair as alt
import matplotlib.pyplot as plt

current_path = os.path.dirname(os.path.abspath(__file__))


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
        for counter, line in enumerate(fasta_lines):
            if line[0] == ">": 
                if (counter != 0):
                    sequence += (sequence_i.rstrip())
                    sequence_i = ""
            elif line[0] != ">":
                sequence += line.rstrip()
    return sequence

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
        if genome_filename.endswith('fa') or genome_filename.endswith('fna') or genome_filename.endswith('fasta'):
            fullpath = Path(genomes_directory) / Path(genome_filename)
            genomeDict[fullpath] = __fasta_to_genome(fullpath)
            if len(genomeDict) == 35:
                break

    return genomeDict

def find_pair(genomeDict: dict, genomeA: str, genomeB: str):
    genomePair = []
    for genPath in genomeDict.keys():
        if genomeA in str(genPath) or genomeB in str(genPath):
            genomePair.append(genPath)
    print(genomePair)
    return genomePair


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
    
    a = set(genome1_kmers)
    b = set(genome2_kmers)

    intersection = len(a.intersection(b))
    union = len(a.union(b))
    
    return round((intersection / union)*100, 4)


def run_dnadiff(genomePair: list, outputDir: str):
    print('Running dnadiff')
    subprocess.run(['dnadiff', '-p', Path(outputDir)/Path('dnadiff_out'), genomePair[0], genomePair[1]])
    #return Path(outputdir + '/dnadiff_out.report')
    similarity=-1.0
    with open(outputDir+ '/dnadiff_out.report') as dnadiff_file:
        for line in dnadiff_file:
            if 'AvgIdentity' in line:
                similarity = line.rstrip().split()[1]
                break; # only first occurrence
    return similarity


# def parse_dnadiff(dnadiff_outdir: str) -> float:
#     similarity=-1.0
#     with open(dnadiff_outdir) as dnadiff_file:
#         for line in dnadiff_file:
#             if 'AvgIdentity' in line:
#                 similarity = line.rstrip().split()[1]
#                 break; # only first occurrence
#     return similarity

def create_genome_pairs(genomeDict: dict):
    '''
    Creates specific genome pairs for simulated similarity. First genome is always 
    the unchanged FASTA file (original), and the second is the simulated FASTA. 
    '''
    genomePairs = []
    #count = 0
    for key in genomeDict.keys():
        pair = [None, None]
        if int(re.sub('\D', '', str(key.name))) is not 100:
            #count += 1
            pair[1] = key
            genomePairs.append(pair)
        else:
            file100 = key

    for p in genomePairs:
        p[0] = file100

    return genomePairs


def create_random_pairs(genomeDict: dict):
    '''
    Creates and returns a random list of genome pairs given a dictionary of genomes. 
    '''
    done = []
    genomePairs = []

    return genomePairs


def plot_method_comparison(genomeDict: dict, kmerLength: int, outputDir: str):
    '''
    dnadiff vs jaccard scatter plot

    INPUT:
        genomeDict: dictionary of genomes (must be pre-populated)
        kmerLength: k-mer length for jaccard calculation based on k-mer set comparisons
        outputDir: where to save dnadiff files and plots

    OUTPUT:
        None. Show or save plots.

    '''
    print('Comparing dnadiff vs. Jaccard index')
    compareDict = {'genome pair': [], 'jaccard value': [], 'dnadiff value': [], 'jaccard time': [], 'dnadiff time':[], 'speed factor':[]}
    
    # create pairs from dictionary
    #genomePairs = list(combinations(genomeDict.keys(), 2))
    genomePairs = create_genome_pairs(genomeDict)

    # run jaccard and dnadiff for all pairs
    for num, pair in enumerate(genomePairs):
        compareDict['genome pair'].append(num+1)

        # run dnadiff
        startTime = time.time()
        dnadiffVal = run_dnadiff(pair, outputDir)
        dnadiffTime = time.time() - startTime
        compareDict['dnadiff value'].append(dnadiffVal)

        # run jaccard
        startTime = time.time()
        jaccardVal = run_jaccard(pair, kmerLength)
        jaccardTime = time.time() - startTime
        compareDict['jaccard value'].append(jaccardVal)


        # append runtimes and calculate speed factor
        compareDict['jaccard time'].append(jaccardTime)
        compareDict['dnadiff time'].append(dnadiffTime)
        compareDict['speed factor'].append(1/(jaccardTime/dnadiffTime))

    # Display table
    compare_df = pd.DataFrame.from_dict(compareDict)
    print(compare_df)

    # Plot dnadiff vs. jaccard vals
    chart = alt.Chart(compare_df).mark_circle(size=70).encode(
        alt.X('jaccard value:Q',
            scale=alt.Scale(domain=[0, 100])
        ),
        alt.Y('dnadiff value:Q',
            scale=alt.Scale(domain=[0, 100])),
        color='genome pair:N',
        tooltip=['jaccard value', 'dnadiff value']
    ).show()
    # filename = outputDir + 'dnadiff_vs_dnadiff_' + str(kmerLength) + '-mer.png'
    # chart.save(filename)
    # print(f'Scatterplot stored in {filename}') 


def plot_runtime(genomePair: list, kmers: list, outputDir: str):

#   # TODO calculate average
    plotting_dict = {'Method': [], 'k-mer length': [], 'Runtime(s)': []}
    for k in kmers:
        plotting_dict['k-mer length'].append(k)
        # Run and time jaccard
        plotting_dict['Method'].append('Jaccard')
        startTime = time.time()
        run_jaccard(genomePair, k)
        jaccardTime = time.time() - startTime
        plotting_dict['Runtime(s)'].append(jaccardTime)

        plotting_dict['k-mer length'].append(k)
        # Run and time dnadiff
        plotting_dict['Method'].append('dnadiff')
        startTime = time.time()
        run_dnadiff(genomePair, outputDir)
        dnadiffTime = time.time() - startTime
        plotting_dict['Runtime(s)'].append(dnadiffTime)

    # PLOTTING
    runtime_df = pd.DataFrame.from_dict(plotting_dict)
    print(runtime_df)
    alt.Chart(runtime_df).mark_bar().encode(
        x='k-mer length:O',
        y='Runtime(s):Q',
        color='Method:N'
).properties(width=400).show()

def plot_comparison_with_kmers(genomeDict: dict, kmers: list, outputDir: str):
    '''
    Dnadiff vs Jaccard plot for multiple k-mer values 
    '''
    plotting_dict = {'kmer length': [], 'jaccard value': [], 'dnadiff value': []}
    genomePairs = create_genome_pairs(genomeDict)
    # Calculations
    for k in kmers:
        for pair in genomePairs:
            plotting_dict['kmer length'].append(k)
            plotting_dict['jaccard value'].append(run_jaccard(pair, k))
            plotting_dict['dnadiff value'].append(run_dnadiff(pair, outputDir))

    # Plotting
    compare_df = pd.DataFrame.from_dict(plotting_dict)
    chart = alt.Chart(compare_df).mark_circle(size=50).encode(
        alt.X('jaccard value:Q',
            scale=alt.Scale(domain=[0, 100])
        ),
        alt.Y('dnadiff value:Q',
            scale=alt.Scale(domain=[0, 100])),
        color='kmer length:N',
        tooltip=['jaccard value', 'dnadiff value']
    ).show()

def plot_simulated_percentages(genomes_directory: Path, original_filename: str, kmer_min: int, kmer_max: int, increment: int, outputDir: str):
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
    genomeDict = populate_genome_dict(genomes_directory) # All fasta files in a directory
   
    genomePair = [ extract_genome_path(genomes_directory, original_filename), None]

    for key in genomeDict.keys():
        for k in range(int(kmer_min),int(kmer_max)+1,int(increment)):
            percent_str = re.sub('\D', '', str(key.name)) # remove anything  that isn't a number
            genomePair[1] = key
            plotting_dict['Simulated Percentage'].append(float(percent_str))
            jaccard = run_jaccard(genomePair, k)
            plotting_dict['Jaccard Index'].append(jaccard)
            plotting_dict['K-mer Length'].append(k)

    simPercent_df = pd.DataFrame.from_dict(plotting_dict)
    print(simPercent_df)
           
    #ALTAIR PLOT
    chart = alt.Chart(simPercent_df).mark_line(point=True).encode(
        alt.X('Simulated Percentage:Q',
            scale=alt.Scale(zero=True)
        ),
        y='Jaccard Index:Q',
        color='K-mer Length:N'
    ).properties(width=700).show()      
    # filename = outputDir + 'jaccard_vs_simulated_' + str(kmer_max) + '-mer.png'
    # chart.save(filename)
    # print(f'Line plot stored in {filename}')   


def plot_dnadiff_vs_simulated(genomeDict: dict, outputDir: str):
    '''
    Plots dnadiff ANI value vs. simulated ANI.
    '''
    print('PLOTTING DNADIFF VS. SIMULATED')

    plotting_dict = {'Simulated ANI': [], 'dnadiff ANI': []}
    genomePair = [None, None]

    # First genome in pair is always the one with 100 in the name
    for key in genomeDict.keys():
        if int(re.sub('\D', '', str(key.name))) is 100:
            genomePair[0] = key

    for key in genomeDict.keys():
        # remove anything  that isn't a number
        simANI = re.sub('\D', '', str(key.name))
        genomePair[1] = key
            
        plotting_dict['Simulated ANI'].append(simANI)
        dnadiffVal = run_dnadiff(genomePair, outputDir)
        plotting_dict['dnadiff ANI'].append(dnadiffVal)

    # Display table
    dnadiff_df = pd.DataFrame.from_dict(plotting_dict).sort_values(by='Simulated ANI')
    print(dnadiff_df)

    # Plot dnadiff vs. jaccard vals
    chart = alt.Chart(dnadiff_df).mark_circle(size=70).encode(
        alt.X('Simulated ANI:Q',
            scale=alt.Scale(zero=True)
        ),
        alt.Y('dnadiff ANI:Q',
            scale=alt.Scale(zero=True)),
        tooltip=['Simulated ANI', 'dnadiff ANI']
    ).properties(width=650).show()


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

    ## ALTAIR PLOT
    chart = alt.Chart(kmer_df).mark_line(
        point=True
    ).encode(
        #x='K-mer Size:Q',
        alt.X('K-Mer Size:Q',scale=alt.Scale(zero=True)),
        y='Jaccard Index:Q'
        #alt.Y('Jaccard Index:Q', scale=alt.Scale(domain=(0, 1.1),clamp=True))
    ).properties(width=500)
    filename = 'jaccard_vs_kmer-size.png'
    chart.save(filename)
    print(f'Line plot stored in {filename}')
    


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
    return parser.parse_args(argv)

def main():
    arguments = parseArgs(argv=sys.argv[1:])
    #plot_simulated_percentages(arguments.genome_directory, 'genome_100.fa', 7, 10, 1, arguments.output_dir)
    genomeDict = populate_genome_dict(arguments.genome_directory)
    plot_simulated_percentages(arguments.genome_directory, 'genA_100ANI.fa', 5, 12, 1, arguments.output_dir)
    #plot_method_comparison(genomeDict, 6, arguments.output_dir)
    #plot_dnadiff_vs_simulated(genomeDict, arguments.output_dir)
    #plot_comparison_with_kmers(genomeDict, [6,7,8,9,10], arguments.output_dir)
    #plot_runtime(genomePair=find_pair(genomeDict, 'Blessica', 'D29'),kmers=[6,7,8,9,10],outputDir=arguments.output_dir)

if __name__ == "__main__":
    main()
