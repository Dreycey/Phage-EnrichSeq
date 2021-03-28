# This script has 3 main jobs:
# 1) parses output files from abundance estimation tools
# 2) analyzes the data from those files
# 3) displays relevant information in a readable manner
import sys
import csv
from re import search
from tqdm import tqdm


abundance_dict = {}
blast_dict = {}
abundance_list = []

def parseBrackenFile(bracken_file):
    bracken_contents = open(bracken_file).readlines()[1:]
    for line in tqdm(bracken_contents):
        line = line.strip("\n").split()

        species_name = line[2].lower()
        abundance = line[8]
        abundance_dict[species_name] = [abundance]

def parseBlastFile(blast_file):
    ## blast output columns:
    # qseqid: query (e.g., unknown gene) sequence id
    # sseqid: subject (e.g., reference genome) sequence id
    # sblastnames: Subject Blast Name(s), separated by a ';'   (in alphabetical order)
    # score: Raw score
    # evalue: number of expected hits of similar quality (score), small e-value - better match
    # pident: percentage of identical matches
    # length: alignment length (sequence overlap)
    # mismatch: number of mismatches
    # gapopen: number of gap openings
    # gaps: Total number of gaps
    # sseq: Aligned part of subject sequence
    blast_contents = open(blast_file).readlines()
    #blast_dict

#def parseBowtieFile(bowtie_file):


def saveInfoToFile(outfile):
    output = open(outfile, "w")
    # save names, taxid and abundances of dict contents to file
    if '.csv' in outfile:
        with open(outfile, mode='w') as csv_file:
            fieldnames = ['phage_name', 'bracken_abundance', 'blast_abundance']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()

            for key in abundance_dict:
                writer.writerow({'phage_name': key, 'bracken_abundance': abundance_dict[key][0]})


def main():
    """ controls the script """

    ### SCRIPT INPUT
    bracken_outfile = sys.argv[1]
    blast_outfile = sys.argv[2]
    outfile = sys.argv[3]

    #blast_contents = open(blast_outfile).readlines()

    # TEST BRACKEN PARSE AND DISPLAY
    parseBrackenFile(bracken_outfile)
    for key in abundance_dict:
        print(f'phage: {key} | abundance: {abundance_dict[key]}')
    saveInfoToFile(outfile)

if __name__ == "__main__":
    main()



