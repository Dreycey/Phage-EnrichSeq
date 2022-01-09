import sys
import os
import csv
import numpy as np
import pandas as pd
import argparse
from typing import List
from dna import DNA
from pathlib import Path
current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_path}/../PathOrganizer_module")
from PathOrganizer import PathOrganizer, PathErrors, DuplicateGenomeError

''' A class that compares between a list of genomes using an adjacency matrix
    of jaccard similarity indices '''


class GenomeCompare:
    def __init__(self, inputfile, genome_dir, kmer_length = 20):
        self.dnaList: List = self.parse_input_file(inputfile, genome_dir, int(kmer_length))
        self.adjacencyMatrix = np.zeros((len(self.dnaList), len(self.dnaList)), float)
        self.clusters: dict = {}
        self.abundances: dict = {}


    def parse_input_file(self, inputfile, genome_dir, kmer_length):
        '''
            DESCRIPTION:
                Parses a line-separated file to extract taxon IDs and creates DNA objects using them
            INPUT:
                line-separated file containing taxon IDs
            OUTPUT:
                List of DNA objects
        '''
        dnaList = []
        print(f"Extracting taxids from LSV file: {inputfile}")
        taxidfile = open(inputfile, 'r')
        #with open(inputfile) as taxidfile:
        for taxid in taxidfile:
            #taxid = taxidfile.readlines()
            taxid = taxid.rstrip()
            print(taxid)
            dnaList.append(DNA(taxid, genome_dir, kmer_length))

        return dnaList

    
    def create_adjacency_matrix(self):
        ''' 
            DESCRIPTION:
                Creates and returns a 2D array representing similarities of DNA objects based
                jaccard index
            INPUT:
                -
            OUTPUT:
                2D array of float values (jaccard indices)
        '''
        print("Creating adjacency matrix from found DNA...")

        # update matrix values
        for i in range(len(self.dnaList)):
            for j in range(len(self.dnaList)):
                self.adjacencyMatrix[i,j] = self.dnaList[i].calc_jaccard(self.dnaList[j])

        return self.adjacencyMatrix

    def prune_adj_matrix(self, threshold):
        ''' 
            DESCRIPTION:
                Retrieves DNA objects that are genetically similar based on some threshold
            INPUT:
                Threshold value between 0.0 and 1.0 (default=0.0)
            OUTPUT:
                A dictionary of all clusters (key: cluster name, value: list of dna objects in cluster). 
                Clusters are represented by the set data structure. Can have a cluster of size=1.
        '''
        print(f"Generating clusters based on threshold value {threshold}")
        threshold = float(threshold)
        clusters = {}
        cluster_num = 0
        added_dna = set() # keep track of already clustered dna objects
        for i in range(len(self.dnaList)):  
            if self.dnaList[i] not in added_dna:
                cluster = set() # initialize a new cluster
                cluster_num += 1
                cluster_name = "C_" + str(cluster_num)
                cluster.add(self.dnaList[i]) # add current dna obj to the cluster set
                added_dna.add(self.dnaList[i]) 
            
                for j in range(i,len(self.dnaList)):
                    if self.adjacencyMatrix[i,j] > threshold and self.dnaList[j] not in added_dna:
                        cluster.add(self.dnaList[j])
                        added_dna.add(self.dnaList[j])
                
                clusters[cluster_name] = cluster # add the cluster(s) to the dictionary
        
        self.clusters = clusters
        return clusters
    

    def estimate_abundances(self, inputfile):
        '''
        DESCRIPTION:


        INPUT:
            merge_overlap_filter/merge_overlap_out_refined.csv

        OUTPUT:
            dictionary (key: cluster name, value: relative abundance)

        TODO: Improve time complexity
        '''
        print("Calculating cluster abundances...")
        self.abundances = {}
        with open(inputfile, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                taxid = int(row[0]) if row[0] != 'UNK' else str(row[0])
                abundance = float(row[1])
                if taxid == 'UNK':
                    self.abundances[taxid] = abundance
                    continue
                for clusterName,dnaSet in self.clusters.items():
                    for dnaObj in dnaSet:
                        if dnaObj.taxid == taxid:
                            if clusterName in self.abundances.keys():
                                self.abundances[clusterName] += abundance
                            else:
                                self.abundances[clusterName] = abundance
                       

    
    def display_adjacency_matrix(self):
        # TODO: make pretty
        # TODO: raise exception if empty?
        print(self.adjacencyMatrix)


    def output_to_file(self, file_path, isDNA=True):
        '''
        DESCRIPTION:
            Prints cluster information to a csv file
        INPUT:
            Prefix for file to create
        OUTPUT:
            CSV file name
            CSV file example:
                C_1, taxid1,
                C_1, taxid3
                C_2, taxid2
                C_3, taxid4
                C_3, taxid5
        '''
        input_dict = self.clusters if isDNA else self.abundances
        # loop through object and save to CSV
        file_out = file_path+"clusters.csv" if isDNA else file_path+"abundances.csv"
        with open(file_out, 'w') as csvfile:
            writer = csv.writer(csvfile)
            for key, value in input_dict.items():
                if isDNA:
                    for dna in value:
                        writer.writerow([key, dna.taxid])
                else:
                    writer.writerow([key, value])
    

# Create an argparse.Namespace object from input args.
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
    parser.add_argument("-i", "--input", help="the input lsv file from merge overlap", required=True)
    parser.add_argument("-o", "--output_dir", help="the output directory", required=True)
    parser.add_argument("-g", "--genome_directory", help="path to the directory of genomes", required=True)
    parser.add_argument("-k", "--kmer_length", help="desired k-mer length for comparison", required=False)
    parser.add_argument("-th", "--threshold", help="threshold value of similarity for pruning", required=True)
    return parser.parse_args(argv)


def main():
    ''' Main should take an LSV file as an argument and parse the taxids from it
        to create DNA objects '''
    print("Running Genome Comparison clustering...")
    arguments = parseArgs(argv=sys.argv[1:])

    # find out if single genome or more.
    single_genome = False
    if not os.path.exists(arguments.input + "merge_overlap_out_filtered_genomes.lsv"): single_genome = True

    # instantiate clustering object.
    if single_genome:
        genomeCompareObj = GenomeCompare(arguments.input + "merge_overlap_out_genomes.lsv", 
                                         arguments.genome_directory, 
                                         arguments.kmer_length)
    else:
        genomeCompareObj = GenomeCompare(arguments.input + "merge_overlap_out_filtered_genomes.lsv", 
                                         arguments.genome_directory, 
                                         arguments.kmer_length)
    # find related genomes.
    genomeCompareObj.create_adjacency_matrix()
    genomeCompareObj.display_adjacency_matrix()
    genomeCompareObj.prune_adj_matrix(arguments.threshold)

    # get renewed abundances
    if single_genome:
        genomeCompareObj.estimate_abundances(arguments.input + "merge_overlap_out.csv")
    else: 
        genomeCompareObj.estimate_abundances(arguments.input + "merge_overlap_out_refined.csv")

    # save information to output files
    genomeCompareObj.output_to_file(arguments.output_dir, True)
    genomeCompareObj.output_to_file(arguments.output_dir, False)
    


if __name__ == "__main__":
    main()