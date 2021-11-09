import sys
import numpy as np
import pandas as pd
import argparse
from typing import List
from dna import DNA
from pathlib import Path

''' A class that compares between a list of genomes using an adjacency matrix
    of jaccard similarity indices '''


class GenomeCompare:
    def __init__(self, dnaList):
        ''' Initializes GenomeCompare object given a list of DNA objects '''
        self.dnaList: List = dnaList # TODO: implement mechanism to prevent duplicates
        self.adjacencyMatrix = np.zeros((len(dnaList), len(dnaList)), float)

    
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
        # update matrix values
        for i in range(len(self.dnaList)):
            for j in range(len(self.dnaList)):
                self.adjacencyMatrix[i,j] = self.dnaList[i].calc_jaccard(self.dnaList[j])
        return self.adjacencyMatrix

    def prune_adj_matrix(self, threshold = 0.0):
        ''' 
            DESCRIPTION:
                Retrieves DNA objects that are genetically similar based on some threshold
            INPUT:
                Threshold value between 0.0 and 1.0 (default=0.0)
            OUTPUT:
                A dictionary of all clusters (key: cluster name, value: list of dna objects in cluster). 
                Clusters are represented by the set data structure. Can have a cluster of size=1.
        '''
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
        
        return clusters
    
    
    def display_adjacency_matrix(self):
        # TODO: make pretty
        if self.adjacencyMatrix == None:
            # TODO: raise exception
            print("Cannot compare --no DNA objects provided.")
        else:
            print(pd.DataFrame(self.adjacencyMatrix))


    def output_to_file(self, file_prefix):
        '''
        Prints cluster information to a csv file

        INPUT:
            Prefix for file to create
        
        OUTPUT:
            CSV file name
            CSV file example:
                C_1, taxid1, taxid3
                C_2, taxid2
                C_3, taxid4, taxid5
        '''
        pass

# Create an argparse.Namespace object from input args.
def parseArgs(argv=None) -> argparse.Namespace:
    '''
    This method takes in the arguments from the command and performs
    parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    # group.add_argument("-v", "--verbose", action="store_true")
    # parser.add_argument("-i", "--input", help="the input lsv file with phage names", required=True)
    parser.add_argument("-o", "--output_prefix", help="the output file prefix", required=True)
    # parser.add_argument("-f", "--fasta", help="input fasta path", required=True)
    parser.add_argument("-g", "--genome_directory", help="path to the directory of genomes", required=True)
    # parser.add_argument("-t", "--threads", help="number of threads", required=True)
    return parser.parse_args(argv)


def main():
    ''' Main should take an LSV file as an argument and parse the taxids from it
        to create DNA objects '''
    print("Running Genome Comparison clustering...")
    arguments = parseArgs(argv=sys.argv[1:])

    GENOME_DIR = arguments.genome_directory


if __name__ == "__main__":
    main()