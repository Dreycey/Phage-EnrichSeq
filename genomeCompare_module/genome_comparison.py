import numpy as np
import pandas as pd
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


    def prune_adj_matrix(self, threshold = 0.0):
        ''' 
            DESCRIPTION:
                Retrieves DNA objects that are genetically similar based on some threshold
            INPUT:
                Threshold value between 0.0 and 1.0
            OUTPUT:
                A dictionary of all clusters (key: cluster name, value: list of dna objects in cluster). 
        '''
        clusters = {}
        for i in range(len(self.dnaList)-1):
            cluster = set()
            cluster_name = "C_" + str(i+1)
            for j in range(i+1,len(self.dnaList)):
                if self.adjacencyMatrix[i,j] > threshold:
                    cluster.add(self.dnaList[i])
                    cluster.add(self.dnaList[j])
                    clusters[cluster_name] = cluster
        return clusters
    
    
    def display_adjacency_matrix(self):
        # TODO: make pretty
        if self.adjacencyMatrix == None:
            # TODO: raise exception
            print("Cannot compare --no DNA objects provided.")
        else:
            print(pd.DataFrame(self.adjacencyMatrix))

