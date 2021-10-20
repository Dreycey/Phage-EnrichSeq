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

