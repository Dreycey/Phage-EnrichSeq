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
        self.dnaList: List = dnaList 

    
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
        pass


    def prune_adj_matrix(self, threshold = 0.0):
        ''' 
            DESCRIPTION:
                Retrieves DNA objects that are genetically similar 
            INPUT:
                Threshold value between 0.0 and 1.0
            OUTPUT:
                List of Sets. Sets represent clusters of similar DNA objects
        '''
        pass
    

