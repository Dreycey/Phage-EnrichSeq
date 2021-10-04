import numpy as np
import pandas as pd
from typing import List
from dna import DNA
from pathlib import Path

''' A class that compares between a list of genomes using an adjacency matrix
    of jaccard similarity indices '''


class GenomeCompare:
    def __init__(self, dna_list):
        ''' Initializes GenomeCompare object given a list of DNA objects '''
        self.dna_list: List = dna_list # should set be used instead of list?
    

