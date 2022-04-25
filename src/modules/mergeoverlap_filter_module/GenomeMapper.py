"""
Module for mapping reads to a genome. This module acts as a wrapper to an underlying
read mapper, allowing for flexible swapping. 

Classes
    1. GenomeMapper - abstract class for defining each mapping object. 
    2. MinimapMapperWithInfo - concrete class using minimap as a mapper. This class
                               also saved inforrmation about the mapped reads.
    3. MinimapMapper - concrete class using minimap as a mapper.

Methods
    N/A
"""
from abc import ABC, abstractmethod
from typing import Dict, Union, List, Tuple, Optional
import mappy as mp




class GenomeMapper(ABC):
    """ This acts as an adapter for mapping reads to a genome """

    def __init__(self, name):
        self._genome_name = name
 
    @abstractmethod
    def map_fasta_read(self, sequence):
        """
        this method is used to map an individual read.
        """
        pass

    @property
    def genome_name(self):
        """
        getter for the genome name. 
        """
        return self._genome_name

    @genome_name.setter
    def genome_name(self, new_name):
        """ set the genome name. """
        self._genome_name = new_name

    @abstractmethod
    def does_read_map(self, sequence) -> bool:
        """
        checks if the input read maps to the given genome.
        """
        pass 

class MinimapMapperWithInfo(GenomeMapper):
    """ mapping object adapter for minimap """

    def __init__(self, name, file_path):
        self._genome_name = name
        self.index_object = mp.Aligner(file_path, best_n=1)
        self.minimap_out = {}

    def map_fasta_read(self, sequence):
        """
        this method is used to map an individual read.
        """
        #print("{}\t{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.mapq, hit.mlen))
        return self.index_object.map(sequence)

    def does_read_map(self, sequence) -> bool:
        """
        checks if the input read maps to the given genome.
        if it does, get information about the mapped read. 
        """
        try:
            first_result = next(self.map_fasta_read(sequence))
            if self.genome_name in self.minimap_out:
                self.minimap_out[self.genome_name]["mapq"].append(first_result.mapq)
                #self.minimap_out[genomeName]["readmaps"].add((first_result.r_st, first_result.r_en))
                self.minimap_out[self.genome_name]["readmaps"].append((first_result.r_st, first_result.r_en))
                self.minimap_out[self.genome_name]["readcount"] += 1
            else:
                self.minimap_out[self.genome_name] = {} # create if first time
                self.minimap_out[self.genome_name]["mapq"] = [first_result.mapq]
                #self.minimap_out[genomeName]["readmaps"] = set((first_result.r_st, first_result.r_en))
                self.minimap_out[self.genome_name]["readmaps"] = [(first_result.r_st, first_result.r_en)]
                self.minimap_out[self.genome_name]["readcount"] = 1
        except StopIteration:
            return False
        return True

class MinimapMapper(GenomeMapper):
    """ mapping object adapter for minimap """

    def __init__(self, name, file_path):
        self._genome_name = name
        self.index_object = mp.Aligner(file_path, best_n=1)

    def map_fasta_read(self, sequence):
        """
        this method is used to map an individual read.
        """
        #print("{}\t{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.mapq, hit.mlen))
        return self.index_object.map(sequence)

    def does_read_map(self, sequence) -> bool:
        """
        checks if the input read maps to the given genome.
        """
        try:
            first_result = next(self.map_fasta_read(sequence))
        except StopIteration:
            return False
        return True