#! usr/bin/python3
"""
DESCRIPTION:
    This script is used for pulling in the names of specific genomes from a line
    delimited files, which correspond to the names in a multi-fasta, then using those
    genomes for more fine grain classification. 

USAGE:
    python mergeoverlap.py <Phage Name line-delimited file> <Multi Fasta file> <out file>

EXAMPLE:

TODO: 
    1. Make sure the genomic mapping can be parrallized. 
        1.A - see if threads can be used for individual mapping.
    2. Look for dipps in coverage along the genome.
        2.A - perhaps plot this as well! 
"""
# std packages
from abc import ABC, abstractmethod
import sys
from typing import Dict, List
import csv
import os
import argparse
import numpy as np
# non-std packages
import mappy as mp
import matplotlib.pyplot as plt
# in house packages
sys.path.append("/Users/dreyceyalbin/Desktop/Phage-EnrichSeq/PathOrganizer_module")
from PathOrganizer import PathOrganizer, PathErrors, DuplicateGenomeError




# Changing to directory of script.
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)

# GLOBALS.
PATH = os.path.dirname(os.path.abspath(__file__))

# Create an argparse.Namespace object from input args.
def parseArgs(argv=None) -> argparse.Namespace:
    """
    This method takes in the arguments from the command and performs
    parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    """
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-i", "--input", help="the input lsv file with phage names", required=True)
    parser.add_argument("-o", "--output_prefix", help="the output file prefix", required=True)
    parser.add_argument("-f", "--fasta", help="input fasta path", required=True)
    parser.add_argument("-g", "--genome_directory", help="path to the directory of genomes", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=True)
    return parser.parse_args(argv)

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

class GenomeTestSet:
    """
    This class is used to test the a folder of simulated sequences
    """
    
    def __init__(self, line_seperated_genomes,genome_directory):
        """ initialize all params """
        # input attributes
        self.input_genomes = [genome_name.strip("\n").split(" ")[2] 
                              for genome_name in open(line_seperated_genomes).readlines()]
        self.object_PathOrganizer = PathOrganizer(genome_directory)
        # primary attributes
        self.genomes = {} # this directory will hold all of the genomes
        self.genomeMap: Dict[str, MinimapMapper] = {} # this dictionary hold minmap indexes
        self.simAbundance = {} # this dictionary hold simulated abundance amounts
        # add genomes
        self.__addGenomes()
        # keep results saved
        self.resultDict = {} #if empty, can't plot
        self.minimap_out = {} # this holds info yeilded from minimap2

    def __addGenomes(self):
        """ add genomes to to genomes attr """
        print("Creating indexes for minimap2")
        for genome_to_grab in self.input_genomes:
            print(f"genome name: {genome_to_grab}")
            file_path = self.object_PathOrganizer.genome(genome_to_grab)
            if file_path != None: # TODO: IF THIS IS NONE THEN THERE'S A PROBLEM FINDING GENOMES!!
                self.genomes[genome_to_grab] = self.parseFasta(file_path)[1][0]
                self.genomeMap[genome_to_grab] = MinimapMapperWithInfo(name=genome_to_grab,
                                                            file_path=str(file_path))
            else:
                print(f"FIX THIS: there's a problem finding the genome for {genome_to_grab}")

    @staticmethod
    def parseFasta(fasta_path):
        """
        DESCRIPTION - parses a fasta or multifasta file.
        INPUT - 1. fasta path
        OUTPUT - 1. name (ordered list); 2. genome sequence (ordered list)
        """
        seq_names, sequences = [], []
        print(f"fasta file: {fasta_path}")
        with open(fasta_path) as fasta_file:
            fasta_input = fasta_file.readlines()
            sequence_i = ""
            for counter, line in enumerate(fasta_input):
                if (line[0] == ">"):
                    seq_name = line[1:]
                    seq_names.append(seq_name.strip("\n"))
                    if (counter != 0):
                        sequences.append(sequence_i.strip("\n"))
                        sequence_i = "" # renew sequence
                else:
                    sequence_i += line.strip("\n")
            sequences.append(sequence_i.strip("\n"))
        return seq_names, sequences     

    def checkSeqFile(self, input_fasta):
        """
        iterate through fasta and generate results
        INPUT
            fasta file path
        OUTPUT
            {Blessia : 0.3, D29: 0.4, ... Persues: 0.3 }
        """
        # init
        seq_names, seqs = self.parseFasta(input_fasta)
        total_reads = len(seqs)
        genome_count: Dict[str,int] = self.count_reads_mapping_per_genome(seqs)
        # normalize the results
        for genome in genome_count.keys():
            genome_count[genome] = genome_count[genome] / total_reads
        self.resultDict = genome_count
        # return dictionary
        for minimap_obj in self.genomeMap.values():
            for gen, minimap_output in minimap_obj.minimap_out.items():
                self.minimap_out[gen] = minimap_output
        return genome_count

    def count_reads_mapping_per_genome(self, seqs) -> Dict[str,int]:
        """
        Method returns a dictionary counter for the 
        number of reads mapping per genome. 
        """
        genome_count: Dict[str,int] = {}
        for sequence in seqs:
            genome = self.__findSeq(sequence)
            if genome in genome_count.keys():
                genome_count[genome] += 1
            else:
                genome_count[genome] = 1
        return genome_count

    def __findSeq(self, input_seq):
        """ 
        INPUT 
            1. input read
        OUTPUT
            which genome
            if 1. no genome matching, output 'UNK'.
               2. more than one, output 'MULTIPLE'.
               3. exactly one, output the name of the corresponding genome.
        """
        genome_from = ""
        # find genome that read is in
        for genome_name in self.genomes.keys():
            minimap_mapper: MinimapMapper = self.genomeMap[genome_name]
            read_maps_to_genome: bool = minimap_mapper.does_read_map(input_seq)
            if read_maps_to_genome:
                if (genome_from == ""):
                    return genome_name
                else:
                    return "MULTIPLE"
        if genome_from == "":  # if nothing, return UNK
            return "UNK"
        return genome_from

    def overlapMerge(self, genome_name):
        """ extend reads to assess % genome covered

        Algorithm:
            Example:
                1. input
                (50,150), (0,100), (200,300), (149,249)
                2. sort by starting reference index
                (0,100), (50,150), (149,249), (200,300)
                3. merge
                    iteration one:
                        p1 = (0,100)
                        p2 = (50, 150)
                        if p1[1] >= p2[0]:
                            merge
                            p1 = merged
                            p2 = move forward
                        else:
                            p1 = p2
                            p2 = move forward

        INPUT:
            self.minimap_out[genomeName]["readmaps"] (set object)
        OUTPUT:
            return
                1. merged set
                2. % of genome with reads mapped
        """
        mergedSet = set()
        #for genome_name in self.genomes.keys():
        # run some type sorting algorithm

        #1. input
        mapped_reads = self.minimap_out[genome_name]["readmaps"] #[(50,120),(110,220),(150,200)]
        print(f"before merge overlap: {len(mapped_reads)}")
        #2. sort the reads by starting index
        sorted_mapped_reads = sorted(mapped_reads, key=lambda x: x[0])

        #3. merge overlapping regions
        p1 = 0 # pointer 1
        p2 = 1 # pointer 2
        p1_start, p1_end = sorted_mapped_reads[p1]
        while p2 < len(sorted_mapped_reads):
            p2_start, p2_end = sorted_mapped_reads[p2]
            if p1_end >= p2_start:  # if overlap,  then merge
                if p1_end <= p2_end:
                    p1_end = p2_end
            else:                   # if no overlap, set p1=p2 and p2++
                mergedSet.add((p1_start, p1_end))
                p1_start, p1_end = sorted_mapped_reads[p2] #notice: using p2

            p2 += 1
        # add last interval to to the merged set
        mergedSet.add((p1_start, p1_end))

        print(f"after merge overlap: {len(mergedSet)}")
        return mergedSet

    def calcGenomePercetage(self, genome_name, merged_set):
        """
        This method calculate the percentage of the genome
        covered by the overlapping intervals after merging.

        input
            overlap merge set
        output
            percentage of overlap for the mergedset
        """
        # initialize
        overall_map_length = 0
        genome_length = len(self.genomes[genome_name])
        # calculate # of bases with mapped reads
        for interval in merged_set:
            interval_size = abs(interval[0] - interval[1])
            overall_map_length += interval_size
        # calculate % of genome with mapped reads
        print(f"overlapp length: {overall_map_length}")
        overlap_percent = overall_map_length / genome_length

        return overlap_percent

    def plotResult(self, plotTitle, out=None, result=None):
        """
        Plot bar plot of the results from the simulation
        """
        if len(self.resultDict.keys()) == 0:
            print("must create a result by running checkSeqFile")
            exit(1)
        # plot color
        PLOTCOLOR='orange'
        TITLESIZE=40
        AXISSIZE=25
        # create figs
        fig, (ax2) = plt.subplots(1, 1)
        fig.suptitle(plotTitle, size=TITLESIZE)
        fig.set_size_inches(20,10)
        # create bottom plot
        ax2.bar(self.resultDict.keys(),
                list(self.resultDict.values()),
                color=PLOTCOLOR,
                edgecolor='black')
        ax2.set_ylabel('Estimated Abundances using Minimap2',size=AXISSIZE)
        ax2.set_xticklabels(self.resultDict.keys(), fontsize=15)
        # save the plot
        plt.savefig(out, dpi=300)

    def saveResultAsCSV(self, csvout="testing.csv"):
        """
        returns a csv that can be used in regression testing.
        """
        with open(csvout, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.resultDict.keys())
                writer.writeheader()
                writer.writerow(self.resultDict)

    def print_minimap2output(self, out_prefix):
        """ prints output in minimap_out datastructure """
        for genomeName in self.minimap_out:
            print(genomeName)
            avg_mapq = np.average(self.minimap_out[genomeName]["mapq"])
            print(f"average map quality: {avg_mapq}")
            readcount = self.minimap_out[genomeName]["readcount"]
            print(f"read count: {readcount}")
            plt.hist(self.minimap_out[genomeName]["mapq"])
            plt.title(f"mapq Histogram for {genomeName}")
            plt.savefig(out_prefix+"_"+genomeName+".png", dpi=300)
            plt.clf()
            overlappMerge = self.overlapMerge(genomeName)
            percenOverlap = self.calcGenomePercetage(genomeName, overlappMerge)
            print(f"percent overlap: {percenOverlap}")

def main():
    print("RUNNING THE MERGE OVERLAP FILTER")
    arguments = parseArgs(argv=sys.argv[1:])

    #TESTING.
    output_testing = open(arguments.output_prefix+"_testingArguments", "w")

    # Running the algorithm.
    GENOME_DIR = arguments.genome_directory
    genomeTestObj = GenomeTestSet(line_seperated_genomes=arguments.input, genome_directory=GENOME_DIR)
    genomeTestObj.checkSeqFile(arguments.fasta)
    genomeTestObj.plotResult("Estimated Abundances Using Raw Read Mapping", out=arguments.output_prefix+".png")
    genomeTestObj.saveResultAsCSV(arguments.output_prefix+".csv")
    genomeTestObj.print_minimap2output(arguments.output_prefix)

    # save to file
    output_testing.write("RESULT: \n"+str(genomeTestObj.resultDict))

if __name__ == "__main__":
    main()
