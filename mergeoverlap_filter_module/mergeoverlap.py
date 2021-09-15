#! usr/bin/python3
"""
DESCRIPTION:
    This script is used for pulling in the names of specific genomes from a line
    delimited files, which correspond to the names in a multi-fasta, then using those
    genomes for more fine grain classification. 

USAGE:
    python mergeoverlap.py <Phage Name line-delimited file> <Multi Fasta file> <out file>

EXAMPLE:

"""
# std packages
from abc import ABC, abstractmethod
import sys
from typing import Dict, List
import csv
import os
import argparse
# non-std packages
import mappy as mp
import matplotlib.pyplot as plt
# in house packages




# Changing to directory of script.
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)

# GLOBALS
PATH = os.path.dirname(os.path.abspath(__file__))

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
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-i", "--config", help="input config file", required=True)
    parser.add_argument("-c", "--coverage", help="coverage for output", required=True)
    parser.add_argument("-o", "--output", help="the output file prefix", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=True)
    parser.add_argument("-rn", "--readnum", help="number of reads to simulate", required=True)
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
            if self.genome_name in first_result:
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
    
    def __init__(self):
        """ initialize all params """
        self.genomes = {} # this directory will hold all of the genomes
        self.genomeMap: Dict[str,MinimapMapper] = {} # this dictionary hold minmap indexes
        self.simAbundance = {} # this dictionary hold simulated abundance amounts
        # add genomes
        self.__addGenomes()
        # keep results saved
        self.resultDict = {} #if empty, can't plot
        # new features
        self.minimap_out = {} # this holds info yeilded from minimap2
        self.genomes_to_view = self.parseGenomeList(genomes_to_view) # this will hold the genome names parse

    def __addGenomes(self):
        """ add genomes to to genomes attr """
        #  parse multifasta  (ASSUMES A GIANT MULTIFASTA FILE)
        print("Creating indexes for minimap2")
        genomeName_list, genome_list = self.parseFasta(self.configFilePath) #TODO: FIX
        for indx in tqdm(range(len(genomeName_list))):
            genome = genome_list[indx]
            genomeName = genomeName_list[indx]
            kracken_out_list = ["Ryadel", "Blessica", "D29", "Paphu", "Perseus"]
            for name in kracken_out_list:
                if  name in genomeName: #in ["Ryadel", "Blessica", "D29", "Paphu", "Perseus"]:
                    print(f"genome name: {genomeName}")
                    # reformat
                    genome_name = genomeName.strip("\n").split(",")[0]
                    genome = genome.strip("\n")
                    # add to the genomes dict
                    self.genomes[genome_name] = genome
                    # create a tempory file for mappy
                    tempfile = open(genome_name+".fa","w")
                    tempfile.write(">"+genome_name+"\n"+genome)
                    # create minimap object
                    self.genomeMap[genome_name] = MinimapMapperWithInfo(name=genome_name, 
                                                                        file_path=genome_name+".fa")
                    # delete the temp file
                    tempfile.close()
                    #os.remove("tempfileDELETE.fa")
    
    def parseGenomeList(self, input_file=None):
        """ parses names seperate by line """
        if (input_file != None):
            genomes_to_view = []
            file_lines = open(input_file).readlines()
            for line in file_lines:
                line = line.strip("\n")
                if len(line) > 0:
                    genomes_to_view.append(line)
            return genomes_to_view
        return self.genomes_to_view

    def print_minimap2output(self):
        """ prints output in minimap_out datastructure """

        for genomeName in self.minimap_out:
            print(genomeName)
            avg_mapq = np.average(self.minimap_out[genomeName]["mapq"])
            print(f"average map quality: {avg_mapq}")
            readcount = self.minimap_out[genomeName]["readcount"]
            print(f"read count: {readcount}")
            #plt.hist(self.minimap_out[genomeName]["mapq"])
            #plt.show()
            overlappMerge = self.overlapMerge(genomeName)
            percenOverlap = self.calcGenomePercetage(genomeName, overlappMerge)
            print(f"percent overlap: {percenOverlap}")

    def __findSeq(self, input_seq):
        """ 
        INPUT 
            1. input read
        OUTPUT
            which genome
            if no genome matching OR more than one, output null 
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

    @staticmethod
    def parseFasta(fasta_path):
        """
        DESCRIPTION - parses a fasta or multifasta file.
        INPUT - 1. fasta path
        OUTPUT - 1. name (ordered list); 2. genome sequence (ordered list)
        """
        seq_names, sequences = [], []
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
        # return dictionary
        self.resultDict = genome_count
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

    def saveResultsAsCSV_1(self, csvout="testing.csv"):
        """
        returns a csv that can be used in regression testing.
        """
        with open(csvout, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.resultDict.keys())
                writer.writeheader()
                writer.writerow(self.resultDict)

    def saveResultAsCSV_2(self, outfilename):
        """
        Save the output as a CSV file.
        """
        if len(self.resultDict.keys()) == 0:
            print("must create a result by running checkSeqFile")
            exit(1)
        # write to a CSV
        with open(outfilename, 'w', newline='') as csvfile:
            fieldnames = ['genome name', 'predicted abundance']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for genome_name, abundance in self.resultDict.items():
                writer.writerow({'genome name': genome_name,
                                 'predicted abundance': abundance})

def main():
    """ Runs the testing script """
    # input
    config_path = sys.argv[1]
    fastaFile = sys.argv[2]
    # Create object
    config_object = GenomeTestSet(config_path) 
    config_object.checkSeqFile(fastaFile) 
    config_object.plotResult("Simulation Validation", out="simTest44.png")
    config_object.saveResultsAsCSV()

if __name__ == "__main__":
    main()
