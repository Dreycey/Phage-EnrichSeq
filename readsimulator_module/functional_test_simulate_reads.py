#! usr/bin/python3
"""
This script is used for functional testing the read simulator.
This tests that the output from simulation matches the expected amounts
of DNA per genome. Black Box Functional testing

USAGE:
python scriptname.py <config_file> <simulated_read_file>

EXAMPLE:
python test_simulate_reads.py simulate_genomes.config simulatedgenomes_illumina.fa
"""
# std packages
from abc import ABC, abstractmethod
import sys
from typing import Dict, List
import csv
# non-std packages
import mappy as mp
import matplotlib.pyplot as plt
# in house packages
from simulate_reads import SimulationConfig




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

@dataclass
class GenomeName:
    """ this data class hold the names corresponding to genomes. """
    EMPTY: str = "EMPTY"
    UNKOWN: str = "UNK"
    MULTIPLE: str = "MULTIPLE"

class GenomeTestSet:
    """
    This class is used to test the a folder of simulated sequences
    """
    
    def __init__(self, configFile):
        """ initialize all params """
        self.configFilePath = configFile
        self.config: SimulationConfig = SimulationConfig(configFile)
        self.genomes = {} # this directory will hold all of the genomes
        self.genomeMap: Dict[str,MinimapMapper] = {} # this dictionary hold minmap indexes
        self.simAbundance = {} # this dictionary hold simulated abundance amounts
        # add genomes
        self.config_info = self.__parseConfig()
        self.__addGenomes()
        # keep results saved
        self.resultDict = {} #if empty, can't plot

    def __parseConfig(self):
        """ parse the contents of the config file """       
        return self.config.config_array

    def __addGenomes(self):
        """ add genomes to to genomes attr """
        for genome_pathAndAbundance in self.config_info:
            file_path, amount = genome_pathAndAbundance
            genome_name, genome_list = self.parseFasta(file_path)
            genome_name = genome_name[0].strip("\n").split(",")[0]
            genome = genome_list[0].strip("\n")
            # TODO: making attributes below, is this the best way??
            self.genomes[genome_name] = genome # add to the genomes dict
            self.genomeMap[genome_name] = MinimapMapper(genome_name, file_path) # create minimap object/index
            self.simAbundance[genome_name] = amount # simulated amounts
    
    def __findSeq(self, input_seq):
        """ 
        DESCRIPTION:
            This method finds a corresponding genome for an input read.
        INPUT 
            1. input read
        OUTPUT
            which genome
            if no genome matching OR more than one, output null 
        """
        genome_that_read_maps_to: str = GenomeName.EMPTY
        # find genome that read is in
        for genome_name in self.genomes.keys():
            minimap_mapper: MinimapMapper = self.genomeMap[genome_name]
            read_maps_to_genome: bool = minimap_mapper.does_read_map(input_seq)
            if read_maps_to_genome:
                if (genome_that_read_maps_to == GenomeName.EMPTY):
                    genome_that_read_maps_to = genome_name
                else:
                    return GenomeName.MULTIPLE
        # if nothing, return UNK
        if genome_that_read_maps_to == GenomeName.EMPTY:
            return GenomeName.UNKOWN
        return genome_that_read_maps_to

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
        genome_count = self.count_reads_mapping_per_genome(seqs)
        # normalize the results
        for genome in genome_count.keys():
            genome_count[genome] = genome_count[genome] / total_reads
        # return dictionary
        self.resultDict = genome_count
        return genome_count

    def count_reads_mapping_per_genome(self, seqs):
        """
        Method returns a dictionary counter for the 
        number of reads mapping per genome. 
        """
        genome_count = {}
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
        # create figs
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.suptitle(plotTitle)
        fig.set_size_inches(20,10)
        # create top plot
        values = [float(i) for i in self.simAbundance.values()]
        ax1.bar(self.simAbundance.keys(), 
                values, 
                color=PLOTCOLOR,
                edgecolor='black')
        ax1.set_ylabel('True Simulated Abundaces')
        # create bottom plot
        ax2.bar(self.resultDict.keys(), 
                list(self.resultDict.values()), 
                color=PLOTCOLOR,
                edgecolor='black')
        ax2.set_ylabel('Estimated Abundances using Minimap2')
        # save the plot
        plt.savefig(out, dpi=300)

    def saveResultsAsCSV(self, csvout="testing.csv"):
        """
        returns a csv that can be used in regression testing.
        """
        with open(csvout, 'w') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.resultDict.keys())
                writer.writeheader()
                writer.writerow(self.resultDict)

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
