#! usr/bin/python3

# std packages
import sys
# non-std packages
#from Bio import pairwise2
import mappy as mp
import matplotlib.pyplot as plt
import csv
# in house packages
from simulate_reads import parseConfig as parse_config

"""
This script tests that the output from simulation matches the expected amounts
of DNA per genome. Black Box testing

USAGE:
python scriptname.py config_file simulated_read_file

EXAMPLE:
python simulate_reads_test.py simulate_genomes.config simulatedgenomes_illumina.fa 
"""


class GenomeTestSet:
    """
    This class is used to test the a folder of simulated sequences
    """
    
    def __init__(self, configFile):
        """ initialize all params """
        self.configFilePath = configFile
        self.genomes = {} # this directory will hold all of the genomes
        self.genomeMap = {} # this dictionary hold minmap indexes
        self.simAbundance = {} # this dictionary hold simulated abundance amounts

        # add genomes
        self.config_info = self.__parseConfig()
        self.__addGenomes()

        # keep results saved
        self.resultDict = {} #if empty, can't plot

    def __parseConfig(self):
        """ parse the contents of the config file """
        return parse_config(self.configFilePath) 

    def __addGenomes(self):
        """ add genomes to to genomes attr """
        for genome_path in self.config_info:
            amount = genome_path[1]
            file_path = genome_path[0]
            #print(f"parsing the genome: {file_path}")
            genome_name, genome_list = self.parseFasta(file_path)
            # reformat
            genome_name = genome_name[0].strip("\n").split(",")[0]
            genome = genome_list[0].strip("\n")
            # add to the genomes dict
            self.genomes[genome_name] = genome
            # create minimap object
            self.genomeMap[genome_name] = mp.Aligner(file_path, best_n=1)
            # simulated amounts
            self.simAbundance[genome_name] = amount

    def mapper_1(self, genomeName, read):
        """
        This mapper uses a wrapper to minimap2

        OUTPUT
            True, if read in genome
            False, else
        """
        for hit in self.genomeMap[genomeName].map(read): # traverse alignments
            #print(f"length of read: {len(read)}, {counter}")
            #print("{}\t{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.mapq, hit.mlen))
            return True
        return False
    
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
        for genome_name, genome_seq in self.genomes.items():
           # print(f"genome: {genome_name}")
            if (self.mapper_1(genome_name, input_seq)):
                if (genome_from == ""):
                    genome_from = genome_name
                else:
                    continue #return "UNK"
            else:
                continue
        # if nothing, return UNK
        if genome_from == "":
            return "UNK"
        else:
            return genome_from

    def parseFasta(self, fasta_path):
        """
        DESCRIPTION
            parses a fasta or multifasta file
        
        INPUT
            fasta
        OUTPUT
            name (list)
            genome sequence (list)
        """
        # parse the input fasta
        fasta_input = open(fasta_path).readlines()
        seq_names = []
        seqs = []
        counter = 0
        seq = ""
        for line in fasta_input:
            if line[0] == ">":
                seq_name = line[1:]
                seq_names.append(seq_name)
                if counter != 0:
                    seqs.append(seq.strip("\n"))
                    seq = "" # renew sequence
            else:
                seq += line.strip("\n")
            counter += 1
        seqs.append(seq.strip("\n"))

        return seq_names, seqs        

    def checkSeqFile(self, input_fasta):
        """
        iterate through fasta and generate results

        INPUT
            fasta file path

        OUTPUT
            {   
                Blessia : 0.3
                D29: 0.4
                Persues: 0.3
                ...
                UNK: 0.01
            }
        """
        # init
        seq_names, seqs = self.parseFasta(input_fasta)
        genome_count = {}
        total_reads = 0
        # loop through fasta seqs
        for sequence in seqs:
            genome = self.__findSeq(sequence)
            if genome in genome_count.keys():
                genome_count[genome] += 1
            else:
                genome_count[genome] = 1
            total_reads += 1

        # normalize the results
        for genome, read_count in genome_count.items():
            genome_count[genome] = genome_count[genome] / total_reads
        # return dictionary
        self.resultDict = genome_count
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
        fig.suptitle(plotTitle, size = 30)
        fig.set_size_inches(20,10)
        # create top plot
        values = [ float(i) for i in self.simAbundance.values()]
        ax1.bar(self.simAbundance.keys(), 
                values, 
                color=PLOTCOLOR,
                edgecolor='black')
        ax1.set_ylabel('True Simulated Abundaces', size = 20)
        # create bottom plot
        ax2.bar(self.resultDict.keys(), 
                list(self.resultDict.values()), 
                color=PLOTCOLOR,
                edgecolor='black')
        ax2.set_ylabel('Estimated Abundances using Minimap2')
        # save the plot
        plt.savefig(out, dpi=300)

    def saveResultAsCSV(self, outfilename):
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
    config_object.plotResult("Pipeline Result", out="pipelineBench.png") 
    config_object.saveResultAsCSV("pipelineBench")

if __name__ == "__main__":
    main()
