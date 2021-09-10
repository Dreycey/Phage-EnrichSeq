#! usr/bin/python3

# std packages
import sys
import csv
import os
# non-std packages
import mappy as mp
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

# in house packages

"""
This wrapper impliments using minimap2 to measure abundance compared to
a config file used for simulations.

USAGE:
python minimap2wrapper.py config_file simulated_read_file outputprefix

EXAMPLE:
python minimap2_module/minimap2wrapper.py simulate_genomes.config simulatedgenomes_illumina.fa run_1

OUTPUT 1: (CSV)
genome_name,abundance

OUTPUT 2: (PNG)
This script also outputs a png file.
"""

###
# Stand alone methods
###

def parse_config(configfile):
    """
    This file parses the input config file
    """
    output_array = []
    config_open = open(configfile).readlines()
    for line in config_open:
        if line[0] != "#":
            line_split = line.replace(" ", "").strip("\n").split(",")
            output_array.append(tuple(line_split))

    return output_array

###
# Classes
###
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
        self.minimap_out = {} # this holds info yeilded from minimap2

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
        """
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


        """
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
                    tempfile.write(">"+genome_name)
                    tempfile.write("\n")
                    tempfile.write(genome)
                    # create minimap object
                    self.genomeMap[genome_name] = mp.Aligner(genome_name+".fa", best_n=1)
                    # delete the temp file
                    tempfile.close()
                    #os.remove("tempfileDELETE.fa")

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
            #print(f" Start mapped read: {hit.r_st}, end: {hit.r_en}, mapq: {hit.mapq}")
            #print(f"{genomeName}")

            if genomeName in self.minimap_out.keys():
                self.minimap_out[genomeName]["mapq"].append(hit.mapq)
                #self.minimap_out[genomeName]["readmaps"].add((hit.r_st, hit.r_en))
                self.minimap_out[genomeName]["readmaps"].append((hit.r_st, hit.r_en))
                self.minimap_out[genomeName]["readcount"] += 1
            else:
                self.minimap_out[genomeName] = {} # create if first time
                self.minimap_out[genomeName]["mapq"] = [hit.mapq]
                #self.minimap_out[genomeName]["readmaps"] = set((hit.r_st, hit.r_en))
                self.minimap_out[genomeName]["readmaps"] = [(hit.r_st, hit.r_en)]
                self.minimap_out[genomeName]["readcount"] = 1

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
        genome_from = []
        # find genome that read is in
        for genome_name, genome_seq in self.genomes.items():
            if (self.mapper_1(genome_name, input_seq)):
                genome_from.append(genome_name)
            else:
                continue

        # if nothing, return UNK
        if (len(genome_from) == 0):
            return ["UNK"]
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
        print("Mapping reads")
        for sequence in tqdm(seqs):
            genome_list = self.__findSeq(sequence)
            for genome in genome_list:
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


"""
python minimap2wrapper.py phageMulti.fa simulatedgenomes_illumina.fasta run_1
"""
def main():
    """ Runs the minimap2 script """
    # input
    config_path = sys.argv[1]
    fastaFile = sys.argv[2]
    outfileprefix = sys.argv[3]
    # Create object
    config_object = GenomeTestSet(config_path)
    config_object.checkSeqFile(fastaFile)
    config_object.plotResult("Simulation Validation", out=outfileprefix+".png")
    config_object.saveResultAsCSV(outfileprefix+".csv")
    config_object.print_minimap2output()

if __name__ == "__main__":
    main()
