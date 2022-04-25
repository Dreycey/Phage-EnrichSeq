#! usr/bin/python3
"""
DESCRIPTION:
    This script is used for pulling in the names of specific genomes from a line
    delimited files, which correspond to the names in a multi-fasta, then using those
    genomes for more fine grain classification. 

USAGE:
    python mergeoverlap.py <Phage Name line-delimited file> <Multi Fasta file> <out file>
"""
# std packages
from abc import ABC, abstractmethod
from pathlib import Path
import sys
from typing import Dict, Union, List, Tuple, Optional
import csv
import os
import argparse
import numpy as np
import pickle as pickle
from matplotlib import pyplot
# non-std packages
import matplotlib.pyplot as plt
from numpy import unique
from numpy import where
# in house packages
current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_path}/../PathOrganizer_module")
from PathOrganizer import PathOrganizer, PathErrors, DuplicateGenomeError
from src.modules.mergeoverlap_filter_module.ClusteringModel import TrueGenomeFinder, ClusteringModel, get_true_positive, get_filtered_genomes
from src.modules.mergeoverlap_filter_module.GenomeMapper import GenomeMapper, MinimapMapperWithInfo, MinimapMapper




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
    parser.add_argument("-p", "--plot_results", help="plots extra results", action='store_true', required=False)
    parser.add_argument("-t", "--threads", help="number of threads", required=True)
    return parser.parse_args(argv)

class GenomeTestSet:
    """
    This class is used to test the a folder of simulated sequences
    """
    
    def __init__(self, line_seperated_genomes, genome_directory):
        """ initialize all params """
        # input attributes
        self.input_taxids = [tax_id.strip("\n") for tax_id in open(line_seperated_genomes).readlines()]
        self.object_PathOrganizer = PathOrganizer(genome_directory)
        # primary attributes
        self.genomes = {} # this directory will hold all of the genomes
        self.mappingAdapter = MinimapMapperWithInfo # pointer, instantiated per genome.
        self.genomeMap: Dict[str, GenomeMapper] = {} # this dictionary hold minmap indexes
        self.simAbundance = {} # this dictionary hold simulated abundance amounts
        # add genomes
        self.__addGenomes()
        # keep results saved
        self.resultDict = {} #if empty, can't plot
        self.minimap_out = {} # this holds info yeilded from minimap2

    def __addGenomes(self):
        """ add genomes to to genomes attr """
        print("Creating indexes for minimap2")
        for genome_taxid in self.input_taxids:
            print(f"genome NCBI tax id: {genome_taxid}")
            file_path: Path = self.object_PathOrganizer.genome(genome_taxid)
            if file_path != None: # TODO: IF THIS IS NONE THEN THERE'S A PROBLEM FINDING GENOMES!!
                self.genomes[genome_taxid] = self.parseFasta(file_path)[1][0]
                self.genomeMap[genome_taxid] = self.mappingAdapter(name=genome_taxid,
                                                                     file_path=str(file_path))
            else:
                print(f"FIX THIS: there's a problem finding the genome for {genome_taxid}")

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
                if (genome_from == ""): # first genome with match wins
                    return genome_name
                else:
                    return "UNK"
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
                writer = csv.writer(csvfile)
                #writer.writeheader()
                for genome_name, abundance in self.resultDict.items():
                    writer.writerow([genome_name, abundance])

    def save_features(self):
        """ this method saves the features of each genome"""
        for genomeName in self.minimap_out:
            overlapMerge = self.overlapMerge(genomeName)
            percentOverlap = self.calcGenomePercetage(genomeName, overlapMerge)
            print(f"percent overlap: {percentOverlap}")
            self.minimap_out[genomeName]["overlapMerge"] = overlapMerge
            self.minimap_out[genomeName]["percentOverlap"] = percentOverlap       

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
            self.save_features()

def print_values_for_mappedinfo(dataset: Dict[str,Dict[str, Union[set, int, float]]]) -> Tuple[List[float]]:
    """ 
    takes in the pickle output from the MergeOverlap and prints metrics 
    
    :input (Dict[str,Dict[str, Union[set, int, float]]]):
        The input is a nested dicitonary containing information from the
        primary data structure generated during MergeOverlap. This is
        put into a pickle file and used for data analysis here.
        
    :output (tuple[List[float]):
        1. np.array(y_vector) - a vector containing the names for the phages
        2. np.array(x_vector) - a vector containing a 2D vector of [mapq, overlapmerge percentage]
    """
    y_vector = []
    x_vector = []
    for index, genome in enumerate(dataset.keys()):
        print("\n" + genome)
        y_vector.append(genome)
        x_vector.append([])
        for metric in dataset[genome].keys(): 
            if metric in ['mapq', 'readmaps', 'readcount']:
                print(f" {metric}: {np.average(dataset[genome][metric])}")
                if (metric == 'mapq'):
                    x_vector[index].append(np.average(dataset[genome][metric]))
            elif metric == "overlapMerge":
                print(f" {metric}: {len(dataset[genome][metric]) }")
            elif metric == "percentOverlap":
                print(f" {metric}: {dataset[genome][metric]}")
                x_vector[index].append(dataset[genome][metric])                
            else:
                print(f" {metric}: {dataset[genome][metric]}")
    return np.array(y_vector), np.array(x_vector)

def plot_scatter_plot(x_vector, y_vector, outfile):
    """
    plots a scatter plot of the output values. 
    """
    # plot a scatterplot of the data
    plt.figure(figsize=(8,6))
    plt.scatter(x_vector[:,0], x_vector[:,1], s=300,color="red")
    plt.title("Scatter Plot with annotations",fontsize=15)
    plt.ylabel("Overlap Percentage",fontsize=15)
    plt.xlabel("Average Mapq",fontsize=15)
    for i, label in enumerate(y_vector):
        plt.annotate(label, (x_vector[:,0][i], x_vector[:,1][i]))
    plt.savefig(outfile, dpi=300)
    plt.close()

def scatter_of_filtered(clusters, x_vector, y_vector, outfile):
    """
    plots a scatterplot of the filtered genomes. 
    """
    # retrieve unique clusters
    clusters_set = unique(clusters)
    # create scatter plot for samples from each cluster
    for cluster in clusters_set:
        # get row indexes for samples with this cluster
        row_ix = where(clusters == cluster)
        # create scatter of these samples
        pyplot.scatter(x_vector[row_ix, 0], x_vector[row_ix, 1], s=300)
        for row in enumerate(row_ix):
            x = list(x_vector[row_ix, 0])[0]
            y = list(x_vector[row_ix, 1])[0]
            for i in range(len(x)):
                plt.annotate(y_vector[row_ix[0][i]], (x[i], y[i]))
    plt.savefig(outfile, dpi=300)
    plt.close()




def main():
    print("RUNNING THE MERGE OVERLAP FILTER")
    arguments = parseArgs(argv=sys.argv[1:])

    # Running the algorithm.
    GENOME_DIR = arguments.genome_directory
    genomeTestObj = GenomeTestSet(line_seperated_genomes=arguments.input, genome_directory=GENOME_DIR)
    genomeTestObj.checkSeqFile(arguments.fasta)
    genomeTestObj.saveResultAsCSV(arguments.output_prefix+".csv")
    
    if (arguments.plot_results):
        genomeTestObj.plotResult("Estimated Abundances Using Raw Read Mapping", out=arguments.output_prefix+".png")
        genomeTestObj.print_minimap2output(arguments.output_prefix)
    else:
        genomeTestObj.save_features()

    # save the pickle output.
    predicted_abundances = genomeTestObj.minimap_out
    with open(f'{arguments.output_prefix}_minimap_out.pickle', 'wb') as handle:
        pickle.dump(predicted_abundances, handle)

    # save to original genomes to LSV (note: reduntant with kraken lsv)
    with open(arguments.output_prefix+"_genomes.lsv", "w") as genomes_csv:
        for genome_name in predicted_abundances.keys():
            genomes_csv.write(genome_name + "\n")

    # if more than one genome, use GMM to see if really in sample
    if len(predicted_abundances.keys()) > 2:
        # use unsupervised model to obtain true genomes in sample
        y_vector, x_vector = print_values_for_mappedinfo(genomeTestObj.minimap_out)
        clusters = ClusteringModel().model_predict(x_vector)
        # if (arguments.plot_results):
        #     plot_scatter_plot(x_vector, y_vector, f"{arguments.output_prefix}_scatterplot.png")
        #     scatter_of_filtered(clusters, x_vector, y_vector, f"{arguments.output_prefix}_clusters_scatterplot.png")
        true_pos_genome = get_true_positive(y_vector, x_vector)
        print(f"true positive genome: {true_pos_genome}")
        true_genomes = get_filtered_genomes(true_positive_name=true_pos_genome, 
                                            model=clusters, 
                                            name_list=y_vector)
        # add filtered genomes to output file. 
        with open(arguments.output_prefix+"_filtered_genomes.lsv", "w") as filtered_genomes:
            for genome_name in true_genomes:
                filtered_genomes.write(genome_name + "\n")

        # use filtered genomes (recalcuting abundances affter filtering)
        if os.path.exists(arguments.output_prefix+"_filtered_genomes.lsv"):
            print("THE FILE EXISTS!!!!!!!! contents:")
            with open(arguments.output_prefix+"_filtered_genomes.lsv", "r") as filtered_genomes:
                print(filtered_genomes.readlines())


        genomeTestObj2 = GenomeTestSet(line_seperated_genomes=Path(arguments.output_prefix+"_filtered_genomes.lsv"),
                                    genome_directory=GENOME_DIR)
        genomeTestObj2.checkSeqFile(arguments.fasta)
        genomeTestObj2.saveResultAsCSV(arguments.output_prefix+"_refined.csv")
        if (arguments.plot_results):
            plot_scatter_plot(x_vector, y_vector, f"{arguments.output_prefix}_scatterplot.png")
            scatter_of_filtered(clusters, x_vector, y_vector, f"{arguments.output_prefix}_clusters_scatterplot.png")
            genomeTestObj2.plotResult("Estimated Abundances Using Raw Read Mapping", 
                                      out=arguments.output_prefix+"_refined.png")

if __name__ == "__main__":
    main()
