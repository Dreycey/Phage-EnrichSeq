#! usr/bin/python3
"""
Description:
    This script simulates a fasta file with different percentages of 
    nucleotides from a list of input genomes in a config file.

Usage:
    python simulate_reads.py -i <input config file> -o <output fasta file> -rn <integer number of reads> -t <threads>
Usage Example:
    python simulate_reads.py -i simulate_genomes.config -o simulate_genomes_1000reads -rn 1000 -t 4
"""
# std library
from pickle import LIST
import sys
import subprocess
import os
from pathlib import Path
from dataclasses import dataclass
# non-std library
import argparse
import ntpath
from typing import List, Optional, Tuple, Union
import pandas as pd

# PRINTING LOGO
## URL: https://patorjk.com/software/taag/#p=display&f=Stop&t=EnrichSeq%0AGenomeSIm
LOGO = (
"""
 _______             _       _        _                   
(_______)           (_)     | |      | |                  
 _____   ____   ____ _  ____| | _     \ \   ____ ____     
|  ___) |  _ \ / ___) |/ ___) || \     \ \ / _  ) _  |    
| |_____| | | | |   | ( (___| | | |_____) | (/ / | | |    
|_______)_| |_|_|   |_|\____)_| |_(______/ \____)_|| |    
                                                   |_|    
  ______                                  _   _____       
 / _____)                                | | (_____)      
| /  ___  ____ ____   ___  ____   ____    \ \   _   ____  
| | (___)/ _  )  _ \ / _ \|    \ / _  )    \ \ | | |    \ 
| \____/( (/ /| | | | |_| | | | ( (/ / _____) )| |_| | | |
 \_____/ \____)_| |_|\___/|_|_|_|\____|______(_____)_|_|_|
""")
print(LOGO)

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
    parser.add_argument("-o", "--output", help="the output file prefix", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=True)
    parser.add_argument("-rn", "--readnum", help="number of reads to simulate", required=True)
    return parser.parse_args(argv)

class ReadSimError(Exception):
    """Base class for exceptions occuring in read simulation"""
    pass

class FastaTooSmall(ReadSimError):
    """
    Exception raised when the simulated fasta does not have enough sequences
    to support the simulation.

    Attributes:
        message -- This is the message displayed.
    """

    def __init__(self, message):
        self.message = message

@dataclass
class SimulationConfig:
    """ This datastructure holds the information from a config file,
    and simualates the reads for that config. """

    config_path: Path
    def __post_init__(self):
        self.config_array: List[Tuple] = self.parseConfig()

    def parseConfig(self) -> List[Tuple]:
        """
        This file parses the input config file.
        Input: 
            A config file seperated CSV
            Example:
                # <genome path> <relative abundance>
                ../genomes/Blessica.fasta, 0.15
                ../genomes/D29_genome.fa, 0.50
                ../genomes/perseus_genome.fa, 0.35
        Output:
            An array with the information in a tuple:
                [(../genomes/Blessica.fasta, 0.15),
                ../genomes/D29_genome.fa, 0.50),
                ../genomes/perseus_genome.fa, 0.35)]
        """
        output_array = []
        config_open = open(self.config_path).readlines()
        for line in config_open:
            if (line[0] != "#") and (len(line.replace(" ","")) > 1):
                line_split = line.replace(" ", "").strip("\n").split(",")
                output_array.append(tuple(line_split))
        return output_array

    def simulate_reads(self, threads, num_reads=1000000) -> List[str]:
        """
        This method runs the bash script for simulating genomes.
        Makes a sub command to a lower level bash script. 
        INPUT:
              1. Number of reads to simulated per genome. [Default: 3000000]
              2. Number of threads.
        OUTPUT:
            returns a list of directory names for the processed reads. 
        """
        directories = []
        for genome in self.config_array:
            genome_path = genome[0]
            directory_name = ntpath.basename(genome_path).strip(".fa").strip(".fasta")
            directory_name = Path(directory_name + "_simulatedreads")
            # TODO: Create a data structure object for running sub commands.
            subprocess.call(f"bash simulate_reads.sh {genome_path} simulatedgenomes {directory_name} {threads} {num_reads}", shell=True)
            directories.append(directory_name)
            print(f"working, {directory_name}")
        return directories

class ReadSegmenter:
    """ 
    This class is defined by a type of input read file. 
    The purpose is to split the reads into different proportions.
    """

    def __init__(self, config_array, sim_read_path):
        self.config_array = config_array
        self.simulated_read_path: List[str] = sim_read_path

    def create_simulated_fasta(self, seqtype, read_number, output_fasta_prefix):
        """
        This method stores the genomes into a specified fasta.
        TODO: This should use a generator.
        """
        type_dict = {"illumina" : "/simulatedgenomes_illumina.fa",
                    "nanopore" : "/simulatedgenomes_aligned_reads.fasta",
                    "pacbio" : "/simulatedgenomes_pacbio.fq"}
        # add each genome to simulated file 
        with self.open_fasta_file(output_fasta_prefix, seqtype) as file_out:
            for index, path_and_abundance in enumerate(self.config_array):
                # create a dataframe.
                fasta_path = Path(str(self.simulated_read_path[index]) + type_dict[seqtype])  # TODO: Use path 
                seqs_dataframe: pd.DataFrame = self.parseFastaToDB(fasta_path)
                # add sub sample of data to file.
                sim_abundance = path_and_abundance[1]
                reads_for_genome = float(read_number) * float(sim_abundance)
                self.check_readcount_allowed(reads_for_genome, fasta_path)
                print(f"genome: {self.config_array[index]}, number of reads: {reads_for_genome}")
                seqs_dataframe = seqs_dataframe.sample(n=int(reads_for_genome))
                self.addDBtoFile(seqs_dataframe, file_out, Path(self.config_array[index][0]).name)
                del seqs_dataframe 

    def check_readcount_allowed(self, reads_for_genome, fasta_path): 
        """
        This method ensure that the read count needed in the partition
        is big enough for the simulation. For example, if there were 10 reads
        simulated for a genome and the fraction needed is 100, then there needs to
        be reads simulated.
        
        TODO: so the file is recreated with 100 over the need count. 
        """
        num_reads_in_fasta = self.fasta2readcount(fasta_path)
        if (reads_for_genome > num_reads_in_fasta):
            ERROR_TEXT = f"There were not enough reads in the simulated"
            ERROR_TEXT += f"sequences file, therefor the alloted number of"
            ERROR_TEXT += f"reads is too much. Number of reads needed is {reads_for_genome}"
            ERROR_TEXT += f"and the number of reads in the fasta is: {num_reads_in_fasta}"
            raise FastaTooSmall(ERROR_TEXT)

    def fasta2readcount(self, fasta_path):
        """
        Counts the number of read names in a fasta file
        """
        read_count = 0
        with open(fasta_path) as fasta:
            fasta_lines = fasta.readlines()
            for line in fasta_lines:
                read_count += 1 if (line[0] == ">") else 0
        return read_count

    def open_fasta_file(self, output_fasta_prefix, seqtype):
        """
        Create the output fasta file.
        """
        filenameout = Path(f"{output_fasta_prefix}_{seqtype}.fa")
        if os.path.exists(filenameout):
            os.remove(filenameout) 
        return open(filenameout, "a+")

    def parseFastaToDB(self, fasta_path: Path) -> pd.DataFrame:
        """
        Turn simulated reads file into pandas dataframe
        """
        fasta_path = Path(fasta_path)
        name_list, seq_list = self.fasta2lists(fasta_path)
        # Make pandas dataframe    
        fa_seqs = {'names': name_list, 'seqs': seq_list}
        dataframe = pd.DataFrame(data=fa_seqs)
        return dataframe

    def fasta2lists(self, fasta_path) -> Tuple[List[str]]:
        """
        Takes a fasta path and returns a list of names and sequences.
        """
        name_list, seq_list = [], []
        with open(fasta_path) as fasta_file:
            fasta_files_array = fasta_file.readlines()
            for line_index, line in enumerate(fasta_files_array):
                if line[0] == ">":
                    name = line.strip(">")
                    seq = fasta_files_array[line_index+1]
                    name_list.append(name)
                    seq_list.append(seq)
        return name_list, seq_list

    def addDBtoFile(self, inputDF: pd.DataFrame, outputfile, file_name) -> None:
        """
        This takes an input pandas datafram with col1 as names and col2 as seqs, 
        then converts this data into an output fasta file.
        """
        for index, row in inputDF.iterrows():
            name, seq = row['names'], row['seqs']
            outputfile.write(">" + name.strip("\n") + "|" + file_name + "\n" + seq)

# MAIN
def main():
    # get the arguments ready.
    argslist = sys.argv[1:] 
    args = parseArgs(argslist)
    # simulate reads.
    sim_config = SimulationConfig(args.config)
    simulated_reads_path = sim_config.simulate_reads(threads=args.threads, num_reads=args.readnum)
    # segment reads.
    read_simulator = ReadSegmenter(sim_config.config_array, simulated_reads_path)
    for read_type in ["illumina"]: #, "nanopore", "pacbio"]:
        read_simulator.create_simulated_fasta(seqtype=read_type, 
                                              read_number=args.readnum, 
                                              output_fasta_prefix=args.output)
if __name__ == "__main__":
    main()