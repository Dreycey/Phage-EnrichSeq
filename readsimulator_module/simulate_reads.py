#! usr/bin/python3

# std library
import sys
import subprocess
import os
# non-std library
import argparse
import ntpath
import pandas as pd


ABOUT = (
"""
Description:

This script simulates a fasta file with different percentages of 
nucleotides from a list of input genomes in a config file.

Usage Example:

python simulate_reads.py -i simulate_genomes.config  -c 30 -o simulatedgenomes -t 4
"""
)

# printing the logo; font: 'stop' 
# URL: https://patorjk.com/software/taag/#p=display&f=Stop&t=EnrichSeq%0AGenomeSIm
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


####
# GLOBALS
####
PATH = os.path.dirname(os.path.abspath(__file__))

#####
# SCRIPT BEGINS HERE
#####
def parseArgs(argv=None):
    """
    This method takes in the arguments from the command and performs
    parsing.
    """
    parser = argparse.ArgumentParser(description=ABOUT)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")
    group.add_argument("-q", "--quiet", action="store_true")
    parser.add_argument("-i", "--config", help="input config file", required=True)
    parser.add_argument("-c", "--coverage", help="coverage for output", required=True)
    parser.add_argument("-o", "--output", help="the output file prefix", required=True)
    parser.add_argument("-t", "--threads", help="number of threads", required=True)

    return parser.parse_args(argv)

def parseConfig(configfile):
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

def simulate_reads(genome_array, args):
    """
    This method runs the bash script for simulating genomes.
    """
    directories = []
    for genome in genome_array:
        genome_path = genome[0]
        directory_name = ntpath.basename(genome[0]).strip(".fa").strip(".fasta")
        directory_name = directory_name + "_simulatedreads"
        subprocess.call(f"bash simulate_reads.sh {genome_path} simulatedgenomes {directory_name} {args.threads}", shell=True)
        directories.append(directory_name)

    return directories

def parseFastaToDB(fasta_path):
    """
    Turn simulated reads file into pandas dataframe
    """
    line_index = 0
    fasta_open = open(fasta_path).readlines()
    # parse the fasta file
    name_list, seq_list = [], []
    for line in fasta_open:
        if line[0] == ">":
            name = line.strip(">")
            seq = fasta_open[line_index+1]
            name_list.append(name)
            seq_list.append(seq)

        line_index += 1
    # Make pandas dataframe    
    fa_seqs = {'names': name_list, 'seqs': seq_list}
    dataframe = pd.DataFrame(data=fa_seqs)

    return dataframe

def addDBtoFile(inputDF, outputfile):
    """
    This takes an input pandas datafram with col1 as names and col2 as seqs, 
    then converts this data into an output fasta file.
    """
    for index, row in inputDF.iterrows():
        name, seq = row['names'], row['seqs']
        outputfile.write(">" + name)
        outputfile.write(seq)

def create_simulated_fasta(seqtype, args, simulated_read_path, config_array, output_fasta):
    """
    This method stores the genomes into a specified fasta.
    """
    type_dict = {"illumina" : "/simulatedgenomes_illumina.fa",
                 "nanopore" : "/simulatedgenomes_aligned_reads.fasta",
                 "pacbio" : "/simulatedgenomes_pacbio.fq"}
    
    print(simulated_read_path, config_array)
    #total_nucleotides = args.coverage * 
    total_reads = 100000
    # get output file ready
    filenameout = f"{output_fasta}_{seqtype}.fa"
    if os.path.exists(filenameout):
        os.remove(filenameout) 
    file_out = open(filenameout, "a+")
    # add each genome to simulated file 
    index = 0
    for genome in config_array:
        # get all simulated reads from file, put into DS
        data = simulated_read_path[index] + type_dict[seqtype]
        seqs_dataframe = parseFastaToDB(data)
        # take percent_reads and add to output file
        percent_reads = total_reads * float(genome[1])
        seqs_dataframe = seqs_dataframe.sample(frac = 1)
        addDBtoFile(seqs_dataframe.head(int(percent_reads)), file_out)
        # increment
        index += 1

###
# MAIN
###
if __name__ == "__main__":

    # get the arguments ready 
    argslist = sys.argv[1:] # init argv in case not testing
    args = parseArgs(argslist)
   
    # read config, simulate reads
    config_array = parseConfig(args.config)
    simulated_reads = simulate_reads(config_array, args)

    # create new fasta with simulated genomes
    create_simulated_fasta("illumina", args, simulated_reads, config_array, args.output)
    create_simulated_fasta("nanopore", args, simulated_reads, config_array, args.output)
    create_simulated_fasta("pacbio", args, simulated_reads, config_array, args.output)

#    if args.quiet:
#        print(answer)
#    elif args.verbose:
#        print(f"{args.x} to the power {args.y} equals {answer}")
#    else:
#        print(f"{args.x}^{args.y} == {answer}")


