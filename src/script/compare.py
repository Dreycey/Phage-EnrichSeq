#! usr/bin/python3
"""
co=ompare karaken to merge overlap


Usage:
python compare.py <simulated_fasta> <CSV>
"""
import sys
import json 

def make_dir():
    ncbi2taxid = "database_prev/krakenDB/seqid2taxid.map"
    dictionary = {}
    with open(ncbi2taxid) as ncbi2taxid:
        lines = ncbi2taxid.readlines()
        for line in lines:
            ncbi, taxid = line.replace("\n","").split("\t")
            dictionary[ncbi] = taxid
    return dictionary

def parse_simulated_fasta(fasta_path):
    """                          
    Description
    -----------
        This function parses the simulated true fasta file used
        for testing.
        
    Input
    -----------
        1. Path to a singular simulated fasta file
        
    Output
    -----------
        Abundance Dictionary containing:
            1. taxid abundances
            abundances_dict = {'tax_id_1' : 33,
                               'tax_id_2' : 33,
                               ...,
                               'tax_id_N' : 33,
                              }
    """
    abundances = {'taxid_abundance' : {}}
    total_counts = 0
    dictionary = make_dir()
    with open(fasta_path, "r") as fasta_opened:
        line = fasta_opened.readline()
        while (line):
            if (line[0] == ">"):
                ncbi_id = line[1:].split("|")[0].split("-")[0]
                # dictioinary from FastViromeExplorer, obtained above
                if ncbi_id in dictionary:
                    name = dictionary[line[1:].split("|")[0].split("-")[0]]
                else:
                    print(f"DICT DOES NOT HAVE NAME {ncbi_id}")
                    #exit(1)
                    name = "unk"
                # increment abundance level per read count
                if name in abundances['taxid_abundance']:
                    abundances['taxid_abundance'][name] += 1
                else:
                    abundances['taxid_abundance'][name] = 1
                total_counts += 1
                # print update
                if (total_counts % 1000000 == 0):
                    print(f"{total_counts} lines parsed")
            line = fasta_opened.readline()
            
    #turn counts into abundances - normalize
    for tax_id in abundances['taxid_abundance'].keys():
        abundances['taxid_abundance'][tax_id] /= total_counts #normalize
    return abundances

def parse_file(file_in, delim=""):
    """ parses an input filie """
    name_list = {}
    lines = open(file_in).readlines()
    for line in lines:
        if delim == "":
            split_on_delim = line.strip("\n")
            name_list[split_on_delim] = 1
        else:
            split_on_delim = line.split(delim)[0]
            name_list[split_on_delim] = 1
    return name_list

def compare_dictionaries(l1, l2):
    """ count differences """
    count = 0
    print(l1, l2)
    for name in l1.keys():
        if name not in l2:
            count += 1
            print(name)
    print(f"this many false negatives {count}")
    print(f"truth is {len(l1.keys())} and we ound {len(l2.keys())}")
    count = 0  
    for name in l2.keys():
        if name not in l1:
            count += 1
    print(f"this many false positives {count}")

if __name__ == "__main__":

    if (len(sys.argv) == 3):
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    else:
        print(__doc__)
        exit(1)

    d1=parse_simulated_fasta(file1)['taxid_abundance']
    d2=parse_file(file2, delim=",")
    print("BEFORE")
    compare_dictionaries(d1, d2)
