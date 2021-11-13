"""
This script adds a corresponding taxid to the header for an NCBI genome.

YOU AREN'T USING THE SCRIPT CORRECTLY, CORRECT USAGE IS BELOW:
(if occuring while performing automated build, email us! Problems were not observed during testing.)


USAGE:
    python ncbi2krakenHeader.py path/to/nucl_gb.accession2taxid path/to/genome_directory/
NOTE: 
    puts everything into memory!
INPUT:
    1. NCBI genbank file path
OUTPUT:
    1. The file is modified
"""
import sys
import os
from tqdm import tqdm
from pathlib import Path




def get_acc2tax_map(genome_directory):
    """
    modifes each genome in a directory, assumes is NCBI
    and adds taxid.
    """
    print("    Obtaining NCBI accessions")
    accession_dict = {}
    for filename in os.listdir(genome_directory):
        # READ HEADER AND FIND TAX ID
        full_path = Path(genome_directory) / Path(filename)
        with open(full_path, "r") as fasta_opened:
            file_line = fasta_opened.readline()
            accession_id = file_line.split(" ")[0].strip(">")
            accession_dict[accession_id] = 0
    return accession_dict

def accession2taxid(genbank_file_path, accession_dict, delimiter="\t"):
    """
    This function takes in a genabk file and
    returns a dictoinary mappinig the accession id
    to a NCBI tax id. 
    """
    accession2taxid = {}
    print("    Starting to read mapping file. Takes a couple of minutes... ")
    with open(genbank_file_path) as gb_file:
        counter = 0
        while True:
            line = gb_file.readline()
            # If line is empty then end of file reached
            if not line:
                break
            line = line.split(delimiter)
            accessonAndVersion = line[1]
            taxid = line[2]
            if accessonAndVersion in accession_dict.keys():
                accession2taxid[accessonAndVersion] = taxid
            counter += 1
            if (counter % 10000000) == 0:
                print(f"        working on line #: {counter}")
    return accession2taxid

def file2correctHeader(acc2taxid, genome_directory):
    """
    modifes each genome in a directory, assumes is NCBI
    and adds taxid.
    """
    print("    Modifying each NCBI file in the genome specified directory")
    for filename in os.listdir(genome_directory):
        full_path = Path(genome_directory) / Path(filename)
        # READ HEADER AND FIND TAX ID
        with open(full_path, "r") as fasta_opened:
            file_lines = fasta_opened.readlines()
            accession_id = file_lines[0].split(" ")[0].strip(">")
            try:
                tax_id = acc2taxid[accession_id]
                print(f"Found for {tax_id}")
            except:
                print(f"cant find mapping for: {accession_id}")
                # os.remove(full_path) # DELETE FILE.
                continue

        os.remove(full_path) # DELETE FILE.

        # MAKE NEW FILE WITH TAXD IN HEADER
        with open(full_path, "w") as fasta_opened:
            header = file_lines[0].strip('\n')
            fasta_opened.write(f"{header} |kraken:taxid|{tax_id} \n")
            for line in file_lines[1:]:
                fasta_opened.write(line)

def main():
    print("Adding taxids to the NCBI fasta genomes.")

    if len(sys.argv) != 3:
        print(__doc__)
        exit(1)

    # input arguments
    genbank_accession_map = sys.argv[1]
    genome_directory = sys.argv[2]

    # create mapping dictionary
    acceptable_names = get_acc2tax_map(genome_directory)
    print("    Done!")
    acc2taxid = accession2taxid(genbank_accession_map, acceptable_names)
    print("    Done!")
    # modify each file in directory
    file2correctHeader(acc2taxid, genome_directory)
    print("    Done!")

if __name__ == "__main__":
    main()

