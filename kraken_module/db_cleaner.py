#! usr/bin/python3
"""
DESCRIPTION:
    The purpose of this script is to clean the reference genome
    database.
USAGE:
    python db_cleaner.py <genome directory>
    python db_cleaner.py genome_directory/ 
TODO:
    TIME COMPLEXITY: O(N+m)
    1. Obtain mapping from taxid2pathAndLength -> {XXXXX : [Path(), 22222]} T.C.=O(N)
        1.5 Per taxid, store path to longest genome.
    2. Create hash for Paths. O(m); where m=# of taxids
    3. During deletion, if Path not in hash, delete. T.C.=O(N)
"""
# STD
import os
import sys
from pathlib import Path
from typing import Dict, List
from tqdm import tqdm
import multiprocessing as mp
# in house packages
current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f"{current_path}/../PathOrganizer_module")
from PathOrganizer import PathOrganizer, PathErrors, DuplicateGenomeError




class DB_Duplicate_Cleaner(PathOrganizer):
    """ This class deletes a file if it's tax id has bee seen already"""

    def get_genome_length(self, fasta_path_in):
        """
        given the path to a fasta, this function gets the 
        genome's length
        """
        genome_seq = ""
        with open(fasta_path_in) as fasta:
            for line_index, line in enumerate(fasta.readlines()):
                if line[0] == ">" and (line_index == 0):
                    continue
                elif (line[0] == ">") and (line_index > 0):
                    break
                else:
                    genome_seq += line.strip("\n")
        return len(genome_seq)


def main():
    # INPUT ARGUMENTS
    genome_dir_path = sys.argv[1]

    # get the paths to all genomes
    full_genome_paths = [Path(genome_dir_path) / Path(fasta_file)
                         for fasta_file in os.listdir(genome_dir_path)]

    # obtain genomes to keep (keep the ones with longest genome.)
    print("FINDING FIILES TO KEEP")
    db_cleaner = DB_Duplicate_Cleaner(genome_dir_path)
    taxid2PathAndLength: Dict[int, List[str]] = {}
    for genome_path in tqdm(full_genome_paths): # O(N)
        tax_id = int(db_cleaner.get_fasta_taxid(genome_path))
        genome_length =  db_cleaner.get_genome_length(genome_path)
        if tax_id in taxid2PathAndLength.keys():
            original_length = taxid2PathAndLength[tax_id][1]
            if genome_length > original_length:
                taxid2PathAndLength[tax_id] = [genome_path, genome_length]
        else:
            taxid2PathAndLength[tax_id] = [genome_path, genome_length]

    # get paths to keep 
    files_to_keep = {}
    for taxid, pathAndLength in taxid2PathAndLength.items():
        path_of_genome = pathAndLength[0]
        files_to_keep[path_of_genome] = 1

    # delete files not in files_to_keep
    print("DELETNG DUPLICATE FILES")
    duplicate_count = 0
    for genome_path in tqdm(full_genome_paths):
        if genome_path not in files_to_keep.keys():
            duplicate_count += 1
            os.remove(f"{genome_path}")

    print(f"SUCCESS! Deleted {duplicate_count} duplicate(s)")

if __name__ == "__main__":
    main()