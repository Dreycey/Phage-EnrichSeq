"""
This script makes config files that can be used by the read simulator
for simulating metagenomic mixtures.

Config file output:
    # <genome path> <relative abundance>
    ../genomes/Blessica.fasta, 0.15

Usage:
    python config_maker.py <path to fasta directory> <genome number>
Example:
    python config_maker.py ../../../database/ref_genomes/ 2000
"""
import os
import sys
import random




class GenomePathObject:
    """ data structure to hold/add/get genome paths """

    def __init__(self):
        self.genomes = []

    def get_paths(self):
        """
        Description:                                                            
            This command              
        Input (parameters):                                                     
            1. N - number of genomes to grab.                                   
        Output (returns):                                                       
            1. List with genome paths
        """
        return self.genomes

    def add(self, genome_path):
        """
        Description:                                                            
            This command adds a genome to the interal data structure                                           
        Input (parameters):                                                     
            1. N - number of genomes to grab.                                   
        Output (returns):                                                       
            1. void
        """
        self.genomes.append(genome_path)
    
class ConfigMaker:

    def __init__(self, genomes_directory):
        self.genomes_directory = genomes_directory

    def choose_N_genomes(self, N, randomize:bool=True) -> GenomePathObject:
        """
        Description:
            This command retrieves N random genome paths from the DB.
        Input (parameters): 
            1. N - number of genomes to grab.
            2. randomize [DEFAULT = True] - randomize genomes choosen.
        Output (returns):
            1. 
        """
        genomes = GenomePathObject() 
        list_of_genomes = os.listdir(self.genomes_directory)
        if randomize: random.shuffle(list_of_genomes)
        if N > len(list_of_genomes):
            print(f"requested {N} genomes but only {len(list_of_genomes)} in DB")
            exit(1)
        i, j = 0, 0
        while (j < N):
            if i >= len(list_of_genomes):                                       
                print(f"Max number of genome paths without spaces is: {j} \n")  
                exit(1) 
            genome_path = list_of_genomes[i]
            if not (' ' in genome_path): # no spaces allowed
                genomes.add(genome_path)
                j += 1
            i += 1
        return genomes

    def create_config(self, N):
        """
        Description:
            creates a config file from input parameters.
        Input (parameters): 
            1. N - number of genomes to grab.
        Output (returns):
            1. returns path to the config file.
        """
        n_genomes = self.choose_N_genomes(N)
        for genome in n_genomes.get_paths():
            print(f"{self.genomes_directory}{genome}, {round(1/N, 5)}")

def main():
    # configuration
    if len(sys.argv) == 3:
        genome_directory = sys.argv[1]
        number_of_genomes = int(sys.argv[2])
    else:
        print("\nERROR: input not correct. refer to the below usage: \n\n")
        print(__doc__)
        exit(1)
    # build config
    config = ConfigMaker(genome_directory)
    config.create_config(number_of_genomes)
    
if __name__ == "__main__":
    main()
