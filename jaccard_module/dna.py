import sys
import re

class DNA:
    def __init__(self, name, fasta_file, kmer_len):
        '''  '''
        self.name: str = name
        self.genome: str = fasta_file
        self.kmers: List = self.create_kmers(self.genome, kmer_len)

    def create_kmers(self, genome, kmer_len = 20):
        ''' generates k-mers given a dna sequence and specified k-mer length'''
        kmers = []
        num_kmers = len(genome) - kmer_len + 1

        for i in range(num_kmers):
            kmer = genome[i:i + kmer_len]
            kmers.append(kmer)

        return kmers


    def fasta_to_genome(self, file):
        ''' extracts just the genome from a given FASTA file
            Assumptions: 1) File name is the genome name
                         2) Only single FASTA file is passed (no multi-fasta) '''

        # FIRST: get genome name from file path

        # SECOND: search file for that name and get genome

        fasta_lines = open(file).readlines()
        i=0
        # search for name in file
        # for line in fasta_lines:
        #     if ">" in line and re.search(name, line, re.IGNORECASE):
        #         print(f"Line #{}{name}")
        for line in fasta_lines:
            if ">" in line and re.search(name, line, re.IGNORECASE):
                genome = name+"_sequence"
                break


    def calc_jaccard(self, dna2):
        ''' calculates Jaccard similarity index between this DNA object and another '''
        a = set(self.kmers)
        b = set(dna2.kmers)

        intersection = len(a.intersection(b))
        union = len(a.union(b))

        return round(intersection / union, 4)


    def calc_minhash():
        ''' calculates probabilistic jaccard using MinHash algorithm between
            two DNA objects '''
        return 0


    def compare_results(jaccard, minhash):
        """
        Description:
            Compares the jaccard to the minhash output
        Output:
            float: aboslute difference
        TODO:
            Think of metric - percent diff.
        """
        return abs(jaccard - minhash)


    def print(self):
        print(f'Name: {self.name} \nGenome: {self.genome} \n#k-mers generated: {len(self.kmers)}')
