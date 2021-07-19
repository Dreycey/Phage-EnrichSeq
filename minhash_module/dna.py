import sys
import re

class DNA:
    name = ""
    genome = ""  # string or FASTA file?
    kmers = []


    def __init__(self, name, fasta_file):
        '''  '''
        self.name = name
        #self.genome = self.fasta_to_genome(fasta_file) # TODO: extract from FASTA file
        self.genome = fasta_file
        self.kmers = self.create_kmers(self.genome, 3)

    def create_kmers(self, genome, kmer_len = 20):
        ''' generates k-mers given a dna sequence and specified k-mer length'''
        kmers = []
        num_kmers = len(genome) - kmer_len + 1

        for i in range(num_kmers):
            kmer = genome[i:i + kmer_len]
            kmers.append(kmer)

        return kmers


    def fasta_to_genome(self):
        ''' extracts the genome from a fasta file '''



    def print(self):
        print(f'Name: {self.name} \nGenome: {self.genome} \n#k-mers generated: {len(self.kmers)}')


    def fasta_to_genome(file):
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



def calc_jaccard(dna1, dna2):
    ''' calculates Jaccard similarity index between two DNA objects '''
    a = set(dna1.kmers)
    b = set(dna2.kmers)

    intersection = len(a.intersection(b))
    union = len(a.union(b))

    return intersection / union



def calc_minhash():
    ''' calculates probabilistic jaccard using MinHash algorithm between
        two DNA objects '''
    return 0


def compare_results(jaccard, minhash):
    return abs(jaccard - minhash)



def main():
    ''' testing dna class methods '''

    # Test extraction of genomes from files
    # string1 = ">Mycobacterium Blessica alsdkflsdk"
    # substring1 = "blah"
    # if re.search(substring1, string1, re.IGNORECASE):
    #     print("FOUND")
    # else:
    #     print("NOT FOUND")

    # Get genome from fasta files (how many files?)
    dna1 = DNA("Blessica", "ACTGAATTTCG")
    dna2 = DNA("Ryadel", "ACTGTTTCCAG")
    dna3 = DNA("D29", "AAAACCCCTTTTGGG")
    file1 = sys.argv[1]
    #file2 = sys.argv[2]



    # Test k-mer building
    #dna1.print()
    print(dna1.kmers)
    print(dna2.kmers)
    print(dna3.kmers)

    # Test Jaccard similarity
    print(format(calc_jaccard(dna1,dna2),".4f"))
    print(format(calc_jaccard(dna1,dna3),".4f"))
    print(format(calc_jaccard(dna2,dna3),".4f"))

    # Test MinHash algorithm


    # Compare Jaccard and MinHash

if __name__ == "__main__":
    main()
