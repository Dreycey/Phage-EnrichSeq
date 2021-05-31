class DNA:
    name = ""
    genome = ""  # string or FASTA file?
    kmers = []


    def __init__(self, name, genome):
        self.name = name
        self.genome = genome # TODO: extract from FASTA file
        self.kmers = self.create_kmers(self.genome)

    def create_kmers(self, genome, kmer_len = 5):
        ''' generates k-mers given a dna sequence and specified k-mer length'''
        kmers = []
        num_kmers = len(genome) - kmer_len + 1

        for i in range(num_kmers):
            kmer = genome[i:i + kmer_len]
            kmers.append(kmer)

        return kmers


    def print(self):
        print(f'Name: {self.name} \nGenome: {self.genome} \n#k-mers generated: {len(self.kmers)}')



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


def extract_genome(file, name):
    ''' extracts just the genome from a given FASTA file '''
    fasta_lines = open(file).readlines()
    # search for name in file


def __main__():
    ''' testing dna class methods '''
    dna1 = DNA("Blessica", "ACTGAATTTCG")
    dna2 = DNA("Ryadel", "ACTGTTTCCAG")

    # Test k-mer building
    #dna1.print()
    print(dna1.kmers)
    print(dna2.kmers)

    # Test Jaccard similarity
    print(calc_jaccard(dna1,dna2))

    # Test MinHash algorithm


    # Compare Jaccard and MinHash

__main__()
