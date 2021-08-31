import sys
import re
import ntpath

class DNA:
    def __init__(self, name, fasta_file, kmer_len):
        ''' Initializes dna object with name, reference genome file, and kmer size '''
        self.name: str = name
        
        # validate file input
        if self.validate_file_extension(fasta_file):
            self.genome: str = self.fasta_to_genome(fasta_file)
            self.kmers: List = self.create_kmers(self.genome, kmer_len)
        else:
            return


    def set_genome(self, genome_str, kmer_len):
        self.genome = genome_str
        self.kmers: List = self.create_kmers(self.genome, kmer_len)


    def create_kmers(self, genome, kmer_len = 20):
        ''' Generates k-mers given a dna sequence and specified k-mer length'''
        kmers = []
        num_kmers = len(genome) - kmer_len + 1

        for i in range(num_kmers):
            kmer = genome[i:i + kmer_len]
            kmers.append(kmer)

        return kmers


    def fasta_to_genome(self, fasta_file):
        ''' Extracts just the genome from the genome member variable, which should be a file path
            Assumptions: 1) File name is the genome name
                         2) Only single FASTA file is passed (no multi-fasta) '''

        # FIRST: get genome name from file path
        file_name = ntpath.basename(fasta_file)
        # file_name = self.validate_file_extension(ntpath.basename(fasta_file))


        # SECOND: search file for that name and get genome
        fasta_lines = open(fasta_file).readlines()
        genome = ""
        # Make sure desired phage name is in the file. If not, assign null value to genome variable
        # if ">" in fasta_lines[0] and not re.search(self.name, fasta_lines[0], re.IGNORECASE):
        #     self.genome = None
        # If name is found (i.e. this is the correct file), get genome
        if ">" in fasta_lines[0] and re.search(self.name, fasta_lines[0], re.IGNORECASE):
            #self.genome = "" # clear current value for genome
            for line in fasta_lines:
                if ">" not in line:
                    genome += line.rstrip()

        return genome


    def validate_file_extension(self, file):
        if not file.endswith((".fasta",".fa",".fna", ".fastq")):
            raise ValueError("Must be a FASTA or FASTQ file")
            return False
        else:
            return True


    def calc_jaccard(self, dna2):
        ''' Calculates Jaccard similarity index between this DNA object and another '''
        a = set(self.kmers)
        b = set(dna2.kmers)

        intersection = len(a.intersection(b))
        union = len(a.union(b))

        return round(intersection / union, 4)


    def calc_minhash():
        ''' Calculates probabilistic jaccard using MinHash algorithm between
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
