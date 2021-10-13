import sys
import re
import ntpath
from typing import List
from pathlib import Path

class DNA:
    def __init__(self, name, fasta_file, kmer_len):
        ''' Initializes dna object with name, reference genome file, and kmer size '''
        self.name: str = name
        self.fasta_file: Path = Path(fasta_file) if self.validate_file_extension(fasta_file) else None
        self.genome: str = self.fasta_to_genome(fasta_file) if self.fasta_file != None else None
        self.kmers: List = self.create_kmers(self.genome, kmer_len) if self.genome != None else []


    def create_kmers(self, genome, kmer_len = 20) -> List:
        ''' Generates k-mers given a dna sequence and specified k-mer length'''
        kmers = []
        g_len = len(genome)
        # assigns 1 kmer if the length passed is bigger than genome length
        if g_len < kmer_len or kmer_len <= 0:
            num_kmers = 1
            kmer_len = g_len
        else:
            num_kmers = g_len - kmer_len + 1

        for i in range(num_kmers):
            kmer = genome[i:i + kmer_len]
            kmers.append(kmer)

        return kmers

    # TODO: add multifasta logic
    def fasta_to_genome(self, fasta_path) -> str:
        ''' 
            DESCRIPTION:
                Extracts just the genome from the genome member variable, which should be a file path.
                Currently assumes only single FASTA file is passed (no multi-fasta). 
            
            INPUT:
                FASTA/FASTQ file (.fa, .fna, .fasta, .fastq)
            
            OUTPUT:
                Genome in string format
        '''
        sequence = ""
        with open(fasta_path) as fasta_file:
            fasta_lines = fasta_file.readlines()
            sequence_i = ""
            # If name is found (i.e. this is the correct file), get genome
            for counter, line in enumerate(fasta_lines):
                if line[0] == ">" and re.search(self.name, line, re.IGNORECASE):
                    if (counter != 0):
                        sequence += (sequence_i.rstrip())
                        sequence_i = ""
                elif line[0] != ">":
                    sequence += line.rstrip()

        return sequence


    def validate_file_extension(self, file) -> bool:
        if str(file).endswith((".fasta",".fa",".fna", ".fastq")):
            return True
        else:
            return False


    def calc_jaccard(self, dna2) -> float:
        ''' Calculates Jaccard similarity index between this DNA object and another '''
        a = set(self.kmers)
        b = set(dna2.kmers)

        intersection = len(a.intersection(b))
        union = len(a.union(b))

        return round(intersection / union, 8)


    def calc_minhash():
        ''' Calculates probabilistic jaccard using MinHash algorithm between
            two DNA objects '''
        raise NotImplimentedError("Minhash function should be implemented first.")


    def compare_results(jaccard, minhash) -> float:
        """
        Description:
            Compares the jaccard to the minhash output
        OUTPUT:
            float: aboslute difference
        TODO:
            Think of metric - percent diff.
        """
        return abs(jaccard - minhash)
    

    def __eq__(self, other):
        if isinstance(other, DNA):
            return (self.name == other.name and self.genome == other.genome)
        else:
            return False
    

    def __hash__(self):
        return hash((self.name, self.genome))

