"""
Mutations
"""
import sys
import random




def openFasta(fasta_path): 
    """
    Description: opens a fasta file
    Input: str - path to fasta file containing tail fibers
    """
    fastaName2Seq = {}
    with open(fasta_path, "r") as opened_fasta:
        fasta_line = opened_fasta.readline()
        line_counter = 0
        while fasta_line:
            if fasta_line[0] == ">":
                if line_counter > 0: fastaName2Seq[fasta_name] = fasta_seq
                fasta_name = fasta_line.strip("\n").strip(">")
                fasta_seq = ""
            else:
                fasta_seq += fasta_line.strip("\n").strip(" ")
            fastaName2Seq[fasta_name] = fasta_seq
            fasta_line = opened_fasta.readline()
            line_counter += 0

    return fastaName2Seq

def genome_mutate(genome_path, ani):
    """ creates a genome with an ANI """
    seq2genome = openFasta(genome_path)
    genome = list(seq2genome.values())[0]
    mutated_genome = [i for i in genome]
    genome = list(seq2genome.values())[0]
    mutation_dictionary = {}
    alphabet = {'A' : 'T', 'C' : 'G', 'T' : 'A', 'G' : 'T'}
    num_edits = int(len(genome) * (1-ani))
    edits_added = 0
    while (edits_added < num_edits):
        index2edit = random.choice(range(0, len(genome)))
        if index2edit not in mutation_dictionary:
            curr_nucleotide = genome[index2edit]
            # print(f"for index: {index2edit}")
            # print(f"before: {curr_nucleotide}")
            new_nucleotide = random.choice([i for i in alphabet.keys() \
                                          if i != curr_nucleotide])
            # print(f"After: {new_nucleotide}")
            mutated_genome[index2edit] = new_nucleotide
            edits_added += 1
            mutation_dictionary[index2edit] = new_nucleotide
        else:
            continue
    return ''.join(mutated_genome)

def print_genome(genome: str):
    """ prints genome with 80 lines each """
    print(">misc_" + str(random.choice(range(0, 1000))))
    for index in range(0, len(genome),80):
        if (index+80) < len(genome):
            print(genome[index:index+80])
        else:
            print(genome[index:len(genome)])


if __name__ == "__main__":
    genome_path = sys.argv[1]
    ani = float(sys.argv[2])
    mutated_genome = genome_mutate(genome_path, ani)
    print_genome(mutated_genome)
    print()