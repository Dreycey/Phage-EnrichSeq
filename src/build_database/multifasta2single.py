"""
DESCRIPTION:
    This script converts a multi fasta into many individual 
    fasta files. 
USAGE:
    python multifasta2single.py <multifasta file> <output directory>
EXAMPLE:
    python multifasta2single.py ref_genomes/actinoReformatted.fa ref_genomes/
"""
import sys
from pathlib import Path
from tqdm import tqdm




def multiFasta2fasta(multifasta_path):
    """
    DESCRIPTION - parses a fasta or multifasta file.
    INPUT - 1. fasta path
    OUTPUT - 1. name (ordered list); 2. genome sequence (ordered list)
    """
    seq_names, sequences = [], []
    with open(multifasta_path) as  multi_fasta_file:
        fasta_input =  multi_fasta_file.readlines()
        sequence_i = ""
        for counter, line in enumerate(fasta_input):
            if (line[0] == ">"):
                seq_name = line[1:]
                seq_names.append(seq_name.strip("\n"))
                if (counter != 0):
                    sequences.append(sequence_i.strip("\n"))
                    sequence_i = "" # renew sequence
            else:
                sequence_i += line.strip("\n")
        sequences.append(sequence_i.strip("\n"))
    return seq_names, sequences

def addSeqsToFiles(output_directory, seq_names, sequences):
    """
    DESCRIPTION: takes sequence names, sequences, and add to output directory
    """
    for seq_index, sequence_name in tqdm(enumerate(seq_names)):
        full_path = Path(output_directory) / Path(sequence_name)
        with open(full_path.with_suffix(".fna"), "w") as fasta_out:
            fasta_out.write(">"+sequence_name + "\n" + sequences[seq_index])

def main():
    if (len(sys.argv) != 3):
        print(__doc__)
        exit(1)
    else:
        multifasta_path = sys.argv[1]
        output_directory = sys.argv[2]
        seq_names, sequences = multiFasta2fasta(multifasta_path)
        addSeqsToFiles(output_directory, seq_names, sequences)

if __name__ == "__main__":
    main()