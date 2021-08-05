import sys
sys.path.append(".")
from dna import DNA

# create objects

test_dna = DNA("Blessica", "ACTGAATTTCG")
compare_dna = DNA("Ryadel", "ACTGTTTCCAG")

def test_init():
    ''' Test constructor '''
    dna1 = DNA("Phage1", "ZZZZZZZ")
    assert dna1.name == "Phage1"
    assert dna1.genome == "ZZZZZZZ"
    assert dna1.kmers == ['ZZZ', 'ZZZ', 'ZZZ', 'ZZZ', 'ZZZ']


def test_create_kmers():
    ''' Test kmer generation method '''
    expected = ['ACT', 'CTG', 'TGA', 'GAA', 'AAT', 'ATT', 'TTT', 'TTC', 'TCG']
    actual = test_dna.create_kmers(test_dna.genome, 3)
    assert actual == expected


def test_calc_jaccard():
    ''' Test jaccard calculation '''
    jaccard_val = 0.2857
    assert test_dna.calc_jaccard(compare_dna) == jaccard_val


def test_fasta_to_genome():
    ''' Test method that extracts a genome from a fasta file '''
    # TODO: compare genome values 
