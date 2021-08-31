import sys
sys.path.append(".")
from dna import DNA

# create objects
TEST_KMER_LEN = 3
test_dna = DNA("Blessica", "ACTGAATTTCG", TEST_KMER_LEN)
compare_dna = DNA("Ryadel", "ACTGTTTCCAG", TEST_KMER_LEN)

def test_init():
    ''' Test constructor '''
    dna1 = DNA("Phage1", "ZZZZZZZ", TEST_KMER_LEN)
    assert dna1.name == "Phage1"
    assert dna1.genome == "ZZZZZZZ"
    assert dna1.kmers == ['ZZZ', 'ZZZ', 'ZZZ', 'ZZZ', 'ZZZ']


def test_create_kmers():
    ''' Test kmer generation method '''
    expected = ['ACT', 'CTG', 'TGA', 'GAA', 'AAT', 'ATT', 'TTT', 'TTC', 'TCG']
    actual = test_dna.create_kmers(test_dna.genome, TEST_KMER_LEN)
    assert actual == expected


def test_calc_jaccard():
    ''' Test jaccard calculation '''
    # EXPECTED JACCARD VALUE SETUP
    set1 = set(test_dna.kmers)
    set2 = set(compare_dna.kmers)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    expected_jaccard_val = round(intersection / union, 4)
    #expected_jaccard_val = 0.2857

    assert test_dna.calc_jaccard(compare_dna) == expected_jaccard_val


def test_fasta_to_genome():
    ''' Test method that extracts a genome from a fasta file '''
    ''' test metrics: 1) partial genome comparison (first 100 chars)
                      2) char value at specific location
                      3) string length comparison '''
    expected_substring = "GTCTCCGAGCGATCTATCCACGACCAATTTGACATGGGTGCGCCGTTTGTAAAGGCCGTGGACAAAGCAGAACCCCCGGCACCGAGGGGGGCCGGGGGCCA"
    expected_val_at = expected_substring[50] # should be A if blessica
    expected_genome_len = 71470 # found by counting base pairs in blessica reference file

    blessica = DNA("Blessica", "/Users/latifa/Genomics/genomes/Blessica.fasta", 20)
    actual_genome = blessica.fasta_to_genome()
    actual_substring = blessica.genome[0:100]
    actual_val_at = actual_substring[50]
    actual_genome_len = len(blessica.genome)

    #assert actual_substring == expected_substring
    assert actual_val_at == expected_val_at
    assert expected_genome_len == actual_genome_len
