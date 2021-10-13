import sys
import pytest
sys.path.append(".")
from dna import DNA
from genome_comparison import GenomeCompare
from pathlib import Path
from pytest import mark


''' Pytest fixture setup '''
@pytest.fixture
def default_kmer_len():
    default_kmer_len = 20
    return default_kmer_len


@pytest.fixture
def genome_ref_path():
    genome_ref_path = "/Users/latifa/GitHub/Phage-EnrichSeq/genomeCompare_module/"
    return genome_ref_path


@pytest.fixture
def genome_compare_obj(default_kmer_len, genome_ref_path):
    ''' set up a GenomeCompare object that has 3 DNA objects '''
    blessica = DNA("Blessica", genome_ref_path + "Blessica.fasta", default_kmer_len)
    d29 = DNA("D29", genome_ref_path + "D29.fasta", default_kmer_len)
    ryadel = DNA("Ryadel", genome_ref_path + "Ryadel.fasta", default_kmer_len)
    phages = [blessica, d29, ryadel]

    genome_compare_obj = GenomeCompare(phages)
    return genome_compare_obj


def test_init():
    # TODO: check that dna objects are added to the dnaList
    pass


def test_create_adjacency_matrix(genome_compare_obj):
    ''' TODO:   1) check expected matrix values (with 3 phages)
                2) what happens when there is some error with dna object? (e.g. no kmers/genome)
    '''
    expected = [[],[]]
    actual = genome_compare_obj.create_adjacency_matrix()
    assert expected == actual


def test_prune_adj_matrix(genome_ref_path, default_kmer_len, genome_compare_obj):
    ''' TODO:   1) check that clusters returned contain the expected dna objects '''
    dict_ = {"C_1": [DNA("Blessica", genome_ref_path + "Blessica.fasta", default_kmer_len), DNA("Ryadel", genome_ref_path + "Ryadel.fasta", default_kmer_len)]}
    assert genome_compare_obj.prune_adj_matrix(0.5) == dict_