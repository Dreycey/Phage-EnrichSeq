import sys
import pytest
import numpy as np
sys.path.append(".")
from dna import DNA
from genome_comparison import GenomeCompare
from pathlib import Path
from pytest import mark


''' Global variables '''
GENOME_REF_PATH = "/Users/latifa/GitHub/Phage-EnrichSeq/genomeCompare_module/"

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


''' TEST LIST CREATION '''
def test_list_one():
    return [DNA("A", GENOME_REF_PATH + "A.fasta", 3), 
            DNA("B", GENOME_REF_PATH + "B.fasta", 3),
            DNA("C", GENOME_REF_PATH + "C.fasta", 3),
            DNA("D", GENOME_REF_PATH + "D.fasta", 3),
            DNA("E", GENOME_REF_PATH + "E.fa", 3)]

def test_list_two():
    return [DNA("A", GENOME_REF_PATH + "A.fasta", 3), 
            DNA("B", GENOME_REF_PATH + "B.fasta", 3),
            DNA("C", GENOME_REF_PATH + "C.fasta", 3)]


def test_list_three():
    return [DNA("C", GENOME_REF_PATH + "C.fasta", 7), 
            DNA("D", GENOME_REF_PATH + "D.fasta", 7)]


def test_init():
    # TODO: check that dna objects are added to the dnaList
    pass


@pytest.mark.parametrize("dnaList, expected", 
                        [
                            (test_list_one(), np.array([[1.0, 0.91891892, 0.4, 0.375, 0.41666667],
                                                       [0.91891892, 1.0, 0.39583333, 0.36956522, 0.38297872],
                                                        [0.4, 0.39583333, 1.0, 0.37777778, 0.33333333],
                                                        [0.375, 0.36956522, 0.37777778, 1.0, 0.76470588],
                                                        [0.41666667, 0.38297872, 0.33333333, 0.76470588, 1.0]])), 
                            (test_list_two(), np.array([[1.000000,0.91891892,0.4],[0.91891892,1.000000,0.39583333],[0.4,0.39583333,1.000000]])),
                            (test_list_three(), np.array([[1.0,0.0],[0.0,1.0]]))
                        ])
def test_create_adjacency_matrix(dnaList, expected):
    ''' TODO:   1) check expected matrix values (with 3 phages)
                2) what happens when there is some error with dna object? (e.g. no kmers/genome)
    '''
    # create GenomeCompare object out of the parametrized dnaList
    genComp = GenomeCompare(dnaList)
    actual = genComp.create_adjacency_matrix()
    assert np.array_equal(expected, actual)

# assert lengths of sets
# TODO: parametrize ??
# def test_prune_adj_matrix(genome_ref_path, default_kmer_len, genome_compare_obj):
#     ''' TODO:   1) check that clusters returned contain the expected dna objects '''
#     dict_ = {"C_1": [DNA("Blessica", genome_ref_path + "Blessica.fasta", default_kmer_len), DNA("Ryadel", genome_ref_path + "Ryadel.fasta", default_kmer_len)]}
#     assert genome_compare_obj.prune_adj_matrix(0.5) == dict_