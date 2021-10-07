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
    genome_ref_path = "/Users/latifa/GitHub/Phage-EnrichSeq/jaccard_module/"
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