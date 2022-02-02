#! usr/bin/python3
"""
TESTING.

This file contains both unit testing and regression testing
for the read simulator portion of EnrichSeq. 
"""
# in-house modules/libs
from pathlib import Path
from functional_test_simulate_reads import GenomeTestSet, MinimapMapper
# std pckgs
import os
import sys

# Global Variables
PATH = os.path.dirname(os.path.abspath(__file__))
TEST_CONFIG = Path(PATH) / Path("simulate_genomes.config")
FASTA_FILE = Path(PATH) / Path("simulated_test_reads_illumina.fa")
TEST_FASTA = Path(PATH) / Path("test_fasta.fa")
# # Changing to directory of script.
os.chdir(PATH)

config_object = GenomeTestSet(TEST_CONFIG)
def test_functional_output(data_regression):
    config_object.checkSeqFile(FASTA_FILE)
    data_regression.check(config_object.resultDict)

def test_parseFasta():
    seq_names, sequences = config_object.parseFasta(TEST_FASTA)
    assert len(seq_names) == 2
    assert seq_names[0] == "name 1"
    assert seq_names[1] == "name 2"
    assert sequences[0] == "ATCGATCGATCGTAGCTAGCTGACTGATCGATGC"
    assert sequences[1] == "ATGCATGCATGCTC"