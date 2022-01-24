"""
DESCRIPTION:
    This script contains the functional tests for the MergeOverlap module.
"""
# std packages
from enum import Enum
import pickle
import os
# non-std packages
import pytest
from pathlib import Path
# in house packages
from mergeoverlap import GenomeTestSet




class singleGenomeTest(Enum):
    #directory_path = "./mergeoverlap_filter_module/"
    directory_path = os.path.dirname(os.path.realpath(__file__))
    input = f"{directory_path}/single_test_files/taxid_file.lsv"
    fasta = f"{directory_path}/single_test_files/paired_illumina_2.fa"
    genome_directory = f"{directory_path}/test_genome_dir"
    output_prefix = f"{directory_path}/testing_output/single_genome"
    truth = [("1076136", 1.0)]

class multiGenomeTest(Enum):
    #directory_path = "./mergeoverlap_filter_module/"
    directory_path = os.path.dirname(os.path.realpath(__file__))
    input = f"{directory_path}/multi_test_files/taxid_file.lsv"
    fasta = f"{directory_path}/multi_test_files/3_genomes_sim_10000_illumina.fa"
    genome_directory = f"{directory_path}/test_genome_dir"
    output_prefix = f"{directory_path}/testing_output/multi_genome"
    truth = [("2886930", 0.3),  ("10868", 0.3), ("2681618", 0.3)]

def parse_output_csv(output_csv, delimiter=","):
    """ parses the output csv """
    results_dict = {}
    output_csv = open(output_csv, "r")
    line = output_csv.readline()
    while line:
        line = line.split(delimiter)
        results_dict[line[0]] = float(line[1])
        line = output_csv.readline()
    return results_dict

def test_singlegenome():
    """
    This tests a single genome.
    """
    genomeTestObj = GenomeTestSet(line_seperated_genomes=singleGenomeTest.input.value, 
                                  genome_directory=singleGenomeTest.genome_directory.value)
    genomeTestObj.checkSeqFile(singleGenomeTest.fasta.value)
    genomeTestObj.saveResultAsCSV(singleGenomeTest.output_prefix.value+".csv")
    results = parse_output_csv(singleGenomeTest.output_prefix.value+".csv")

    for taxid, abundance in singleGenomeTest.truth.value:
        assert abs(results[taxid] - abundance) < 0.10

    
def test_multiplegenomes():
    """
    This tests a file with multiple genomes
    """
    genomeTestObj = GenomeTestSet(line_seperated_genomes=multiGenomeTest.input.value, 
                                  genome_directory=multiGenomeTest.genome_directory.value)
    genomeTestObj.checkSeqFile(multiGenomeTest.fasta.value)
    genomeTestObj.saveResultAsCSV(multiGenomeTest.output_prefix.value+".csv")
    results = parse_output_csv(multiGenomeTest.output_prefix.value+".csv")

    for taxid, abundance in multiGenomeTest.truth.value:
        assert abs(results[taxid] - abundance) < 0.10
