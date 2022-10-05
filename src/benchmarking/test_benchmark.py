#!/usr/local/bin/python
import os
import sys
import pytest
import pandas as pd
from pathlib import Path

from benchmark import Benchmarking, SimulatedTruth, Result, EnrichSeqResult


# benchmark_obj = Benchmarking()
# simulated_truth_obj = SimulatedTruth(1, 'num_genomes', '3_genomes_500000')
# enrichseq_obj = EnrichSeqResult()


@pytest.fixture
def all_tools_correct_csv():
    return '/Users/latifa/GitHub/benchmarking-enrichseq/metadata_all_correct.csv'

@pytest.fixture
def missing_truth_csv():
    return '/Users/latifa/GitHub/benchmarking-enrichseq/metadata_missing_truth.csv'


def test_extract_genomes_and_reads():
    '''
    needs: string in a specific format
    cases:  2 exact same format, different length numbers, 
            1 with same keywords but different order,
            1 with one of the keywords, 
            1 with completely incorrect format
    '''

    pass


def test_parse_simulated_fasta():
    pass


def test_calc_precision():
    '''
    needs:  true positives, false positives, true negatives, false negatives (from )
    cases:  
    '''
    pass

def test_calc_recall():
    pass

def test_calc_f1score():
    pass

def test_calc_l2distance():
    pass


def test_parse_enrichseq_results():
    pass


def test_parse_fve_results():
    pass


def test_parse_bracken_results():
    pass


def test_parse_metadata_csv(all_tools_correct_csv):
    benchmark_obj = Benchmarking(all_tools_correct_csv)
    expected_data = [[2, 'num_reads_200genomes', '200_genomes_400000_reads', 'truth', '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/test_2/num_reads_200genomes/200_genomes_400000_reads.fa'],
                     [2, 'num_reads_200genomes', '200_genomes_400000_reads', 'EnrichSeq', '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/test_2/EnrichSeq/num_reads_200genomes/200_genomes_400000_reads/enrichseq/output_files/taxid_abundances.csv'],
                     [2, 'num_reads_200genomes', '200_genomes_400000_reads', 'FastViromeExplorer', '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/test_2/FastViromeExplorer/num_reads_200genomes/200_genomes_400000_reads/FastViromeExplorer-final-sorted-abundance.tsv']]
    expected_df = pd.DataFrame(expected_data, columns=['Trial_Num', 'Experiment', 'Condition', 'Tool', 'File_Path'])
    actual_df = benchmark_obj.parse_metadata_csv()
    pd.testing.assert_frame_equal(expected_df,actual_df)



def test_parse_truths_metadata():
    '''
    goal: make sure the function parses the 'truth' data from metadata csv correctly
    needs: CSV file containing metadata of experiment files
    cases: 
        
    '''
    # assert SimulatedTruth objects in lists have the same values
    pass



def test_parse_results_metadata():
    '''
    goal: make sure the function parses the results' data from metadata csv correctly
    needs: CSV file containing metadata of experiment files
    cases:  1) all fields are present and correctly formatted
        
    '''
    
    pass


def test_assign_truth_objs():
    '''
    cases:  1) all of result's metadata matches truth's metadata (assert == 1.0)
            2) 1/2 of given results don't match truths based on trial num (assert == 0.5)
            3) 1/2 of given results don't match truths based on experiment name (assert == 0.5)
            4) 1/2 of given results don't match truths based on experiment condition (assert == 0.5)
            5) none match based on any or all columns (assert == 0.0)

    '''
    pass