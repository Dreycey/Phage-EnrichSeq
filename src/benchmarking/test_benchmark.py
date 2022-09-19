import os
import pytest
from pathlib import Path

from benchmark import Benchmark, Truth, Result


benchmark_obj = Benchmark()
truth_obj = Truth(1, 'num_genomes', '3_genomes_500000')

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


def test_parse_truths_metadata():
    '''
    goal: make sure the function parses the 'truth' data from metadata csv correctly
    needs: CSV file containing metadata of experiment files
    cases: 
        
    '''
    pass



def test_parse_results_metadata():
    '''
    goal: make sure the function parses the results' data from metadata csv correctly
    needs: CSV file containing metadata of experiment files
    cases: 
        
    '''
    pass