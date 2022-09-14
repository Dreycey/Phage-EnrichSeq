import sys
import os
import argparse
import re
import csv
from pathlib import Path


TRUTH_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/'
RESULTS_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/'
OUTPUT_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/'


def parse_directory_trees(truth_path: str, result_path: str, output_path: str):
    '''
    DESCRIPTION:

    
    INPUT:

    OUTPUT:
        CSV file containing: <trial #>,<tool>,<test>,<subtest>,<path to final result file>
        e.g. 3,EnrichSeq,num_genomes,10_genomes_500000_reads,/projects/benchmarking/results/test_3/EnrichSeq/num_genomes/10_genomes_500000_reads/enrichseq/output_files/taxid_abundances.csv
    '''
    
    # create CSV file to append to
    output_csv_path = Path(output_path) / 'benchmarking_metadata.csv'

    with open(output_csv_path, 'a') as csvfile:
        writer = csv.writer(csvfile) 
        writer.writerow(['Trial No.', 'Experiment', 'Condition', 'Tool', 'Path to file'])

    # parse truth
    truth_rows = parse_truth_filepaths(truth_path)
    write_to_csv(truth_rows, output_csv_path)

    # parse results
    result_rows = parse_results_filepaths(result_path)
    write_to_csv(result_rows, output_csv_path)
    

def parse_truth_filepaths(dir_path: str):
    '''
    DESCRIPTION:
        Get metadata from file paths and write to CSV
        
    INPUT:
        Path to all truth files ('tests' in our case)
        Output CSV file path
    
    OUTPUT:
        creates CSV file
    '''

    rows = []

    for trial in os.listdir(dir_path):
        trial_num = int(''.join(filter(str.isdigit, trial)))
        trial_path = Path(dir_path) / trial
        if os.path.isdir(trial_path):
            for test in os.listdir(trial_path):
                test_path = Path(trial_path) / test
                if os.path.isdir(test_path) and 'config' not in test:
                    for subtest in os.listdir(test_path):
                        subtest_path = Path(test_path) / subtest
                        if os.path.isfile(subtest_path):
                            if subtest.endswith('.fa'):
                                info = [trial_num, test, subtest, 'truth', os.path.abspath(subtest_path)] 
                                rows.append(info)
                 
    return rows



def parse_results_filepaths(results_path: str):
    result_rows = []
    for trial in os.listdir(results_path):
        trial_num = int(''.join(filter(str.isdigit, trial)))
        trial_path = Path(results_path) / trial
        if os.path.isdir(trial_path):
            for tool in os.listdir(trial_path):
                tool_path = Path(trial_path) / tool
                if os.path.isdir(tool_path):
                    #result_rows.append(__parser_selector__(tool_path, trial_num))
                    for test in os.listdir(tool_path):
                        test_path = Path(tool_path) / test
                        if os.path.isdir(test_path) and 'config' not in test:
                            for subtest in os.listdir(test_path):
                                subtest_path = Path(test_path) / subtest
                                if os.path.isdir(subtest_path):
                                    # Now in tool's output directory
                                    result_info = [trial_num, test, subtest]
                                    tool_info = __tool_selector(tool, subtest_path)
                                    if tool_info: # if the info does not come back empty, append the row
                                        result_info.extend(tool_info)
                                        result_rows.append(result_info)
    return result_rows


def __tool_selector(tool_name: str, tool_result_path: Path):
    result_row = []
    if re.search('enrichseq', tool_name, re.IGNORECASE):
        full_result_path = Path(tool_result_path) / 'enrichseq/output_files/taxid_abundances.csv'
        if os.path.exists(full_result_path) and os.path.isfile(full_result_path):
            result_row = ['EnrichSeq', full_result_path]


    elif re.search('fastviromeexplorer', tool_name, re.IGNORECASE):
        full_result_path = Path(tool_result_path) / 'FastViromeExplorer-final-sorted-abundance.tsv'
        if os.path.exists(full_result_path) and os.path.isfile(full_result_path):
            result_row = ['FastViromeExplorer', full_result_path]

    elif re.search('bracken', tool_name, re.IGNORECASE):
        full_result_path = Path(tool_result_path) / 'abundances.bracken'
        if os.path.exists(full_result_path) and os.path.isfile(full_result_path):
            result_row = ['FastViromeExplorer', full_result_path]

    return result_row


def write_to_csv(rows: list, output_file: str):
    with open(output_file, 'a') as csvfile: 
        writer = csv.writer(csvfile) 
        writer.writerows(rows)



def main():
    parse_directory_trees(TRUTH_PATH, RESULTS_PATH, OUTPUT_PATH)


if __name__ == "__main__":
    main()