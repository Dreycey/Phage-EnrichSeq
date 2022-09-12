import sys
import os
import argparse
import glob
import csv
from pathlib import Path


TRUTH_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/'
RESULTS_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/'
OUTPUT_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/'


def parse_directory_trees(dir_path: str, output_path: str):
    '''
    DESCRIPTION:

    
    INPUT:

    OUTPUT:
        CSV file containing: <trial #>,<tool>,<test>,<subtest>,<path to final result file>
        e.g. 3,EnrichSeq,num_genomes,10_genomes_500000_reads,/projects/benchmarking/results/test_3/EnrichSeq/num_genomes/10_genomes_500000_reads/enrichseq/output_files/taxid_abundances.csv
    '''
    
    # create CSV file to append to
    output_csv_path = Path(output_path) / 'benchmarking_metadata.csv'
    #output_csv = open(Path(output_path) / output_csv_filename, 'a')

    # parse truth
    parse_truth_filepaths(dir_path, output_csv_path)

    # parse results
    


def parse_truth_filepaths(dir_path: str, output_file: str):
    '''
    DESCRIPTION:
        Get metadata from file paths and write to CSV
        
    INPUT:
        Path to all truth files ('tests' in our case)
        Output CSV file path
    
    OUTPUT:
        CSV file
    '''

    ## TODO: figure out what to do with config_files. Exclude from CSV now or later?
    truth_paths_dict = get_truth_filepaths(dir_path)
    with open(output_file, 'a') as csvfile: 
        writer = csv.writer(csvfile) 
        for trial in os.listdir(dir_path):
            row = []
            trial_num = int(''.join(filter(str.isdigit, trial)))
            row.append(trial_num)
            trial_path = Path(dir_path) / trial
            if os.path.isdir(trial_path):
                row.append('truth')
                for test in os.listdir(trial_path):
                    test_path = Path(trial_path) / test
                    if os.path.isdir(test_path) and 'config' not in test:
                        row.append(test)
                        for subtest in os.listdir(test_path):
                            subtest_path = Path(test_path) / subtest
                            if os.path.isfile(subtest_path):
                                if subtest.endswith('.fa'):
                                    row.append(truth_paths_dict[subtest])
            print(row)                    
            writer.writerow(row)
    return row



def get_truth_filepaths(dir_path: str) -> dict:
    '''
    DESCRIPTION:
        Recursively finds the truth files given a parent directory, then
        stores the path in a dictionary
    
    INPUT:
        Path to truth directory
    
    OUTPUT:
        Dictionary in the format [file name]:[file path]
    '''
    truth_dict = {}

    for path in Path(dir_path).rglob('*.fa'):
        truth_dict[path.name] = str(path)

    return truth_dict


def parse_enrichSeq_structure():
    pass


def parse_fastviromeexplorer_structure():
    pass


def parse_bracken_structure():
    pass


def parse_args(argv=None) -> argparse.Namespace:
    '''
    DESCRIPTION:
        This method takes in the arguments from the command and performs
        parsing.

    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-t", "--truth_dir", help="path to truth directory", required=True)
    parser.add_argument("-r", "--result_dir", help="path to results directory", required=True)
    parser.add_argument("-o", "--output_dir", help="output CSV path", required=True)
    return parser.parse_args(argv)


def main():
    arguments = parse_args(argv=sys.argv[1:])
    parse_directory_trees(TRUTH_PATH, OUTPUT_PATH)


if __name__ == "__main__":
    main()