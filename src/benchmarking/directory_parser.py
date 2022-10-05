import os
import sys
import argparse
import csv
from pathlib import Path
from enum import Enum
from dynaconf import Dynaconf

TRUTH_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/'
RESULTS_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/'
OUTPUT_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/'

class ToolPathFinder(Enum):
    ENRICHSEQ = 'enrichseq/output_files/taxid_abundances.csv'
    FASTVIROMEEXPLORER = 'FastViromeExplorer-final-sorted-abundance.tsv'
    BRACKEN = 'abundances.bracken'


def parse_directory_trees(truth_path: str, result_path: str, output_prefix: str) -> None:
    '''
    DESCRIPTION:
        Get directory structure metadata from file paths and write to CSV.
        Handles truth and result directory trees
    
    INPUT:
        1. Path to truth directory tree
        2. Path to results directory tree (i.e. the results of the tools that were tested) 
        3. Path to desired place to store CSV file

    OUTPUT:
        CSV file containing: <trial #>,<tool>,<test>,<subtest>,<path to final result file>
        e.g. 3,EnrichSeq,num_genomes,10_genomes_500000_reads,/projects/benchmarking/results/test_3/EnrichSeq/num_genomes/10_genomes_500000_reads/enrichseq/output_files/taxid_abundances.csv
    '''
    
    # create CSV file to append to
    output_csv_path = Path(settings.OUTPUT_DIRECTORY) / Path(output_prefix + '.csv')

    with open(output_csv_path, 'a') as csvfile:
        writer = csv.writer(csvfile) 
        writer.writerow(['Trial_Num', 'Experiment', 'Condition', 'Tool', 'File_Path'])

    # parse truth
    truth_rows = parse_truth_filepaths(truth_path)
    file = write_to_csv(truth_rows, output_csv_path)

    # parse results
    result_rows = parse_results_filepaths(result_path)
    file = write_to_csv(result_rows, output_csv_path)
    
    print(f'Metadata saved to:{file}')


def parse_truth_filepaths(dir_path: str) -> list:
    '''
    DESCRIPTION:
        Get metadata from truth directory paths and store in list
        
    INPUT:
        Path to all truth files (called "tests" in our case)
    
    OUTPUT:
        A list of lists. Each row contains the metadata for each "truth" entry
    '''

    rows = []

    for trial in os.listdir(dir_path):
        try:
            trial_num = int(''.join(filter(str.isdigit, trial)))
        except ValueError as e:
            trial_num = 0
            print(f"Directory does not have an integer, causing error: {e}", file=sys.stderr)
        finally:
            trial_path = Path(dir_path) / trial
            if os.path.isdir(trial_path):
                for test in os.listdir(trial_path):
                    test_path = Path(trial_path) / test
                    if os.path.isdir(test_path) and 'config' not in test:
                        for subtest in os.listdir(test_path):
                            subtest_path = Path(test_path) / subtest
                            if os.path.isfile(subtest_path):
                                if subtest.endswith('.fa'):
                                    info = [trial_num, test, subtest.strip('.fa'), 'truth', os.path.abspath(subtest_path)] 
                                    rows.append(info)
                 
    return rows



def parse_results_filepaths(results_path: str) -> list:
    '''
    DESCRIPTION:
        Get metadata from tools' results directory paths and store in list
        
    INPUT:
        Path to all tools' results files.
        Expected directory structure for results: 
            results/<trial num>/<tool name>/<experiment>/<condition>/
        Example expected structure:
            results/trial_2/EnrichSeq/num_genomes/200_genomes_500000_reads/...
        * up to 9 trials
    
    OUTPUT:
        A list of lists. Each row contains the metadata for each tool's results
    '''
    result_rows = []
    for trial in os.listdir(results_path):
        try:
            trial_num = int(''.join(filter(str.isdigit, trial)))
        except ValueError as e:
            trial_num = 0
            print(f"Directory does not have an integer, causing error: {e}", file=sys.stderr)
        finally:
            trial_path = Path(results_path) / trial
            if os.path.isdir(trial_path):
                for tool in os.listdir(trial_path):
                    tool_path = Path(trial_path) / tool
                    if os.path.isdir(tool_path):
                        #result_rows.append(__parser_selector__(tool_path, trial_num))
                        for test in os.listdir(tool_path):
                            test_path = Path(tool_path) / test
                            if os.path.isdir(test_path):
                                for condition in os.listdir(test_path):
                                    condition_path = Path(test_path) / condition
                                    if os.path.isdir(condition_path):
                                        # Now in tool's output directory
                                        result_info = [trial_num, test, condition]
                                        tool_info = _tool_selector(tool, condition_path)
                                        if tool_info: # if the info does not come back empty, append the row
                                            result_info.extend(tool_info)
                                            result_rows.append(result_info)
    return result_rows


def _tool_selector(tool_name: str, tool_result_path: Path) -> list:
    '''
    DESCRIPTION:
        Helper function to add the correct result file path of each tool (as they are all different)
        uses Enum
        
    INPUT:
        1. The tool name - to use as a selector for choosing the appropriate path
        2. Tool result file path - location of result files of the given tool
    
    OUTPUT:
        A list containing the metadata for the given tool's result for that experiment condition ("subtest")
    '''
    tool_name_uppercase = tool_name.upper()
    for tool_name_enum in ToolPathFinder:
        if tool_name_enum.name == tool_name_uppercase:
            full_result_path = Path(tool_result_path) / tool_name_enum.value
            if os.path.exists(full_result_path) and os.path.isfile(full_result_path):
                result_row = [tool_name_uppercase, full_result_path]

    return result_row



def write_to_csv(rows: list, output_file: str) -> str:
    '''
    DESCRIPTION:
        Writes rows to a CSV file
        
    INPUT:
        1. A list of rows to write to the CSV file.
        2. The file to write to
    
    OUTPUT:
        Returns the CSV file path (string) in which the metadata was stored
    '''
    with open(output_file, 'a') as csvfile: 
        writer = csv.writer(csvfile) 
        writer.writerows(rows)

    return output_file


def parseArgs(argv=None) -> argparse.Namespace:
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
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-t", "--truths_path", help="directory path containing all truth files", required=True)
    parser.add_argument("-r", "--results_path", help="directory path containing tool result files", required=True)
    parser.add_argument("-o", "--output_prefix", help="prefix for metadata csv file", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    settings = Dynaconf(settings_files="settings.toml")
    arguments = parseArgs(argv=sys.argv[1:])
    parse_directory_trees(arguments.truths_path, arguments.results_path, arguments.output_prefix)