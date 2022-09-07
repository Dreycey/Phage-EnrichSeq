import os
from pathlib import Path
from typing import List, Dict

import pprint


# Find and return given file name in a directory structure


# Parse FastViromeExplorer
def parse_fastvirome(result_file_name: str):
    fve_result_dict = {}
    fve_path = Path('/projects/laal5512/results/test_1/FastViromeExplorer/num_genomes/10_genomes_500000_reads/')
    with open(fve_path / Path(result_file_name)) as result_file:
        result_file_lines = result_file.readlines()[1:]
        for line in result_file_lines:
            fve_result_dict[line[0]] = line[3]

    return fve_result_dict




if __name__ == "__main__":
    fve_dict = parse_fastvirome("FastViromeExplorer-final-sorted-abundance.tsv")
    pprint(fve_dict)
    print(f"Sum of estimated abundance = {sum(fve_dict)}")


