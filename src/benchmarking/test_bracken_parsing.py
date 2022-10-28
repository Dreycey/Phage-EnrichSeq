import sys
import argparse
from typing import Tuple, List, Dict, Set

def parse_bracken_file(file) -> set:
    delim = '\t'
    tax_id_column = 1
    abundance_column = 6
    bracken_taxid_set = set()
    try:
        bracken_file_opened = open(file, 'r')
    except OSError as e:
        print(f"Unable to open {file}: {e}", file=sys.stderr)
    else:
        bracken_file_lines = bracken_file_opened.readlines()[1:]
        for line in bracken_file_lines:
            tax_id = line.strip('\n').split(delim)[tax_id_column]
            #abundance_val = line.strip('\n').split(delim)[abundance_column]
            bracken_taxid_set.add(int(tax_id))
        bracken_file_opened.close()
    
    return bracken_taxid_set


def parse_kraken_file(file) -> set:
    delim = '\t'
    species_level_column = 3
    tax_id_column = 4
    kraken_taxid_set = set()
    try:
        kraken_file_opened = open(file, 'r')
    except OSError as e:
        print(f"Unable to open {file}: {e}", file=sys.stderr)
    else:
        kraken_file_lines = kraken_file_opened.readlines()
        for line in kraken_file_lines:
            species_level = line.strip('\n').split(delim)[species_level_column]
            if species_level == 'S':
                tax_id = line.strip('\n').split(delim)[tax_id_column]
                kraken_taxid_set.add(int(tax_id))
        kraken_file_opened.close()
    
    return kraken_taxid_set


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
    parser.add_argument("-k", "--kraken_file", help="kraken report file", required=True)
    parser.add_argument("-b", "--bracken_file", help="bracken output file", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":

    arguments = parseArgs(argv=sys.argv[1:])
    kraken_file = arguments.kraken_file
    bracken_file = arguments.bracken_file

    kraken_set = parse_kraken_file(kraken_file)
    bracken_set = parse_bracken_file(bracken_file)

    if kraken_set == bracken_set:
        print('Sets are equal')
    else:
        taxids_to_check = kraken_set ^ bracken_set
        print(f'Unequal sets. Taxids to check: {taxids_to_check} ')


    
    

