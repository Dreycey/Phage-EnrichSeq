#! usr/bin/python3
"""
DESCRIPTION:
    Python wrapper for dnadiff

USAGE:
    python dnadiff_wrapper.py -q <query fasta file> -r <reference fasta file>
"""

import os
import sys
import subprocess
import argparse



def call_dnadiff(query: str, reference: str, outfile: str):
    #cmd = 'dnadiff'
    #os.system(cmd)
    subprocess.run(['dnadiff', query, reference])


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
    parser.add_argument("-q", "--query", help="query genome file", required=True)
    parser.add_argument("-r","--reference", help="reference genome file", required=True)
    parser.add_argument("-o", "--output_prefix", help="output prefix for report directory", required=False)
    return parser.parse_args(argv)

def main():
    arguments = parseArgs(argv=sys.argv[1:])
    call_dnadiff(arguments.query, arguments.reference, arguments.output_prefix)


if __name__ == "__main__":
    main()