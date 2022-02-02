"""
DESCRIPTION:
    This wrapper script executes the primary function
    of EnrichSeq by sending subcommands to the 
    build script and nextflow.

positional arguments:
  {db_build,enrichseq}  use enrichseq or build db.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose
"""
import argparse
import sys
import os
from enum import Enum
import src.py_modules.utils as utils




CURR_PATH = os.path.realpath(__file__)
CURR_PATH = os.path.dirname(CURR_PATH)

class SubparserNames(Enum):
    db_build = "db_build"
    enirchseq = "enrichseq"

def parseArgs(argv=None) -> argparse.Namespace:
    """
    This method takes in the arguments from the command and performs
    parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    """
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(help='use enrichseq or build db.', dest='sub_parser')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true")

    # if db_build
    db_build_parser = subparsers.add_parser(SubparserNames.db_build.value)
    # if enrichseq
    enrich_parser = subparsers.add_parser(SubparserNames.enirchseq.value)
    enrich_parser.add_argument("-1", "--input_1", help="the input fasta file (single or paired end 1)", required=True)
    enrich_parser.add_argument("-2", "--input_2", help="the input fasta file (single or paired end 2)", required=False)
    enrich_parser.add_argument("-o", "--output", help="the path to the output directory", required=True)
    enrich_parser.add_argument("-kdb", "--kracken_db", help="the path to the kraken database", required=False)
    enrich_parser.add_argument("-gdb", "--genome_db", help="the path to the phage genome directory", required=False)
    enrich_parser.add_argument("-t", "--threads", help="number of threads to use [Default 4]", required=False)
    return parser.parse_args(argv)

def run_dbbuild():
    """
    This method calls a sub process run for the DB build
    """
    print(" \n Building the database \n")
    CMD_list = ["bash", f"{CURR_PATH}/src/build_database/enrichseqDB_install.sh"]
    out = utils.subproc_call(CMD_list, shell=False)
    return out 

def run_enrichseq(primary_args):
    """
    This method calls a sub process to run EnrichSeq
    """
    print(" \n Running Enrichseq \n")
    CMD_list = ["nextflow", f"{CURR_PATH}/src/nextflow/enrichseq.nf"]
    if primary_args.input_2:
        CMD_list += ["--read", "paired"]
        CMD_list += ["--1", primary_args.input_1]
        CMD_list += ["--2", primary_args.input_2]
    else:
        CMD_list += ["--read", "single"]
        CMD_list += ["--fasta", primary_args.input_1]
    CMD_list += ["--toolpath", f"{CURR_PATH}/src/modules/"]
    if primary_args.kracken_db:
        CMD_list += ["--dbdir", primary_args.kracken_db]
    else:
        CMD_list += ["--dbdir", f"{CURR_PATH}/database/krakenDB/"]
    if primary_args.threads:
        CMD_list += ["--threads", primary_args.threads]
    else:
        CMD_list += ["--threads", "4"]
    if primary_args.genome_db:
        CMD_list += ["--genomedir", primary_args.genome_db]
    else:
        CMD_list += ["--genomedir", f"{CURR_PATH}/database/ref_genomes/"]
    CMD_list += ["--workdir", primary_args.output]
    out = utils.subproc_call(CMD_list)
    return out 

def main():
    primary_args = parseArgs(sys.argv[1:])
    if primary_args.sub_parser == SubparserNames.db_build.value:
        run_dbbuild()
    elif primary_args.sub_parser == SubparserNames.enirchseq.value:
        run_enrichseq(primary_args)
    else:
        print(__doc__)

if __name__ == "__main__":
    main()
