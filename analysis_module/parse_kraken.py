"""
This script extracts phage names from the kraken output file (intended file is kraken post-
de novo assembly).
"""
import sys
import csv
import re


data_dict = {}

"""
Extracts phage names from kraken report (only at species level)

Input:
    kraken report file
    
Output: 
    returns list of phages
"""
def parseKrakenFile(kraken_file):
    kraken_reports = open(kraken_file).readlines()
    phages = []

    for line in kraken_reports:
        genome_level = line.strip("\n").split('\t')[3]
        full_name = line.strip("\n").split('\t')[5]
        if genome_level == 'S': ## species have 3 names (e.g. Mycobacterium phage D29)
            phages.append(full_name)
            print(full_name)

    return phages


"""
Cross-references phage names from kraken report with multi-fasta reference genomes

Input:
    (full) phage name from kraken report
    multifasta file path

Output: 
    ?
"""

# def __validatePhageNames(phage_list, multifasta_file):
#


"""
Saves phage names in list to specified output file

Input:
    output file to open and write to
    list to store phage names

Output: 
    none (writes to file)
"""
def saveInfoToFile(outfile, kraken_phages):
    # save only phage names specified output file
    output_file = open(outfile, "w")
    for phage in kraken_phages:
        output_file.write(phage + "\n")
    output_file.close()

def main():
    """ controls the script """
    kraken_phages = []


    ### SCRIPT INPUT
    kraken_file = sys.argv[1]
    outfile = sys.argv[2]


    # TEST
    kraken_phages = parseKrakenFile(kraken_file)
    saveInfoToFile(outfile, kraken_phages)
    # parseConfigFile(config_file)
    # parseBrackenFile(bracken_file)
    # calculateError()
    # saveInfoToFile(outfile)

    # for key in data_dict:
    #     print(f'phage: {key} | config abundance: {data_dict[key][0]} | '
    #           f'bracken abundance: {data_dict[key][1]} | error: {data_dict[key][2]}')


if __name__ == "__main__":
    main()



