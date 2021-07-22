"""
This script extracts phage names from the kraken output file (intended file is kraken post-
de novo assembly).

USAGE:
    python parse_kraken.py <kraken_outfile> <outputfile_taxids> <outputfile_names>
"""
import sys
import csv
import re


"""
Extracts phage names from kraken report (only at species level)

Input:
    kraken report file

Output:
    returns list of phages
"""
def parseKrakenFile(kraken_file):
    kraken_reports = open(kraken_file).readlines()
    phages = {}

    # TODO: Add 'S1' in if statement and don't include duplicates
    for line in kraken_reports:
        genome_level = line.strip("\n").split('\t')[3]
        full_name = line.strip("\n").split('\t')[5]
        # filter by species
        if genome_level == 'S' or genome_level == 'S1':
            taxid = line.strip("\n").split()[4]
            # check that the phage doesn't already exist
            if not checkPhageExists(full_name, phages):
                phages[taxid] = full_name
                print(f'{taxid} {full_name}')

    return phages


"""
Check if phage name already exists in the dictionary

INPUT:
    full name
    phage dictionary

OUTPUT:
    returns true if it alredy exists
            false if not
"""
def checkPhageExists(phage_name, phage_dict):
    shortened_name = phage_name.split()[2]
    for key in phage_dict:
        if shortened_name in phage_dict[key]:
            return True
    return False


"""
Saves phage taxon ids in list to specified output file

Input:
    output file to open and write to
    list to store phage taxon ids

Output:
    none (writes to file)
"""

def saveTaxidToFile(outfile, kraken_phages):
    # save only phage taxids specified output file
    output_file = open(outfile, "w")
    for taxid in kraken_phages:
        output_file.write(taxid + "\n")
    output_file.close()

"""
Saves phage names in list to specified output file

Input:
    output file to open and write to
    list to store phage names

Output:
    none (writes to file)
"""

def saveNameToFile(outfile, kraken_phages):
    # save only phage taxids specified output file
    output_file = open(outfile, "w")
    for taxid in kraken_phages:
        output_file.write(kraken_phages[taxid].strip() + "\n")
    output_file.close()


def main():
    """ controls the script """
    kraken_phages = {}

    ### SCRIPT INPUT
    kraken_file = sys.argv[1]
    outfile_taxid = sys.argv[2]
    outfile_names = sys.argv[3]

    # TEST
    kraken_phages = parseKrakenFile(kraken_file)
    saveTaxidToFile(outfile_taxid, kraken_phages)
    saveNameToFile(outfile_names, kraken_phages)



if __name__ == "__main__":
    main()
