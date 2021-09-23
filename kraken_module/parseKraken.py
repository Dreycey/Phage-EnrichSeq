"""
This script extracts phage names from the kraken output file (intended file is kraken post-
de novo assembly).

USAGE:
    python parse_kraken.py <kraken_outfile> <outputfile_taxids> <outputfile_names>
"""
import sys
import re
from pathlib import Path


def parseKrakenFile(kraken_file) -> dict:
    """
    DESCRIPTION:
        Extracts phage names from kraken report (only at species level)

    INPUT:
        kraken report file

    OUTPUT:
        returns list of phages
    """
    phages = {}
    try:
        kraken_reports = open(kraken_file).readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f'{kraken_file} does not exist.')
    else:
        # TODO: don't include duplicates
        for line in kraken_reports:
            line_array = line.strip("\n").split('\t')
            if line and not line.isspace(): 
                genome_level = line_array[3]
                full_name = line_array[5]
                # filter by species
                if genome_level == 'S' or genome_level == 'S1':
                    taxid = line.strip("\n").split()[4]
                    # check that the phage doesn't already exist
                    if not checkPhageExists(full_name, phages):
                        phages[taxid] = full_name.strip()

    return phages


def checkPhageExists(phage_name, phage_dict) -> bool:
    """
    DESCRIPTION:
        Check if phage name already exists in the dictionary

    INPUT:
        full name
        phage dictionary

    OUTPUT:
        returns true if it alredy exists
                false if not
    """
    if len(phage_name.split()) > 2:
        shortened_name = phage_name.split()[2]
    elif len(phage_name.split()) == 1:
        shortened_name = phage_name
    for key in phage_dict:
        if re.search(shortened_name, phage_dict[key], re.IGNORECASE):
            return True
    return False


def saveTaxidToFile(outfile, kraken_phages) -> Path:
    """
    DESCRIPTION:
        Saves phage taxon ids in list to specified output file

    INPUT:
        output file to open and write to
        list to store phage taxon ids

    OUTPUT:
        none (writes to file)
    """ 

    if not kraken_phages:
        return None
    else:
        output_file = open(outfile, "w")
        for taxid in kraken_phages:
            output_file.write(taxid + "\n")
        output_file.close()
        return Path(outfile)


def saveNameToFile(outfile, kraken_phages) -> Path:
    """
    DESCRIPTION:
        Saves phage names in list to specified output file

    INPUT:
        output file to open and write to
        list to store phage names

    OUTPUT:
        file path or None if no phages found
    """
    if not kraken_phages:
        return None
    else:
        output_file = open(outfile, "w")
        for taxid in kraken_phages:
            output_file.write(kraken_phages[taxid].strip() + "\n")
        output_file.close()
        return Path(outfile)


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
