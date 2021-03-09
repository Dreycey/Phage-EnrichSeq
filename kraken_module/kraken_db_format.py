import sys
from re import search
from tqdm import tqdm

"""
v2021.02.11
Description:
This file parses the names.dmp NCBI file to find tax ids, then modifies phagesDB file to include them.
Extended from findnames.py --> previous version finds taxon id given a specific name

PREREQUISITE:
* Must have installed the phagesDB reference fasta file "Actinobacteriophages-All.fasta"
* Must have run kraken2-build --download-taxonomy --db <DBNAME>. names.dmp file will be in DBNAME/taxonomy/

USAGE:
python kraken_db_format.py <NCBI names.dmp file> <DB FASTA file to be formatted> <Kraken outfile>
e.g. python kraken_db_format.py names.dmp Actinobacteriophages-All.fasta phageDB_out.fa
"""

def parseIntoDictionary(dmp_file):
    """ parse file into data structure"""
    dict_tx_names = {}
    for line in tqdm(dmp_file):
        line = line.strip("\n").replace('\t', '').split("|")
        taxid = line[0]
        species_name = line[1]
        line_type  = line[3]

        if "scientific" in line_type and "phage" in species_name: 
            dict_tx_names[species_name] = taxid   
        elif "scientific" in line_type and "virus" in species_name: 
            dict_tx_names[species_name] = taxid

    return dict_tx_names

def findSpeciesId(tax_dict, name):
    """ finds the taxid """
    for key in tax_dict.keys():
        if search(name.lower(), key.lower()):
            return tax_dict[key]

    # if it doesn't find anything
    return 0


def addTaxIdToFile(tax_dict, fasta_file, outfile):
    """ adds taxid to the phagesDB fasta file for kraken2 format """
    fasta_file_contents = open(fasta_file).readlines()

    # new DB file to add taxon id to
    taxonfile = open(outfile, "w")

    nonfound_count = 0 #counting number of names not found in names.dmp
    for line in tqdm(fasta_file_contents):
        line = line.strip("\n") # will this affect the final file format when using in kraken?
        if '>' in line:
            species_name = line.split(' ')[2].strip(",")
            taxid = findSpeciesId(tax_dict, species_name)
            if taxid != 0:
                phage_name_line = line.split()
                newline = '_'.join(phage_name_line[:3]) + "|kraken:taxid|" + str(taxid) + " " + ' '.join(phage_name_line[3::])
                taxonfile.write(newline + "\n")
            else:
                nonfound_count += 1
                continue
        else:
            taxonfile.write(line + "\n")

    taxonfile.close()
    return outfile, nonfound_count


def main():
    """ controls the script """

    ### SCRIPT INPUT
    dmp_file_path = sys.argv[1]
    db_file_path = sys.argv[2]
    outfile = sys.argv[3]

    ### SCRIPT
    file_name = open(dmp_file_path).readlines()
    tax_map = parseIntoDictionary(file_name)

    # MODIFY DB FILE (ADD TAXON ID)
    print('Creating kraken-compatible database file...')
    f,outnon = addTaxIdToFile(tax_map, db_file_path, outfile)
    print(f'File created: {f}')
    print(f"{outnon} phages did not have tax mappings")

if __name__ == "__main__":
    main()
