import sys
from re import search


"""
v2021.02.11
Description:
This file parses the names.dmp NCBI file to find tax ids, then modifies phagesDB file to include them.
Extended from findnames.py --> previous version finds taxon id given a specific name

PREREQUISITE:
* Must have installed the phagesDB reference fasta file "Actinobacteriophages-All.fasta"
* Must have run kraken2-build --download-taxonomy --db <DBNAME>. names.dmp file will be in DBNAME/taxonomy/

USAGE:
python kraken_db_format.py <NCBI names.dmp file> <DB FASTA file to be formatted> <output destination>
e.g. python kraken_db_format.py krakenDB/taxonomy/names.dmp Actinobacteriophages-All.fasta phage_genomes/
"""

def parseIntoDictionary(dmp_file):
    """ parse file into data structure"""
    dict_tx_names = {}
    for line in dmp_file:
        line = line.strip("\n").replace('\t', '').split("|")
        taxid = line[0]
        species_name = line[1]
        line_type  = line[3]
        if "perseus" in line:
            print(line_type)
        if "scientific" in line_type and "phage" in species_name:
            dict_tx_names[species_name] = taxid

    return dict_tx_names

def findSpeciesId(tax_dict, name):
    """ finds the taxid """
    #return tax_dict[name]
    for key in tax_dict.keys():
        # below could be improved using regex
        #if key.lower().find(name.lower()) != -1:
        #if search(name.lower().replace(" ", ""), key.lower().replace(" ", "")):
        if search(name.lower(), key.lower()):
            return tax_dict[key]

    # if it doesn't find anything
    return 0


def addTaxIdToFile(tax_dict, fasta_file, output_path):
    """ adds taxid to the phagesDB fasta file for kraken2 format """
    fasta_file_contents = open(fasta_file).readlines()
    #names_file_contents = open(names_file).readlines() #list

    # new DB file to add taxon id to
    filename = "phagesDB.fasta"
    taxonfile = open(output_path + filename, "w")

    for line in fasta_file_contents:
        line = line.strip("\n") # will this affect the final file format when using in kraken?
        if '>' in line:
            species_name = line.split(' ')[2]
            taxid = findSpeciesId(tax_dict, species_name)
            phage_name_line = line.split()

            #print(f'Name: {species_name}, taxid: {taxid}') #works
            #newline = line[0] + line[1] + "|kraken:taxid|" + str(taxid)

            newline = ' '.join(phage_name_line[:3:]) + "|kraken:taxid|" + str(taxid) + " " + ' '.join(phage_name_line[3::])
            #print(newline)
            # write this modified line to the new file
            taxonfile.write(newline + "\n")
        else:
            # write the rest of the nucleotide info into the new file (unchanged)
            taxonfile.write(line + "\n")
    #print(fasta_file_contents[100].split(' ')[2]) # prints the phage name

    taxonfile.close()
    return filename


def main():
    """ controls the script """

    ### SCRIPT INPUT
    dmp_file_path = sys.argv[1]
    db_file_path = sys.argv[2]
    output_path = sys.argv[3]

    ### SCRIPT
    file_name = open(dmp_file_path).readlines()
    tax_map = parseIntoDictionary(file_name)


    # FIND NAMES
    #name_file = open(line_sep_path).readlines()

    # for species_name in name_file:
    #     species_name = str(species_name).strip("\n")
    #     taxid = findSpeciesId(tax_map, species_name)
        #print(species_name, taxid)

    # MODIFY DB FILE (ADD TAXON ID)
    print('Creating kraken-compatible database file...')
    f = addTaxIdToFile(tax_map, db_file_path, output_path)
    print(f'File created: {f}')


if __name__ == "__main__":
    main()
