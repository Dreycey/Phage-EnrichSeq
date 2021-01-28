import sys
from re import search


"""
Description:
This file parses the names.dmp NCBI file to find tax ids.

USAGE:
python findnames.py names.dmp name_by_line.txt 
"""

def parseIntoDictionary(file_in):
    """ parse file into data structure"""
    dict_tx_names = {}
    for line in file_in:
        line = line.strip("\n").replace('\t', '').split("|")
        taxid = line[0]
        species_name = line[1]
        line_type  = line[3]
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

def main():
    """ controls the script """

    ### SCRIPT INPUT
    dmp_file_path = sys.argv[1]
    line_sep_path = sys.argv[2]
    
    ### SCRIPT
    file_name = open(dmp_file_path).readlines()
    tax_map = parseIntoDictionary(file_name)

    # FIND NAMES
    name_file = open(line_sep_path).readlines()
    for species_name in name_file:
        species_name = str(species_name).strip("\n")
        taxid = findSpeciesId(tax_map, species_name)
        print(species_name, taxid)

if __name__ == "__main__":
    main()
