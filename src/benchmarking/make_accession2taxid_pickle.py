import pickle


file_with_names = "./abundance.tsv"
accession2taxid_file = "../../database/krakenDB/taxonomy/nucl_gb.accession2taxid"


def get_needed_ncbi_ids():
    """
    gets NCBI ids
    """
    file_with_names_opened = open(file_with_names)
    
    # grab all needed NCBI ids
    ncbi_ids = []
    line = file_with_names_opened.readline()
    line_counter = 0
    while(line):
        if (line_counter > 0): # skip header
            ncbi_id = line.split("\t")[0]
            ncbi_ids.append(ncbi_id)
        line = file_with_names_opened.readline()
        line_counter += 1
        
    file_with_names_opened.close()
    
    return ncbi_ids

ncbi_ids = get_needed_ncbi_ids()


dictionary = {}

# grab taxid mappings
acc2tax_open = open(accession2taxid_file)
line = acc2tax_open.readline()
counter = 0
while(line):
    line = acc2tax_open.readline()
    line = line.split("\t")
    accession = line[1]
    taxid = line[2]
    if accession in ncbi_ids:
        dictionary[accession] = taxid
    if (len(dictionary.keys()) == len(ncbi_ids)):
        break
    counter += 1

    if (counter % 100000 == 0):
        print(f"{counter} lines have been parsed \n")
        print(dictionary)
acc2tax_open.close()

with open('accession2taxid.pickle', 'wb') as handle:
    pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
