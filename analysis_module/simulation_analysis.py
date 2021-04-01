# This script extracts info from: 1) simulation config files, 2) bracken files
# compares results, and saves the information in a CSV file
import sys
import csv
import re


data_dict = {}

def parsePhageFile(phage_file):
    phages = open(phage_file).readlines()

    for line in phages:
        phage_name = line.strip("\n")
        data_dict[phage_name] = []


def parseConfigFile(config_file):
    config_contents = open(config_file).readlines()[1:]

    # search for phage name in the config file
    for line in config_contents:
        line = line.strip("\n").split()
        phage_path = line[0] # path of phage genome is the first string
        sim_abundance = line[1] # simulated abundance is the second string
        for key in data_dict:
            # if the name matches that specified in the txt file, record its simulated abundance
            if re.search(key, phage_path, re.IGNORECASE):
                data_dict[key].append(sim_abundance)


def parseBrackenFile(bracken_file):
    bracken_contents = open(bracken_file).readlines()[1:]
    for line in bracken_contents:
        line = line.strip("\n").split()

        species_name = line[2].lower()
        bracken_abundance = line[8]

        for key in data_dict:
            if key == species_name:
                data_dict[key].append(bracken_abundance)


def calculateError():
    for key in data_dict:
        sim_abundance = float(data_dict[key][0])*100
        reported_abundance = float(data_dict[key][1])*100
        error = round(abs(sim_abundance - reported_abundance)/sim_abundance, 5)
        data_dict[key].append(error)


def saveInfoToFile(outfile):
    # save names, taxid and abundances of dict contents to file
    if '.csv' in outfile:
        with open(outfile, mode='w') as csv_file:
            fieldnames = ['phage_name', 'config_abundance', 'bracken_abundance', 'error']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
            writer.writeheader()

            for key in data_dict:
                writer.writerow({'phage_name': key, 'config_abundance': data_dict[key][0],
                                'bracken_abundance': data_dict[key][1], 'error': data_dict[key][2]})
    else:
        print("Not a CSV file")

def main():
    """ controls the script """

    ### SCRIPT INPUT
    phages = sys.argv[1]
    config_file = sys.argv[2]
    bracken_file = sys.argv[3]
    outfile = sys.argv[4]

    # TEST
    parsePhageFile(phages)
    parseConfigFile(config_file)
    parseBrackenFile(bracken_file)
    calculateError()
    saveInfoToFile(outfile)

    # for key in data_dict:
    #     print(f'phage: {key} | config abundance: {data_dict[key][0]} | '
    #           f'bracken abundance: {data_dict[key][1]} | error: {data_dict[key][2]}')


if __name__ == "__main__":
    main()



