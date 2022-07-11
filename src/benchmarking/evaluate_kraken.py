"""
This script calculates precision, recall and F1 scores for Kraken results.

USAGE:
    python evaluate_kraken.py <kraken output file> <simulation file?>

kraken output file:
    0.01	1	1	S	2927287	      Arthrobacter phage Tallboi
simulation file:
    >NC_024381.1-245872|GCF_000923035.1_ViralProj253340_genomic.fna
"""

import sys
import os
import csv
from pathlib import Path


# Changing to directory of script.
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
#os.chdir("../../")
working_dir = os.getcwd()

def parse_predicted(kraken_file):
    predicted_phages = {}
    prev = None
    count = 0
    with open(kraken_file) as file:
        for line in file:
            count += 1
            kraken_info = line.strip('\n').split('\t')
            if line and not line.isspace(): 
                genome_level = kraken_info[3]
                # Grab phage info at species level only
                if genome_level == 'S':
                    taxid = int(kraken_info[4])
                    full_name = kraken_info[5]
                    predicted_phages[taxid] = {full_name, genome_level}
                    prev = taxid
                elif genome_level == 'S1' and prev is not None:
                    predicted_phages.pop(prev)
                    taxid = kraken_info[4]
                    full_name = kraken_info[5]
                    predicted_phages[int(taxid)] = {full_name, genome_level}
                    prev = None
    print(f'{count} lines read')
    # # clean entries: if S is followed by S1, remove S row
    # for taxid_key in predicted_phages:
    #     # Get next key in Dictionary
    #     next_taxid = None
    #     temp = iter(predicted_phages)
    #     count_removals = 0
    #     for key in temp:
    #         if key == taxid_key:
    #             next_taxid = next(temp, None)
    #             if predicted_phages[taxid_key][1] == 'S' and predicted_phages[next_taxid][1] == 'S1':
    #                 print(f'Removing phage species {taxid_key} due to strain present {next_taxid}...')
    #                 predicted_phages.pop(taxid_key)
    #                 count_removals += 1
    # print(f'Total removals: {count_removals}')

    return predicted_phages


def parse_truth(sim_file):
    true_phages = {}
    files = set()
    with open(sim_file) as file:
        print(f'Opening simulation file: {sim_file}')
        for line in file:
            if line.startswith(">"):
                info = line.strip('\n').split('|')
                files.add(info[1])
    # try:
    #     file = open(sim_file)
    # except FileNotFoundError:
    #     print(f'{sim_file} not found.')
    
    database_path = Path(working_dir) / Path('database/ref_genomes/')
    for phage_file in files:
        try:
            file = open(database_path / Path(phage_file))
        except FileNotFoundError:
            print(f'{phage_file} not found in {database_path}')
        else:
            line = file.readline()
            if line.startswith(">"):
                file_info = line.strip('\n').split('|kraken:taxid|')
                true_phages[int(file_info[1].strip(' '))] = file_info[0]
                file.close()
                continue

        # with open(database_path / Path(phage_file)) as file:
        #     line = file.readline()
        #     if line.startswith(">"):
        #         file_info = line.strip('\n').split('|kraken:taxid|')
        #         true_phages[int(file_info[1].strip(' '))] = file_info[0]

    return true_phages


def get_true_positives(predicted_phages, true_phages):
    """
    True positives: phage in kraken report AND in simulated file
    """
    true_positives = []
    for pred_taxid in predicted_phages.keys():
        if pred_taxid in true_phages.keys():
            true_positives.append(pred_taxid)
    print(f'true positives: {true_positives}')
    return len(true_positives)


def get_false_positives(predicted_phages, true_phages):
    false_positives = []
    for pred_taxid in predicted_phages.keys():
        if (pred_taxid not in true_phages.keys()) and (pred_taxid != "UNK"):
            false_positives.append(pred_taxid)
    return len(false_positives)


def get_false_negatives(predicted_phages, true_phages):
    false_negatives = []
    for pred_taxid in true_phages.keys():
        if pred_taxid not in predicted_phages.keys():
            false_negatives.append(pred_taxid)
    return len(false_negatives)


def calc_precision(true_positives, false_positives):
    """
    precision = true positives / (true positives + false positives)
    """
    precision = true_positives / (true_positives + false_positives)
    return round(precision,2)


def calc_recall(true_positives, false_negatives):
    """
    recall = true positives / (true positives + false negatives)
    """
    recall = true_positives / (true_positives + false_negatives)
    return round(recall,2)


def calc_f1_score(precision, recall):
    """
    f1 = 2 * precision * recall / (precision + recall)
    """
    f1 = 2 * precision * recall / (precision + recall)
    return round(f1,2)


def main():
    """ controls the script """

    ### SCRIPT INPUT
    kraken_file = Path(sys.argv[1])
    sim_file = Path(sys.argv[2])

    print(working_dir)
    predicted_phages = parse_predicted(kraken_file)
    true_phages = parse_truth(sim_file)
    true_pos = get_true_positives(predicted_phages, true_phages)
    print(f'True positives are: {true_pos}')
    false_pos = get_false_positives(predicted_phages, true_phages)
    print(f'False positives are: {false_pos}')
    false_neg = get_false_negatives(predicted_phages, true_phages)
    print(f'False negatives are: {false_neg}')

    precision = calc_precision(true_pos, false_pos)
    recall = calc_recall(true_pos, false_neg)
    f1_score = calc_f1_score(precision, recall)

    print(f'Precision = {precision}')
    print(f'Recall = {recall}')
    print(f'F1 score = {f1_score}')

    # # error handling
    # if (len(true_abundances) == 0): 
    #     print(f"TRUTH EMPTY")
    #     return results_structure


if __name__ == "__main__":
    main()