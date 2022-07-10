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
from pathlib import Path


def parse_predicted(kraken_file):
    predicted = {}
    with open(kraken_file) as file:
        line = file.readline().rstrip()
        kraken_info = line.strip('\t')
        if line and not line.isspace(): 
                genome_level = kraken_info[3]
                # Grab phage info at species level only
                if genome_level == 'S' or genome_level == 'S1':
                    taxid = kraken_info[4]
                    full_name = kraken_info[5]
                    predicted[taxid] = {full_name, genome_level}
                    # TODO: if S is followed by S1, remove S row?
    
    return predicted


def parse_truth(sim_file):
    truth = {}
    with open(sim_file) as file:
        line = file.readline().rstrip()


def get_true_positives(kraken_file, sim_file):
    """
    True positives: phage in kraken report AND in simulated file
    """
    true_positives = []
    
    return len(true_positives)


def get_false_positives(kraken_file, sim_file):
    pass


def get_true_negatives(kraken_file, sim_file):
    pass


def get_false_negatives(kraken_file, sim_file):
    pass

def calc_precision(true_positives, false_positives):
    """
    precision = true positives / (true positives + false positives)
    """
    precision = true_positives / (true_positives + false_positives)
    return precision


def calc_recall(true_positives, false_negatives):
    """
    recall = true positives / (true positives + false negatives)
    """
    recall = true_positives / (true_positives + false_negatives)
    return recall


def calc_f1_score(precision, recall):
    """
    f1 = 2 * precision * recall / (precision + recall)
    """
    f1 = 2 * precision * recall / (precision + recall)
    return f1


def main():
    """ controls the script """

    ### SCRIPT INPUT
    kraken_file = Path(sys.argv[1])
    sim_file = Path(sys.argv[2])


    precision = calc_precision(get_true_positives(kraken_file, sim_file), get_false_positives(kraken_file, sim_file))
    recall = calc_recall(get_true_positives(kraken_file, sim_file), get_false_negatives(kraken_file, sim_file))
    f1_score = calc_f1_score(precision, recall)

    print(f'Precision = {precision}')
    print(f'Recall = {recall}')
    print(f'F1 score = {f1_score}')

    # error handling
    if (len(true_abundances) == 0): 
        print(f"TRUTH EMPTY")
        return results_structure
    
    # classification results
    ## true positives
    true_positives = []
    for pred_taxid in predicted_abundances.keys():
        if pred_taxid in true_abundances.keys():
            true_positives.append(pred_taxid)
    ## false positives
    false_positives = []
    for pred_taxid in predicted_abundances.keys():
        if (pred_taxid not in true_abundances.keys()) and (pred_taxid != "UNK"):
            false_positives.append(pred_taxid)
    
    ## false negatives
    false_negatives = []
    for pred_taxid in true_abundances.keys():
        if pred_taxid not in predicted_abundances.keys():
            false_negatives.append(pred_taxid)

if __name__ == "__main__":
    main()