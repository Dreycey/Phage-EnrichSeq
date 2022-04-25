"""
Description:

Usage:
python db_compare.py ref_genomes/ krakenDB/seqid2taxid.map 
Example:
python db_compare.py ref_genomes/ krakenDB/seqid2taxid.map 

"""
from typing import Dict
import glob
import sys



def fasta2dict(genome_path) -> Dict:
    """
    Description: scrapes the reference db for taaxids
    Input: <str> a path of directory genomes
    Output: <Dict> name2seq
    """
    name2seq : Dict = {}
    seq_name : str = ""
    with open(genome_path) as fasta_dir:
        curr_line = fasta_dir.readline()
        while (curr_line):
            if (">" == curr_line[0]):
                seq_name = curr_line[1:].strip("\n")
                name2seq[seq_name] = ""
            else:
                name2seq[seq_name] += curr_line.strip("\n")
            # read next line below; close
            curr_line = fasta_dir.readline()
    return name2seq

def parse_seqid2taxid(seqid2taxid_path):
    """
    Description: scrapes the seqid2taxid file for taxids
    Input: <str> seqid2taxid
    Output: <Dict> name2seq
    """
    name2seq = {}
    with open(seqid2taxid_path) as seqid2tax_file:
        seqid2tax = seqid2tax_file.readline()
        while (seqid2tax):
            name, taxid = seqid2tax.split("\t")
            name2seq[taxid.strip("\n")] = name
            # read next line below; close
            seqid2tax = seqid2tax_file.readline()
    return name2seq

def fastdir2dict(seqid2taxid_path):
    """
    Description: scrapes a fasta folder for taxidss
    Input: <str> directtory path
    Output: <Dict> Dict
    """
    unitColor = '\033[5;36m\033[5;47m'
    endColor = '\033[0;0m\033[0;0m'
    name2seq: Dict = {}
    fasta_dir = glob.iglob(seqid2taxid_path + "*")
    fasta_dir_p = True
    i = 0
    count = 6875
    for fasta_dir_p in fasta_dir:
        name2seq.update(fasta2dict(fasta_dir_p))
        incre = int(50.0 / count * i)
        if (i != count - 1):
            sys.stdout.write('\r' + '|%s%s%s%s| %d%%' % (unitColor, '\033[7m' + ' '*incre + ' \033[27m', endColor, ' '*(50-incre), 2*incre))
        else:
            sys.stdout.write('\r' + '|%s%s%s| %d%%' % (unitColor, '\033[7m' + ' '*20 + 'COMPLETE!' + ' '*21 + ' \033[27m', endColor, 100) + "\n")
        sys.stdout.flush()
        i+=1
    return name2seq

def get_diffs(kraken_dict, ref_genome):
    """
    get difference summary
    """
    kraken_dict = set(kraken_dict.keys())
    ref_set = set()
    for name in ref_genome.keys():
        try: 
            taxid = name.split("kraken:taxid|")[1].split(" ")[0]
            print(taxid)
            ref_set.add(taxid)
        except:
            continue
    print(f"Out of {len(ref_set)}, {len(ref_set.intersection(kraken_dict))} are shared")
    print(f"")

if __name__ == "__main__":
    ref_seqid2taxid_argpath = sys.argv[1]
    kraken_seqid2taxid_argpath = sys.argv[2]
    ref_names2taxids = fastdir2dict(ref_seqid2taxid_argpath)
    kraken_names2taxids = parse_seqid2taxid(kraken_seqid2taxid_argpath)

    # comparison
    get_diffs(kraken_names2taxids, ref_names2taxids)