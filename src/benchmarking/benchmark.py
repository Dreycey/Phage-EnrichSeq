from functools import total_ordering
import sys
import csv
import json
from importlib.metadata import metadata

#from unittest import result
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Dict


NCBI_TO_TAXID_FILEPATH = '/Users/latifa/GitHub/benchmarking-enrichseq/out_dictionary.txt'
PHAGE_REF_DB_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tools/Phage-EnrichSeq/database/ref_genomes/'
KRAKEN_DB_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tools/Phage-EnrichSeq/database/KrakenDB/'
BENCHMARKING_OUTPUT_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/'

class SimulatedTruth:

    # SimulatedTruth object should only be created if it actually exists in the metadata CSV
    def __init__(self, trial_num, experiment_name, experiment_condition, sequencing_file) -> None:
        self.trial_num = trial_num
        self.experiment_name = experiment_name
        self.experiment_condition = experiment_condition
        self.sequencing_file = sequencing_file    
        self.num_genomes, self.num_reads = self.extract_genomes_and_reads()
        self.abundances_dict = {}
        
        # try:
        #     json_file = open(NCBI_TO_TAXID_FILEPATH)
        # except OSError as e:
        #     print(f"Unable to open {NCBI_TO_TAXID_FILEPATH}: {e}", file=sys.stderr)
        # else:
        #     self.ncbi2taxid = json.load(json_file)
        #     json_file.close()
        #self.abundances_dict = self.parse_simulated_fasta() # taxid -> abundance 
            

    
    def extract_genomes_and_reads(self) -> Tuple[int, int]:
        '''
        DESCRIPTION:

        INPUT:

        OUTPUT:
        '''
        num_genomes, num_reads = 0, 0
        if 'genomes' in self.experiment_condition and 'reads' in self.experiment_condition:
            values = self.experiment_condition.strip('genomes').strip('reads').split('_')
            try:
                num_genomes, num_reads = values[0], values[1]
            except IndexError as ex:
                print(f'Index out of bounds. Possibly caused by incorrect format.')
        else:
            print('Format of experiment condition title is not supported. Example: "10_genomes_1000000_reads".')      

        return num_genomes, num_reads
    

    def parse_simulated_fasta(self, ncbi_to_taxid_dict: Dict) -> Dict[str, float]:
        '''
        DESCRIPTION:
            Parses 'truth' sequencing file and stores tax IDs and corresponding abundances
            in a dictionary
        INPUT:
            -
        OUTPUT:
            Dictionary where <key : value> pair is taxid : abundance
        '''

        total_counts = 0
        not_found = set()
        prev_len = 0
        with open(self.sequencing_file, "r") as fasta_opened:
            line = fasta_opened.readline()
            while (line):
                if (line[0] == ">"):
                    ncbi_id = line[1:].split("|")[0].split("-")[0]
                    # convert NCBI to taxonomic ID
                    if ncbi_id in ncbi_to_taxid_dict: # TODO: this part should work
                        tax_id = ncbi_to_taxid_dict[ncbi_id]
                    else: # TODO: this part is not tested
                        print("NCBI ID not in ncbi2taxid dictionary")
                        file_name = line[1:].split("|")[-1].strip("\n")
                        path = Path(PHAGE_REF_DB_PATH + file_name)
                        try:
                            tax_id = open(path).readline().split("kraken:taxid|")[1].strip("\n").replace(" ", "")
                        except:
                            not_found.add(path)
                            new_len = len(not_found)
                            if new_len != prev_len:
                                print(f"ERROR (parsing):")
                                print(f"\tTool: simulated files")
                                print(f"\t\tnumber of true NOT found (fine if 1): {new_len}")
                                print(f"\t\tNCBI ID: {ncbi_id}")
                                prev_len = new_len
                    if tax_id in self.abundances_dict:
                        self.abundances_dict[tax_id] += 1
                    else:
                        self.abundances_dict[tax_id] = 1
                    total_counts += 1
                line = fasta_opened.readline()
            
        #turn counts into abundances - normalize
        for tax_id in self.abundances_dict.keys():
            self.abundances_dict[tax_id] /= total_counts 

        return self.abundances_dict


class Result:

    def __init__(self, trial_num, experiment_name, experiment_condition, tool_name, result_file, truth_obj, ncbi_to_taxid_dict = {}) -> None:
        self.trial_num = trial_num
        self.tool_name = tool_name
        self.experiment_name = experiment_name
        self.experiment_condition = experiment_condition
        self.num_genomes, self.num_reads = self.extract_genomes_and_reads()
        self.result_file = result_file
        self.truth_obj = truth_obj
        self.abundances_dict = {}
        self.true_positives = []
        self.false_positives = []
        self.false_negatives = []
        self.ncbi_to_taxid_dict = ncbi_to_taxid_dict
    

    def extract_genomes_and_reads(self) -> Tuple[int, int]:
        '''
        DESCRIPTION:

        INPUT:

        OUTPUT:
        '''
        num_genomes, num_reads = 0, 0
        if 'genomes' in self.experiment_condition and 'reads' in self.experiment_condition:
            values = self.experiment_condition.strip('genomes').strip('reads').split('_')
            try:
                num_genomes, num_reads = values[0], values[1]
            except IndexError as ex:
                print(f'Index out of bounds. Possibly caused by incorrect format.')
        else:
            print('Format of experiment condition title is not supported. Example: "10_genomes_1000000_reads".')      

        return num_genomes, num_reads


    def calc_classification_outcomes(self) -> Tuple[List, List, List]: 
        '''
        DESCRIPTION:
            Given the parsed data from the simulated truths and results files, this
            function finds the outcomes (true positives, fasle positives, false negatives)
            required for calculating F1 scores.
        
        INPUT:
            -
        
        OUTPUT:
            Assigns outcomes to the class's member variable lists (also returns the 3 lists)
        
        '''
        pass


    def calc_precision(self) -> float:
        if self.true_positives and self.false_positives and self.false_negatives:
            return len(self.true_positives) / (len(self.true_positives + len(self.false_positives)))
        else: # if outcomes lists are empty
            return 0.0


    def calc_recall(self) -> float:
        if self.true_positives and self.false_positives and self.false_negatives:
            return len(self.true_positives) / (len(self.true_positives + len(self.false_negatives)))
        else: # if outcomes lists are empty
            return 0.0


    def calc_f1score(self) -> float:
        precision = self.calc_precision()
        recall = self.calc_recall()
        return (precision * recall) / (precision + recall)


    def calc_l2distance(self) -> float:
        l2_distance = 0.0
        return l2_distance



class EnrichSeqResult(Result):
    delim = ','

    # (has everything parent class, Result, has)

    def parse_results(self):
        try:
            enrichseq_file_opened = open(self.result_file, 'r')
        except OSError as e:
            print(f"Unable to open {self.result_file}: {e}", file=sys.stderr)
        else:
            enrichseq_file_lines = enrichseq_file_opened.readlines()
            for line in enrichseq_file_lines:
                tax_id, abundance_val = line.strip('\n').split(self.delim)
                self.abundances_dict[tax_id] = abundance_val
            enrichseq_file_opened.close()

        return self.abundances_dict


class FVEResult(Result):
    delim = '\t'
    ncbi_id_column = 0
    count_column = 3
    #ncbi_to_taxid_dict = {}


    def parse_results(self) -> Dict:
        try:
            fve_file_opened = open(self.result_file, 'r')
        except OSError as e:
            print(f"Unable to open {self.result_file}: {e}", file=sys.stderr)
        else:
            total_abundance = 0
            fve_file_lines = fve_file_opened.readlines()[1:] # ignore header
            for line in fve_file_lines:
                ncbi_id = line.strip('\n').split(self.delim)[self.ncbi_id_column]
                counts_val = line.strip('\n').split(self.delim)[self.count_column]

                # Taxon ID lookup
                if ncbi_id in self.ncbi_to_taxid_dict:
                    tax_id = self.ncbi_to_taxid_dict[ncbi_id]
                elif ncbi_id.split(".")[0] in self.ncbi_to_taxid_dict: # TODO: check if this is needed @Dreycey
                    tax_id = self.ncbi_to_taxid_dict[ncbi_id.split(".")[0]]
                else:
                    tax_id = 'unk'
                    if float(counts_val) > 0.01:
                        print(f"ERROR (parsing):")
                        print(f"\tTool: FVE")
                        print(f"\t\tFVE: {ncbi_id} NOT CONVERTED")
                        print(f"\t\tabundances: {float(counts_val)}") 
                # Calculate 
                self.abundances_dict[tax_id] = float(counts_val)
                total_abundance += float(counts_val)
            fve_file_opened.close()

        #self._clean_abundances(total_abundance)
        return self.abundances_dict

    # TODO: May not need this, reading from a different FVE file
    def _clean_abundances(self, total_abundance):
        '''
        DESCRIPTION:
            Cleans the FVEResult object's abundances dictionary of any negligible 
            abundances (based on a threshold)
        '''
        ids_to_del = []
        for tax_id, abundance in self.abundances_dict.items():
            if total_abundance > 0:
                self.abundances_dict[tax_id] /= total_abundance # Normalize
            if self.abundances_dict[tax_id] < 0.001:
                ids_to_del.append(tax_id)

        # delete 0 value tax ids
        for val in ids_to_del:
            del self.abundances_dict[val]

        return self.abundances_dict

class BrackenResult(Result):
    delim = '\t'
    tax_id_column = 1
    abundance_column = 6

    def parse_results(self) -> Dict:
        try:
            bracken_file_opened = open(self.result_file, 'r')
        except OSError as e:
            print(f"Unable to open {self.result_file}: {e}", file=sys.stderr)
        else:
            bracken_file_lines = bracken_file_opened.readlines()[1:]
            for line in bracken_file_lines:
                tax_id = line.strip('\n').split(self.delim)[self.tax_id_column]
                abundance_val = line.strip('\n').split(self.delim)[self.abundance_column]
                self.abundances_dict[tax_id] = abundance_val
            bracken_file_opened.close()
        
        return self.abundances_dict


class Benchmarking:

    def __init__(self, metadata_csv = None) -> None:
        self.metadata_csv = metadata_csv
        self.result_objs = []

        # init metadata
        self.metadata_df = self.parse_metadata_csv()
    
    def parse_metadata_csv(self) -> pd.DataFrame:
        """
        Description:
            This parses the input csv into a pandas Dataframe.
        Errors:
            returns empty dataframe
        """
        try:
            metadata_df = pd.read_csv(self.metadata_csv)
        except OSError as e:
            print(f"Unable to open {self.metadata_csv}: {e}", file=sys.stderr)
            return pd.DataFrame()
        else:
            return metadata_df

    
    def assign_objs(self) -> float:
        """
        Description:
            This creates all of the result objects.
        Errors:
            returns empty array on error
        """
        # TODO: error handeling
        if self.metadata_df.empty:
            print("No metadata to parse results with!")
            return []
        unique_tests = self.metadata_df[["Trial_Num", "Experiment", "Condition"]].drop_duplicates()
        for index, unique_test_set in unique_tests.iterrows():
            trial_num_temp, experiment_temp, condition_num_temp = (unique_test_set.Trial_Num, 
                                                                   unique_test_set.Experiment, 
                                                                   unique_test_set.Condition)
            # filter matching rows (TODO: do in 1 step)
            filtered_data = self.metadata_df.loc[(self.metadata_df["Trial_Num"] == trial_num_temp) &
                                                 (self.metadata_df["Experiment"] == experiment_temp) &
                                                 (self.metadata_df["Condition"] == condition_num_temp)]
            # create a truth object
            truth_row_temp = filtered_data.loc[filtered_data["Tool"] == "truth"].drop(columns=['Tool']).values.tolist()
            # if not truth_row_temp:
            #     print("TRUTH IS EMPTY!")
            #     exit()
            #print(truth_row_temp)
            truth_obj = SimulatedTruth(*truth_row_temp[0]) 
            for index2, sub_row in filtered_data.iterrows():
                tool_name = sub_row["Tool"]
                if (tool_name != "truth"):
                    result_obj = self._return_tool_obj(tool_name, sub_row.values.tolist(), truth_obj, self.ncbi_to_taxid_mapping())
                    self.result_objs.append(result_obj)
            break

        #print(unique_tests)

    def _return_tool_obj(self, tool_name: str, result_array, truth_obj: SimulatedTruth, ncbi_to_taxid_mapping: Dict = {}):
        """
        creates result obj based on tool name.
        """
        name2obj = {
            "ENRICHSEQ" : EnrichSeqResult,
            "FASTVIROMEEXPLORER" : FVEResult,
            "BRACKEN" : BrackenResult
        }
        ToolClass = name2obj[tool_name.upper()]
        return ToolClass(*result_array, truth_obj, ncbi_to_taxid_mapping)


    def ncbi_to_taxid_mapping(self) -> Dict[str, str]: 
        '''
        DESCRIPTION:
            Only needs to be called once. Returns dictionary, which should be passed on to 
            objects that need an NCBI lookup (e.g. FVEResult)
        
        INPUT:
            -
        
        OUTPUT:
            Mapping dictionary where key is NCBI ID and value is taxon ID.
        '''
        ncbi_map_filepath = Path(BENCHMARKING_OUTPUT_PATH + 'ncbi_map.tab')
        
        ncbi_to_taxid_dict = {}
        try:
            lookup_file_opened = open(ncbi_map_filepath)
        except OSError as e:
            print(f'Unable to open {ncbi_map_filepath}: {e}', file=sys.stderr)
        else:
            line = lookup_file_opened.readline()
            while (line):
                line = line.split('\t')
                ncbi_id = line[1]
                tax_id = line[2]
                ncbi_to_taxid_dict[ncbi_id] = tax_id
                line = lookup_file_opened.readline()
            lookup_file_opened.close()

        return ncbi_to_taxid_dict


    # TODO: complete this. @Latifa: what was this for?? Results already has a similar method
    def calculate_outcomes(self):
        '''
        DESCRIPTION:
            Finds true positives, false positives, and false negatives
        '''
        for result_obj in self.result_objs:
            break

    
    def write_intermediate_csv(self, output_file) -> str:
        '''
        DESCRIPTION:
            Writes all data of the result objects to a CSV file
            i.e. [trial_num, experiment, condition, tool, true_positives, true_positives, false_negatives, precision, recall, f1, l2]
        
        INPUT:
            1. The file to write to
        
        OUTPUT:
            Returns the CSV file path (string) in which the data was stored
        '''

        all_result_metrics = [[]]

        with open(BENCHMARKING_OUTPUT_PATH + output_file, 'w') as csvfile: 
            writer = csv.writer(csvfile) 
            writer.writerows(all_result_metrics)

        return BENCHMARKING_OUTPUT_PATH + output_file


    def write_to_csv(self, output_file):
        """
        creates an output csv file
        """
        for result_object in self.result_objs:
            print(result_object.calc_f1score())



    
    def validate_metrics():
        '''
        DESCRIPTION:
            Checks if there are any anomolies in terms of exactly equal precision and recall.

        INPUT:


        OUTPUT:
            Returns the experiment conditions that were flagged as exactly equal.
            Dict of table key = occurrence #, value = table of 'result' data that share these numbers 
                (List of lists, each row in the list is a list of result info)
                or Dict of pandas.DataFrame?

        '''
        pass



if __name__ == "__main__":
    benchmark_obj = Benchmarking('/Users/latifa/GitHub/benchmarking-enrichseq/benchmarking_metadata.csv')
    ncbi2taxid_dict = benchmark_obj.ncbi_to_taxid_mapping()
    
    truth_obj = SimulatedTruth(2, 'num_reads_200genomes', '200_genomes_400000_reads', '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/test_2/num_reads_200genomes/200_genomes_200000_reads.fa')
    truth_obj.parse_simulated_fasta(ncbi2taxid_dict)

    #benchmark_obj.assign_objs() # assign truth objects to result objects
    
    fve_obj = FVEResult(2, 'num_reads_200genomes', '200_genomes_400000_reads', 'FastViromeExplorer', '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/test_2/FastViromeExplorer/num_reads_200genomes/200_genomes_200000_reads/FastViromeExplorer-final-sorted-abundance.tsv', truth_obj,ncbi2taxid_dict)
    print(fve_obj.parse_results()) # works

    # validate with less number of genomes
    truth_obj2 = SimulatedTruth(1, 'num_genomes', '10_genomes_500000_reads', '/Users/latifa/GitHub/benchmarking-enrichseq/tests-ALL/test_1/num_genomes/10_genomes_500000_reads.fa')
    fve_obj2 = FVEResult(1, 'num_genomes', '10_genomes_500000_reads', 'FastViromeExplorer', '/Users/latifa/GitHub/benchmarking-enrichseq/results_simulated/test_1/FastViromeExplorer/num_genomes/10_genomes_500000_reads/FastViromeExplorer-final-sorted-abundance.tsv', truth_obj2,ncbi2taxid_dict)
    print(truth_obj2.abundances_dict)





