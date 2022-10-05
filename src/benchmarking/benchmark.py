import sys
import csv
import argparse
from abc import ABC, abstractmethod
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Dict
from dynaconf import Dynaconf


#NCBI_TO_TAXID_FILEPATH = '/Users/latifa/GitHub/benchmarking-enrichseq/out_dictionary.txt'
PHAGE_REF_DB_PATH = '/Users/latifa/GitHub/benchmarking-enrichseq/tools/Phage-EnrichSeq/database/ref_genomes/'
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
                    if ncbi_id in ncbi_to_taxid_dict:
                        tax_id = ncbi_to_taxid_dict[ncbi_id]
                    else: 
                        # TODO: write to log
                        print(f"NCBI ID {ncbi_id} not in ncbi2taxid dictionary")
                        file_name = line[1:].split("|")[-1].strip("\n")
                        try:
                            path = Path(settings.PHAGE_REF_DB_PATH + file_name)
                        except AttributeError as e:
                            print(f"Unable to open 'PHAGE_REF_DB_PATH': {e}", file=sys.stderr)
                        else:
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


class Result(ABC):

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

        # Run methods for creating data structure (order matters here!)
        self.truth_obj.parse_simulated_fasta(self.ncbi_to_taxid_dict)
        self.parse_results()
        self.calc_classification_outcomes()


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


    def _clean_abundances(self):
        '''
        DESCRIPTION:
            Cleans the Result object's abundances dictionary of any negligible 
            abundances (based on a threshold)
        '''
        ids_to_del = []
        for tax_id in self.abundances_dict:
            if self.abundances_dict[tax_id] < 0.0001:
                ids_to_del.append(tax_id)

        # delete 0 value tax ids
        for val in ids_to_del:
            del self.abundances_dict[val]

        return self.abundances_dict


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
        if self.truth_obj is None:
            print(f'{type(self).__name__} object [{self.trial_num}, {self.num_genomes}, {self.num_reads}] does not have a truth object assigned to it.')
        elif not self.abundances_dict:
            print(f'{type(self).__name__} object [{self.trial_num}, {self.num_genomes}, {self.num_reads}] does not have TAXID : ABUNDANCE dictionary.')
        else:
            for pred_taxid in self.abundances_dict.keys():
                if pred_taxid in self.truth_obj.abundances_dict.keys():
                    self.true_positives.append(pred_taxid)
            
            for pred_taxid in  self.abundances_dict.keys():
                if (pred_taxid not in self.truth_obj.abundances_dict.keys()) and (pred_taxid != "UNK"):
                    self.false_positives.append(pred_taxid)

            for pred_taxid in self.truth_obj.abundances_dict.keys():
                if pred_taxid not in self.abundances_dict.keys():
                    self.false_negatives.append(pred_taxid)

            print(f'{len(self.true_positives)}, {len(self.false_positives)}, {len(self.false_negatives)}')
            return self.true_positives, self.false_positives, self.false_negatives


    def calc_precision(self) -> float:
        try:
            precision = len(self.true_positives) / (len(self.true_positives) + len(self.false_positives))
        except ZeroDivisionError:
            print('Division by zero caught for "precision" metric.')
            return 0.0
        else:
            return precision


    def calc_recall(self) -> float:
        try:
            recall = len(self.true_positives) / (len(self.true_positives) + len(self.false_negatives))
        except ZeroDivisionError:
            print('Division by zero caught for "recall" metric.')
            return 0.0
        else:
            return recall


    def calc_f1score(self) -> float:
        precision = self.calc_precision()
        recall = self.calc_recall()
        #print(f'Precision = {precision}, Recall = {recall}')
        try:
            f1 = 2 * (precision * recall) / (precision + recall)
        except ZeroDivisionError:
            print(f'Zero division in f1 score calculation for {self.tool_name}.')
            return 0.0
        else:
            return f1 


    def calc_l2distance(self) -> float:
        l2_distance = 0.0
        return l2_distance


    @abstractmethod
    def parse_results(self):
        """
        parses the file of a particular type
        """
        pass


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
                self.abundances_dict[tax_id] = float(abundance_val)
            enrichseq_file_opened.close()

        return self.abundances_dict


class FVEResult(Result):
    delim = '\t'
    ncbi_id_column = 0
    count_column = 3


    def parse_results(self, threshold=0.0) -> Dict:
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
                elif ncbi_id.split(".")[0] in self.ncbi_to_taxid_dict: 
                    tax_id = self.ncbi_to_taxid_dict[ncbi_id.split(".")[0]]
                else:
                    tax_id = 'UNK'
                    if float(counts_val) > threshold :
                        print(f"ERROR (parsing):")
                        print(f"\tTool: FVE")
                        print(f"\t\tFVE: {ncbi_id} NOT CONVERTED")
                        print(f"\t\tabundances: {float(counts_val)}") 
                # Calculate 
                self.abundances_dict[tax_id] = float(counts_val)
                total_abundance += float(counts_val)
            fve_file_opened.close()

        for tax_id in self.abundances_dict:
            self.abundances_dict[tax_id] /= total_abundance # Normalize
        self._clean_abundances()
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
                self.abundances_dict[tax_id] = float(abundance_val)
            bracken_file_opened.close()
        
        self._clean_abundances()
        return self.abundances_dict


class Benchmarking:

    def __init__(self, metadata_csv = None) -> None:
        self.metadata_csv = metadata_csv
        self.result_objs = []
        self.ncbi_to_taxid_dict: dict = self.ncbi_to_taxid_mapping()
        # init metadata
        self.metadata_df = self.parse_metadata_csv()
        self.assign_objs()
    

    def parse_metadata_csv(self) -> pd.DataFrame:
        """
        Description:
            This parses the input csv into a pandas Dataframe.
        Errors:
            returns empty dataframe
        """
        try:
            metadata_df = pd.read_csv(self.metadata_csv, skipinitialspace=True)
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
        
        # filter matching rows
        # filtered_data = self.metadata_df.loc[(self.metadata_df["Trial_Num"] == trial_num_temp) &
        #                                          (self.metadata_df["Experiment"] == experiment_temp) &
        #                                          (self.metadata_df["Condition"] == condition_num_temp)]

        # grab unique triplicates for truth-only rows
        truth_rows = self.metadata_df.loc[self.metadata_df["Tool"] == "truth"]
        unique_tests = truth_rows[["Trial_Num", "Experiment", "Condition"]].drop_duplicates()
        #print(f"unique: {unique_tests}")
        for index, unique_test_set in unique_tests.iterrows():
            trial_num_temp, experiment_temp, condition_num_temp = (unique_test_set.Trial_Num, 
                                                                   unique_test_set.Experiment, 
                                                                   unique_test_set.Condition)
            # filter matching rows
            filtered_data = self.metadata_df.loc[(self.metadata_df["Trial_Num"] == trial_num_temp) &
                                                 (self.metadata_df["Experiment"] == experiment_temp) &
                                                 (self.metadata_df["Condition"] == condition_num_temp)]
            print(f"filtered data: {filtered_data}")
            # create a truth object
            truth_row_temp = filtered_data.loc[filtered_data["Tool"] == "truth"].drop(columns=['Tool']).values.tolist()
            # if not truth_row_temp:
            #     print("TRUTH IS EMPTY!")
            #     exit()
            #print(f"\n {truth_row_temp}")
            truth_obj = SimulatedTruth(*truth_row_temp[0]) 
            for index2, sub_row in filtered_data.iterrows():
                tool_name = sub_row["Tool"]
                if (tool_name != "truth"):
                    result_obj = self._return_tool_obj(tool_name, sub_row.values.tolist(), truth_obj, self.ncbi_to_taxid_dict)
                    self.result_objs.append(result_obj)


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
        #ncbi_map_filepath = Path(BENCHMARKING_OUTPUT_PATH + 'ncbi_map.tab')
        ncbi_map_filepath = Path(settings.NCBI_TO_TAXID_PATH)

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

    
    def write_to_csv(self, output_file) -> str:
        '''
        DESCRIPTION:
            Writes all data of the result objects to a CSV file
            i.e. [trial_num, experiment, condition, tool, true_positives, true_positives, false_negatives, precision, recall, f1, l2]
        
        INPUT:
            The file to write to
        
        OUTPUT:
            Returns the CSV file path (string) in which the data was stored
        '''

        all_result_metrics = [[]]

        if self.result_objs:
            for res_obj in self.result_objs:
                row = [res_obj.tool_name, res_obj.trial_num, res_obj.experiment_name, res_obj.experiment_condition]
                row.extend([len(res_obj.true_positives), len(res_obj.false_positives), len(res_obj.false_negatives),
                                            res_obj.calc_precision(), res_obj.calc_recall(), res_obj.calc_f1score(), res_obj.calc_l2distance()])


                all_result_metrics.append(row)
        try:
            csvfile = open(settings.OUTPUT_DIRECTORY + output_file, 'w')
        except AttributeError as e:
            print(f"Unable to find 'OUTPUT_DIRECTORY' in settings: {e}", file=sys.stderr)
        else:
            writer = csv.writer(csvfile) 
            writer.writerows(all_result_metrics)
            csvfile.close()
            return settings.OUTPUT_DIRECTORY + output_file


def metadata_passes(metadata_csv, output_file):
    """
    validates the metadata file and creates a log

    Checks:
        1. Number of unique tripicate columns are equal to 1+#tools_tested
            i.e. 4 (truth, enrichseq, bracken, FVE)
            Note: adds missing files to log. 
        2. 
    """
    return True





def parseArgs(argv=None) -> argparse.Namespace:
    '''
    DESCRIPTION:
        This method takes in the arguments from the command and performs
        parsing.
    INPUT: 
        Array of input arguments
    OUTPUT:
        returns a argparse.Namespace object
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("-f", "--metadata_file", help="metadata CSV file path containing experiment information", required=True)
    parser.add_argument("-o", "--output_prefix", help="prefix for final csv file", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    settings = Dynaconf(settings_files="settings.toml")
    arguments = parseArgs(argv=sys.argv[1:])
    # if len(sys.argv) < 3:
    #     print(f"USAGE:")
    #     print(f"python3 benchmark.py <metadata CSV> <output PREFIX>")
    if metadata_passes(arguments.metadata_file, arguments.output_prefix + ".log"):
        benchmark_obj = Benchmarking(arguments.metadata_file)
        benchmark_obj.write_to_csv(arguments.output_prefix + ".csv")








