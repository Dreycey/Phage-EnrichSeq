
from pathlib import Path
from typing import Tuple, List, Dict

# Takes CSV file with result files

class Truth:
    # trial_num = 1
    # experiment_name = ''
    # experiment_condition = ''
    # num_genomes = 1
    # num_reads = 1000000

    # Truth object should only be created if it actually exists in the metadata CSV
    def __init__(self, trial_num, experiment_name, experiment_condition, sequencing_file) -> None:
        self.trial_num = trial_num
        self.experiment_name = experiment_name
        self.experiment_condition = experiment_condition
        self.sequencing_file = sequencing_file
        self.num_genomes, self.num_reads = self.extract_genomes_and_reads(self.experiment_condition)

    
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
    

    def parse_simulated_fasta(self) -> Dict[int, float]:
        '''
        DESCRIPTION:
            Parses 'truth' sequencing file and stores tax IDs and corresponding abundances
            in a dictionary
        INPUT:
            -
        OUTPUT:
            Dictionary where <key : value> pair is taxid : abundance
        '''
        simulated_data_dict = {}
        # Read from simulated fasta file (member variable)
        # Populate dict with taxids and abundances found in the fasta file
        return simulated_data_dict


class Result:
    # tool_name = ''
    # trial_num = 1
    # experiment_name = ''
    # experiment_condition = ''
    # num_genomes = 1
    # num_reads = 1000000
    #corresponding_truth_obj = Truth()

    def __init__(self) -> None:
        pass
    

    # The following metrics require a truth object to be associated with 
    def calc_precision(self) -> float:
        pass

    def calc_recall(self) -> float:
        pass

    def calc_f1score(self) -> float:
        pass

    def calc_l2distance(self) -> float:
        pass

    def write_to_csv(self) -> str:
        ''' 
        writes benchmarking results to a CSV file (precision, recall, f1 scores and l2 distance)
        '''

class EnrichSeqResult(Result):
    def parse_results(self):
        pass


class FVEResult(Result):
    def parse_results(self):
        pass


class BrackenResult(Result):
    def parse_results(self):
        pass


class Benchmark:
    '''
    Each Result object should be associated with a Truth object 
    '''
    # benchmark_num = 0
    # benchmark_title = ''

    
    def parse_truths_metadata() -> list:
        '''
        parse metadata CSV and store in Truth objects
        '''
        truth_objs = []
        return truth_objs


    def parse_results_metadata() -> list:
        '''
        parse CSV and store in Result objects
        '''
        result_objs = []
        return result_objs

    



# Parse FastViromeExplorer
def parse_fastvirome(result_file_name: str):
    fve_result_dict = {}
   # fve_path = Path('/projects/laal5512/results/test_1/FastViromeExplorer/num_genomes/10_genomes_500000_reads/')

    return fve_result_dict




if __name__ == "__main__":
    fve_dict = parse_fastvirome("FastViromeExplorer-final-sorted-abundance.tsv")



