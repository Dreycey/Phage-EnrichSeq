

### because of some dependencies (i.e. cd to readSimulatorModule then run ../analysis_module/automate_testing.sh)

ANALYSIS_PATH="/Users/latifa/GitHub/Phage-EnrichSeq/analysis_module"
ENRICHSEQ_PATH="/Users/latifa/GitHub/Phage-EnrichSeq/"
SIM_BENCHMARK_PATH="${ENRICHSEQ_PATH}/benchmarking"
if [[ -d ${SIM_BENCHMARK_PATH} ]]; then
	rm -rf ${SIM_BENCHMARK_PATH}
fi

mkdir ${SIM_BENCHMARK_PATH}


# 1. automate simulations
#######################################
# generateSimulations()
# creates simulation FASTA files using the simulate_reads.py in the readSimulatorModule
# Globals:
#   None
# Arguments:
#   directory path of config files
#   output directory name for resulting simulations
# Outputs:
#   creates simulation FASTA files in the specified directory
#######################################
function generateSimulations() {
	#  arguments
	local simScriptPath=$1;
 	local configDir=$2;
	local simOutDir=$3;
	local threads=$4;


	# some input validation
	if [[ ! -d ${configDir} ]]; then
		echo "Not running generateSimulations(). Config directory doesn't exist.";
		exit 1
	fi
	if [[ ! -f ${simScriptPath} ]]; then
		echo "Not running generateSimulations(). Simulation script doesn't exist.";
		exit 2
	fi

	# run script
	echo "Running generateSimulations()";
	if [[ ! -d ${simOutDir} ]]; then
		mkdir ${simOutDir}
	fi
	for config in ${configDir}/*.config
	do
		# extract filename sans extension
		filename=$(basename -s .config ${config})
		# run simulator script for each config file
		python ${simScriptPath} --quiet -i ${config} -c 30 -o ${simOutDir}/${filename} -t ${threads}
	done
	echo "Simulation files generated in ${simOutDir}:"
	ls ${simOutDir}


}

# 2. automate nextflow runs
#######################################
# runEnrichSeq()
# runs enrichseq workflows on all simulated files in specified directory
# Globals:
#   None
# Arguments:
#   directory path of simulated FASTA files
#
# Outputs:
#   creates working directories containing output files of all modules in enrichseq
#######################################

function runEnrichSeqIllumina() {
	#  arguments
	local enrichseqPath=$1;
	local simOutDir=$2;
	local krakenDbPath=$3;
	local threads=$4;


	# input validation
	if [[ ! -f ${enrichseqPath} ]]; then
		echo "Not running runEnrichSeq(). Enrichseq script or path does not exist.";
		exit 3
	fi


	if [[ -z ${krakenDbPath} ]]; then
		echo "Not running runEnrichSeq(). Database path not provided.";
		exit 4
	fi

	echo "Running runEnrichSeqIllumina()";
	for simfile in ${simOutDir}/*_illumina.fa
	do
		# extract filename sans extension
		filename=$(basename -s .fa ${simfile})
		# run enrichseq.nf for each simulation FASTA file
		nextflow ${enrichseqPath} --read single --fasta ${simfile} \
				--workdir ${SIM_BENCHMARK_PATH}/${filename} --dbdir ${krakenDbPath} --threads ${threads}
	done
	echo "EnrichSeq output directories located in ${SIM_BENCHMARK_PATH}"

}

# 3. save enrichseq results to CSV
#######################################
# saveResults()
# extracts results from subdirectories and stores in CSV file
# Globals:
#   None
# Arguments:
#
#
# Outputs:
#
#######################################

function saveResults() {
	# arguments
	local phageNamesFile=$1; # a file containing the names of all phages simulated
	local configDir=$2;
	local resultsDir=$3;

	if [[ ! -f ${phageNamesFile} ]]; then
		echo "Not running saveResults(). Please provide a file with the names of tested phages."
		exit 5
	fi

	if [[ ! -d ${configDir} ]]; then
		echo "Not running saveResults(). Config directory doesn't exist."
		exit 6
	fi

	echo "Running saveResults()"

	if [[ -d ${resultsDir} ]]; then
		rm -rf ${resultsDir}
	fi
	mkdir ${resultsDir}

	# run for all config files in the config directory
	for config in ${configDir}/*.config
	do
		# extract filename sans extension
		filename=$(basename -s .config ${config})


		if [[ ! -d ${SIM_BENCHMARK_PATH}/${filename}_illumina ]]; then
			echo "${filename}_illumina directory does not exist."
			exit 7
		elif [[ ! -f ${SIM_BENCHMARK_PATH}/${filename}_illumina/enrichseq/bracken/bracken_run_orig.bracken ]]; then
			echo "Bracken output does not exist."
			exit 8
		else
			# run simulation analysis python script to save to CSV file
			python ${ANALYSIS_PATH}/simulation_analysis.py ${phageNamesFile} ${config} \
					${SIM_BENCHMARK_PATH}/${filename}_illumina/enrichseq/bracken/*.bracken \
					${resultsDir}/${filename}_analysis.csv
		fi
	done
	echo "CSV files for simulations stored in ${resultsDir}."


#######################################
# Main function runs all of the other methods for the testing
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Calls on functions for automating testing on simulations
#######################################
function main(){
	# input arguments

	local simScriptPath="${ENRICHSEQ_PATH}/readsimulator_module/simulate_reads.py"
	local configDir="/Users/latifa/Genomics/simulation_experiments_march28/sim_configs"
	local simOutDir="${SIM_BENCHMARK_PATH}/simulations"
	local enrichseqPath="${ENRICHSEQ_PATH}/Nextflow/enrichseq.nf"
	local krakenDbPath="${ENRICHSEQ_PATH}/kraken_module/krakenDB"
	local phageNamesFile="${ANALYSIS_PATH}/phage_names.txt"
	local resultsDir="${SIM_BENCHMARK_PATH}/simulation_CSVs"
	local threads=4

	# running underlying methods
	generateSimulations ${simScriptPath} ${configDir} ${simOutDir} ${threads}
	runEnrichSeqIllumina ${enrichseqPath} ${simOutDir} ${krakenDbPath} ${threads}
	saveResults ${phageNamesFile} ${configDir} ${resultsDir}
}

echo "Running the simulation testing script";
main;
