### NOTE: Currently needs to be run in the same directory as the simulate_reads.py script
### because of some dependencies

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

	# some input validation
	if [[ ! -d ${configDir} ]]; then
		echo "Not running generateSimulations(). Config directory doesn't exist.";
	elif [[ ! -f ${simScriptPath} ]]; then
		echo "Not running generateSimulations(). Simulation script doesn't exist.";
	else
		echo "Running generateSimulations()";
		if [[ ! -d ${simOutDir} ]]; then
			mkdir ${simOutDir}
		fi
		for config in ${configDir}/*.config
		do
			# extract filename sans extension
			filename=$(basename -s .config ${config})
			# run simulator script for each config file
			python ${simScriptPath} --quiet -i ${config} -c 30 -o ${simOutDir}/${filename} -t 4
		done
		echo "Simulation files generated:"
		ls ${simOutDir}
	fi
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
# function runEnrichSeq() {
#
# }

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
	local simScriptPath="/Users/latifa/GitHub/Phage-EnrichSeq/readsimulator_module/simulate_reads.py"
	local configDir="/Users/latifa/Genomics/simulation_experiments_march28/sim_configs"
	local simOutDir="/Users/latifa/GitHub/Phage-EnrichSeq/analysis_module/simulations"

	# running underlying methods
	generateSimulations ${simScriptPath} ${configDir} ${simOutDir};

}

echo "Running the BLASTDB build script";
main;
