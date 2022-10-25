

usage() {
    echo; echo "Usage: bash $0 -i -t -o -d"
    echo "  -i      Path to simulated test files"
    echo "  -t      Path to tools (megahit, kraken, bracken)"
    echo "  -d      Path to the kraken2 database created using kraken2Build.sh"
    echo "  -o      Desired location to generate output files"
    echo "  -h, --help      Print this help message out"; echo;
    exit 1;
}


while getopts ":i:t:d:o:" flag; do
    case "${flag}" in
        i)
            truth_dir=${OPTARG}
            ;;
        t)
            tools_path=${OPTARG}
            ;;
        d)
            db_dir=${OPTARG}
            ;;
        o)
            output_dir=${OPTARG}
            ;;
    esac
done

# if [ -z "${truth_dir}" ] || [ -z "${tools_path}" ] || [ -z "${db_dir}" ]; then
#     usage
# fi


function initialize_files() {
    local output_path=$1;

    if [[ -d ${output_path} ]]; then
        if [[ -d ${output_path}/Bracken_Results/ ]]; then
            rm -rf ${output_path}/Bracken_Results/;
        fi
        mkdir ${output_path}/Bracken_Results/;
    else 
        output_path=$(pwd);
    fi

    mkdir ${output_path}/Bracken_Results/;
}


function run_megahit() {
    local tool_path=$1;
    local input_file=$2;
    local output_dir=$3;

    echo "-- RUNNING MegaHIT --";

    if [[ -d ${truth_dir} ]]; then
        #echo ${input_file};
        # bash ${tool_path}/megahit_module/megahitRun.sh --read="single" \
        # 		  --input1=${input_file} \
    	# 		  --threads="4" \
    	# 		  --out=${output_dir};
        megahit -r ${input_file} -t 4 -m 1e9 -o ${output_dir} --out-prefix megahit_out 
    else
        echo "Test directory ${truth_dir} does not exist."
    fi

}


function run_kraken2() {
    local assembled=$1;
    local dbDir=$2;
    local inFasta=$3;

    
    if [[ -d ./kraken_output ]]; then
        mkdir "./kraken_output"
    fi
    local krakenOutFile='.kraken_out/K${kmerLength}L${minimizerLength}S${seed}.kraken'

    if [[ -f ${dbDir}/taxo.k2d ]]; then
        kraken2 --use-names --threads 4 --db ${dbDir} --report kraken.report ${inFasta} > ${krakenOutFile};
    else
        echo "Kraken database must be built first"
    fi
}


function run_bracken() {
    local assembled=$1;
}


function run_pipeline() {
    local truth_dir=$1;
    local result_dir=$2;
    local assembled=$3;

    for trial_dir in $truth_dir/*; do
        trial=$(basename -- ${trial_dir});
        for test_dir in $trial_dir/*; do
            test=$(basename -- ${test_dir});
            for fasta_file in $test_dir/*.fa; do
                # 1. Run MegaHIT
                test_condition="$(basename -- $fasta_file)";
                output_dir="${result_dir}/${trial}/Bracken_Assembled/${test}";
                mkdir -p ${output_dir};
                megahit_output_dir="${output_dir}/${test_condition%.fa}";
                run_megahit ${tools_path} ${fasta_file} ${megahit_output_dir};
            done
        done
    done
}


function main() {
    
    initialize_files ${output_dir};
    run_pipeline ${truth_dir} "${output_dir}/Bracken_Results/" 1;
    #run_pipeline ${truth_dir} ${output_dir}Bracken_Results/ 0 &;
    
}

echo "Running Bracken pipeline";
main;