bash ${params.toolpath}/megahit_module/megahitRun.sh --read=${params.read} \
        		  --input1=${fastafile} \
                  --input2=${fastafile_2} \
    			  --threads=${THREADS} \
    			  --out=${megahitDir}

usage() {
    echo; echo "Usage: bash $0 --input=simulatedgenomes_illumina.fa --refdir=database/ref_genomes/ --kmer_length=35 --min_length=31 --seed=7"
    echo "  --input_dir     Path to simulated test files"
    echo "  --tools_path    Path to tools (megahit, kraken, bracken)"
    echo "  --kraken_db     Path to the kraken2 database created using the config and kraken2Build.sh"
    echo "  --output_path   Desired location to generate output files"
    echo "  -h, --help      Print this help message out"; echo;
    exit 1;
}


function initialize_files() {
    local outputPath=$1;

    if [[ -d ${outputPath} ]]; then
        if [[ -d ${outputPath}/Bracken_Results/ ]]; then
            rm -rf ${outputPath}/Bracken_Results/;
        mkdir ${outputPath}/Bracken_Results/;
    else 
        outputPath=$(pwd);
    fi

    mkdir ${outputPath}/Bracken_Results/;

    return "${outputPath}/Bracken_Results/"
}


function run_megahit() {
    local toolPath=$1;
    local testDir=$2;
    local file_suffix=".fa"
    echo "-- RUNNING MegaHIT --";

    if [[ -d ${testDir} ]]; then
        for dir in ${testDir}/*; do
            if [[ -d $dir ]]; then
                for file in $test_dir/*${file_suffix}; do
                echo; echo $file; echo;
                basefile="$(basename -- $file)";
                test_dir_name="$(basename "${test_dir##*/}")";
                mkdir -p results/${tool_name}/${test_dir_name}/;
                kraken2 --use-names --threads 4 --db tools/${tool_name}/minikraken2_v2_8GB_201904_UPDATE \
                            --report results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}}.report \
                            ${file} > kraken2_benchmarking.log
                                
                done
            fi
        done
    else
        echo "Test directory ${testDir} does not exist."
    fi

}


function run_kraken2() {
    local assembled=$1;
    local dbDir=$2;
    local inFasta=$3;
    local kmerLength=$4;
    local minimizerLength=$5;
    local seed=$6;
    
    if [[ -d ./kraken_output ]]; then
        mkdir "./kraken_output"
    fi
    local krakenOutFile='.kraken_out/K${kmerLength}L${minimizerLength}S${seed}.kraken'

    if [[ -f ${dbDir}/taxo.k2d ]]; then
        kraken2 --use-names --threads 4 --db ${dbDir} --report kraken.report ${inFasta} > ${krakenOutFile}
    else
        echo "Kraken database must be built first"
    fi
}


function run_bracken() {
    local assembled=$1;
}


function main() {
    local dbDir='krakenDB';

    buildDatabase ${dbDir}; ## Which params to pass?
    runMegahit ${};
    runKraken2 ${dbDir} ${inFasta} ${kmerLength} ${minimizerLength} ${seed};
}

echo "Running Bracken pipeline";
main;