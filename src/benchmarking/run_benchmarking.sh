
usage() {
    echo; echo "Usage: bash $0 -i <path/to/truth/files> -o <path/to/desired/output/loc>"
    echo "  -i      Path to simulated test files"
    echo "  -o      Desired location to generate output files"
    echo "  -h, --help      Print this help message out"; echo;
    exit 1;
}


while getopts ":i:t:d:o:" flag; do
    case "${flag}" in
        i)
            truth_dir=${OPTARG}
            ;;
        o)
            output_dir=${OPTARG}
            ;;
    esac
done

function test_bracken_assembled() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="bracken_assembled";
    local file_suffix=".fa";
    local kraken_db="/projects/laal5512/tools/Phage-EnrichSeq/database/krakenDB/";

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            test_name="$(basename -- $test_dir)";
            if [[ -d $test_dir ]]; then
                for fasta_file in $test_dir/*${file_suffix}; do
                    echo; echo $fasta_file; echo;
                    basefile="$(basename -- $fasta_file)";
                    mkdir -p ${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}};

                    # Define output locations for megahit, kraken and bracken
                    megahit_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/megahit";
                    kraken_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/kraken";
                    mkdir ${kraken_output_dir};
                    bracken_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/bracken";
                    mkdir ${bracken_output_dir};

                    # 1. Run megahit
                    echo "-- Running MegaHIT --";
                    megahit -r ${fasta_file} -t 4 -m 1e9 -o ${megahit_output_dir} --out-prefix megahit_out;


                    # 2. Run Kraken2
                    echo "-- Running Kraken2 --";
                    kraken2 --use-names --threads 4 --db ${kraken_db} \
                                --report ${kraken_output_dir}/kraken.report \
                                ${megahit_output_dir}/megahit_out.contigs.fa > ${kraken_output_dir}/kraken2_benchmarking.log

                    # 3. Run Bracken
                    echo "-- Running Bracken --";
                    bracken -d ${kraken_db} \
                            -i ${kraken_output_dir}/kraken.report \
                            -o ${bracken_output_dir}/assembled_abundances.bracken;

                done
            fi
        done
    done
}


function test_bracken_no_assembly() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="bracken_no_assembly";
    local file_suffix=".fa";
    local kraken_db="tools/Phage-EnrichSeq/database/krakenDB/";

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            test_name="$(basename -- $test_dir)";
            if [[ -d $test_dir ]]; then
                for fasta_file in $test_dir/*${file_suffix}; do
                    echo; echo $fasta_file; echo;
                    basefile="$(basename -- $fasta_file)";
                    mkdir -p ${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}};

                    # Define output locations for megahit, kraken and bracken
                    megahit_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/megahit";
                    kraken_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/kraken";
                    mkdir ${kraken_output_dir};
                    bracken_output_dir="${result_dir}/${trial}/${tool_name}/${test_name}/${basefile%${file_suffix}}/bracken";
                    mkdir ${bracken_output_dir};

                    # 1. Run Kraken2
                    echo "-- Running Kraken2 --";
                    kraken2 --threads 4 --db ${kraken_db} \
                                --report ${kraken_output_dir}/kraken.report \
                                ${fasta_file} > ${kraken_output_dir}/kraken2_benchmarking.log

                    # 2. Run Bracken
                    echo "-- Running Bracken --";
                    bracken -d ${kraken_db} \
                            -i ${kraken_output_dir}/kraken.report \
                            -o ${bracken_output_dir}/nonassembled_abundances.bracken;

                done
            fi
        done
    done
}


function test_enrichseq() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="enrichseq";
    local file_suffix=".fa";

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            if [[ -d $test_dir ]]; then
                for file in $test_dir/*${file_suffix}; do
                echo; echo $file; echo;
                basefile="$(basename -- $file)";
                test_dir_name="$(basename "${test_dir##*/}")";
                mkdir -p ${result_dir}/${trial}/${tool_name}/${test_dir_name}/;
                python3 tools/Phage-EnrichSeq/EnrichSeq.py enrichseq \
                            -1 $file \
                            -o ${result_dir}/${trial}/${tool_name}/${test_dir_name}/${basefile%${file_suffix}};
                done
            fi
        done
    done
}

function test_viromexplorer() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="FastViromeExplorer";
    local file_suffix=".fq";

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            if [[ -d $test_dir ]]; then
                for file in $test_dir/*${file_suffix}; do
                echo; echo $file; echo;
                basefile="$(basename -- $file)";
                test_dir_name="$(basename "${test_dir##*/}")";
                mkdir -p ${result_dir}/${trial}/${tool_name}/${test_dir_name}/;
                java -cp tools/FastViromeExplorer/bin FastViromeExplorer \
                                                            -1 $file \
                                                            -i tools/FastViromeExplorer/ncbi-virus-kallisto-index-k31.idx \
                                                            -o ${result_dir}/${trial}/${tool_name}/${test_dir_name}/${basefile%${file_suffix}} \
                                                            -l tools/FastViromeExplorer/ncbi-viruses-list.txt; 
                done
            fi
        done
    done
}

function test_kraken2() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="Kraken2";
    local file_suffix=".fa";
    

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            if [[ -d $test_dir ]]; then
                for file in $test_dir/*${file_suffix}; do
                echo; echo $file; echo;
                basefile="$(basename -- $file)";
                test_dir_name="$(basename "${test_dir##*/}")";
                mkdir -p ${result_dir}/${trial}/${tool_name}/${test_dir_name}/;
                kraken2 --use-names --threads 4 --db tools/${tool_name}/minikraken2_v2_8GB_201904_UPDATE \
                            --report ${result_dir}/${trial}/${tool_name}/${test_dir_name}/${basefile%${file_suffix}}.report \
                            ${file} > kraken2_benchmarking.log
                                
                done
            fi
        done
    done
}


function test_bracken() {
    local truth_dir=$1;
    local result_dir=$2;
    local tool_name="Bracken";
    local dependency_tool_name="Kraken2";
    local file_suffix=".fa";

    echo "TESTING ${tool_name}";
    for trial_dir in $truth_dir/*; do
        trial="$(basename -- $trial_dir)";
        for test_dir in $trial_dir/*; do
            if [[ -d $test_dir ]]; then
                for file in $test_dir/*${file_suffix}; do
                echo; echo $file; echo;
                basefile="$(basename -- $file)";
                test_dir_name="$(basename "${test_dir##*/}")";
                mkdir -p ${result_dir}/${trial}/${tool_name}/${test_dir_name}/;

                bracken -d ${truth_dir}/${dependency_tool_name}/minikraken2_v2_8GB_201904_UPDATE \
                            -i ${result_dir}/${trial}/${dependency_tool_name}/${test_dir_name}/${basefile%${file_suffix}}.report \
                            -o ${result_dir}/${trial}/${tool_name}/${test_dir_name}/${basefile%${file_suffix}}.bracken;
                done
            fi
        done
    done
}

function main() {
    local result_dir=${output_dir}/results2/;
    mkdir ${result_dir};
    #conda activate enrichseq;
    #test_enrichseq;
    #test_kraken2;
    #test_bracken;
    #test_viromexplorer;
    test_bracken_assembled ${truth_dir} ${result_dir}; 
    test_bracken_no_assembly ${truth_dir} ${result_dir}; 
    
}
main;
