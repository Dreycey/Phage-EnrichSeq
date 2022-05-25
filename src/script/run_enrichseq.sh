

function test_enrichseq() {
    path=$1;
    results_dir=$2;
    tool_name="enrichseq";
    file_suffix=".fa";
    #echo "TESTING ${tool_name}";
    for test_dir in $path*; do
        if [[ -d $test_dir ]]; then
            for file in $test_dir/*${file_suffix}; do
               #echo; echo $file; echo;
               basefile="$(basename -- $file)";
               test_dir_name="$(basename "${test_dir##*/}")";
               mkdir -p results/${tool_name}/${test_dir_name}/;
               #echo results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}};
               if [[ ! -d results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}}  ]]; then
                   echo "python3 ./EnrichSeq.py enrichseq -1 $file -o ${results_dir}/${tool_name}/${test_dir_name}/${basefile%${file_suffix}};";
               fi
             done
    	fi
    done
}

function main() {
    #mkdir results/;
    #conda activate enrichseq;
    test_enrichseq $1 $2;
}
main $1 $2
