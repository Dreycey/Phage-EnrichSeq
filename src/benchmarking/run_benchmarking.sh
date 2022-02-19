

function test_enrichseq() {
    tool_name="enrichseq2";
    file_suffix=".fa";
    echo "TESTING ${tool_name}";
    for test_dir in ./tests/*; do
        if [[ -d $test_dir ]]; then
            for file in $test_dir/*${file_suffix}; do
               echo; echo $file; echo;
               basefile="$(basename -- $file)";
               test_dir_name="$(basename "${test_dir##*/}")";
               mkdir -p results/${tool_name}/${test_dir_name}/;
           	   python3 tools/Phage-EnrichSeq/EnrichSeq.py enrichseq \
                           -1 $file \
                           -o results2/${tool_name}/${test_dir_name}/${basefile%${file_suffix}};
    	    done
    	fi
    done
}

function test_viromexplorer() {
    tool_name="FastViromeExplorer2";
    file_suffix=".fq";
    echo "TESTING ${tool_name}";
    for test_dir in ./tests/*; do
        if [[ -d $test_dir ]]; then
            for file in $test_dir/*${file_suffix}; do
               echo; echo $file; echo;
               basefile="$(basename -- $file)";
               test_dir_name="$(basename "${test_dir##*/}")";
               mkdir -p results/${tool_name}/${test_dir_name}/;
           	   java -cp tools/FastViromeExplorer/bin FastViromeExplorer \
                                                        -1 $file \
                                                        -i tools/FastViromeExplorer/ncbi-virus-kallisto-index-k31.idx \
                                                        -o results2/${tool_name}/${test_dir_name}/${basefile%${file_suffix}} \
                                                        -l tools/FastViromeExplorer/ncbi-viruses-list.txt; 
    	    done
    	fi
    done
}

function main() {
    mkdir results2/;
    conda activate enrichseq;
    test_enrichseq;
    conda activate FastViromeExplorer;
    test_viromexplorer;
}
main
