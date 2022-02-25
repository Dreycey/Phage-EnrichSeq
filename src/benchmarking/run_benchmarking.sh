

function test_enrichseq() {
    tool_name="enrichseq";
    file_suffix="fa";
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
                           -o results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}};
    	    done
    	fi
    done
}

function test_viromexplorer() {
    tool_name="FastViromeExplorer";
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
                                                        -o results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}} \
                                                        -l tools/FastViromeExplorer/ncbi-viruses-list.txt; 
    	    done
    	fi
    done
}

function test_kraken2() {
    tool_name="Kraken2";
    file_suffix=".fa";
    echo "TESTING ${tool_name}";
    for test_dir in ./tests/*; do
        if [[ -d $test_dir ]]; then
            for file in $test_dir/*${file_suffix}; do
               echo; echo $file; echo;
               basefile="$(basename -- $file)";
               test_dir_name="$(basename "${test_dir##*/}")";
               mkdir -p results/${tool_name}/${test_dir_name}/;
               kraken2 --use-names --threads 4 --db ${dbDir} --report ${reportFile} ${inFasta} > ${krakenOutFile}
           	   #java -cp tools/FastViromeExplorer/bin FastViromeExplorer \
               #                                         -1 $file \
               #                                         -i tools/FastViromeExplorer/ncbi-virus-kallisto-index-k31.idx \
               #                                         -o results/${tool_name}/${test_dir_name}/${basefile%${file_suffix}} \
               #                                         -l tools/FastViromeExplorer/ncbi-viruses-list.txt; 
    	    done
    	fi
}

function main() {
    mkdir results/;
    conda activate enrichseq;
    #test_enrichseq;
    test_kraken2;

    #conda activate FastViromeExplorer;
    #test_viromexplorer;
}
main
