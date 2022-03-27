test_number=3;
config_file_arr=( 200_genomes_${test_number}.config \
                  400_genomes_${test_number}.config \
                  600_genomes_${test_number}.config );
for config_file in "${config_file_arr[@]}"; do
    python simulate_reads.py -i ENRICHSEQ_BENCHMARKING/test_${test_number}/${config_file} \
                             -o ${config_file%.config}_500000_reads.fa \
                             -rn 500000 \
                             -t 4;
done

