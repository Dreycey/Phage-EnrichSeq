######
# DESCRIPTION:
#    This runs a read simulation script on
#    an entire directory of config files.
#####


num_reads_arr=( 100000 200000 300000 400000 500000 600000 700000 800000 900000 1000000 );
test=3;
num_genomes=800;
for read_count in "${num_reads_arr[@]}"; do
    python simulate_reads.py -i ENRICHSEQ_BENCHMARKING/test_${test}/${num_genomes}_genomes_${test}.config \
                             -o ${num_genomes}_genomes_${read_count}_reads.fa \
                             -rn ${read_count} \
                             -t 4;
done


