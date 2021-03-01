


blastn -db PHAGEGENOMES_renamed.fa -query ../readsimulator_module/sim_illumina.fa \
-out blast_out -outfmt "6 qseqid sseqid sblastnames score evalue pident length mismatch gap open gaps sseq"
