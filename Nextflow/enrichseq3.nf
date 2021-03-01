#!/usr/bin/env nextflow

megahitDir = "megahit_output"
assembly_ch = Channel.fromPath(params.query)

// TODO: change the if statement to a 'when' clause before the script
process runMegahit {
	input:
	path input_fasta from assembly_ch

	output:
	path(megahit_out) into classify_ch

	script:
	megahit_out = "$megahitDir/*.contigs.fa"
	"""
	if [ ! -d "$megahitDir" ]
	then
		megahit -r $PWD/$params.query -t $params.threads -m 1e9 -o $megahitDir --out-prefix megahit
	fi
	"""
}

// this should only run if the databases have not already been set. how to check this?
process prepDatabases {
	script:
	"""
	if [ ! -d "$PWD/phage_genomes" ]
	then
		echo "running create_create_kraken_db.sh" > $PWD/testing.txt
	else
		echo "phage_genomes exists" > $PWD/testing.txt
	fi

	"""
}


process runKraken2 {
	input:
	path assembled_fasta from classify_ch

	script:
	"""
		kraken2 --use-names --threads 4 --db $params.dbDir --report $PWD/kr1.report.txt $assembled_fasta > $PWD/test1.kraken
		#echo $assembled_fasta > $PWD/testing.txt
	"""
}

process runMetaphlan {

}
