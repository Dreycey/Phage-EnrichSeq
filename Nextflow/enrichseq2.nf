#!/usr/bin/env nextflow

megahitDir = "megahit_output"
assembly_ch = Channel.fromPath(params.query)

// change the if statement to a 'when' clause before the script
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

process runKraken2 {
	input:
	path assembled_fasta from classify_ch

	script:
	"""
	echo "$assembled_fasta" > $PWD/testing.txt
	"""
}
