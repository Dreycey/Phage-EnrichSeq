#!/usr/bin/env nextflow

Executor = 'local'

def usage() {
    log.info ''
	log.info 'Usage: nextflow run enrichseq.nf --read single --fasta /Path/to/infile.fasta --workdir /Path/to/working_directory --dbdir /Path/to/databases [--threads 4] [--log=/Path/to/run.log]'
	log.read '  --read       Ready type (single / paired / long)'
	log.info '  --fasta		Path to the input FASTA file. If paired end reads, only use file prefix (i.e. before underscore)! Example: '--fasta fastaA' if 'fastaA_1.fa'/'fastaA_2.fa' [suffix must be "_{1/2}.fa"])'
	log.info '  --workdir	Path to the output working directory'
	log.info "  --dbdir		Path to the classification databases"
    log.info "  --genomedir Path to the genome directory, built running the Krake2Build.sh script"
	log.info "  --threads	Number of threads to use (Default=1)"
	log.info "  --log		Log file (Default=Don't save the log)"
	log.info '  --help		Print this help message out'
	log.info ''
	exit 1
}

if (params.help) {
    usage()
}


if (!params.fasta) {
	log.error "Missing argument for --fasta option"
	usage()
}

if (!params.workdir) {
    log.error "Missing argument for --workdir option"
    usage()
}

if (!params.dbdir) {
    log.error "Missing argument for --dbdir option"
    usage()
}

// Get reads ready based on type
if (params.read == 'single') {
    println "Single end reads entered!"
    fastafile = file(params.fasta)
    fastafile_2 = file(params.fasta) // UNUSED IN THIS CASE.
}
else if (params.read == 'paired') {
    println "Paired end reads entered!"
    fastafile = file(params.fasta+'_1.fa')
    fastafile_2 = file(params.fasta+'_2.fa')
}

// opening files for misc parameters
logfile = file(params.log)
workingDir = file(params.workdir)
databasesDir = file(params.dbdir)
genomeDir = file(params.genomedir)

BASE = fastafile.getName()
THREADS = params.threads
WORKFLOW = "enrichseq"

megahitDir = file("$workingDir/$WORKFLOW/megahit")
krakenDir = file("$workingDir/$WORKFLOW/kraken")
mergeOverlapDir = file("$workingDir/$WORKFLOW/merge_overlap_filter")
genomeCompareDir = file("$workingDir/$WORKFLOW/genome_comparison")
outputDir = file("$workingDir/$WORKFLOW/output_files")


process Create_Working_Directories {
    output:
    stdout create

    """
    mkdir -p $workingDir
    mkdir -p $workingDir/$WORKFLOW

    if [ -d $megahitDir ]; then rm -rf $megahitDir; fi;
    if [ -d $krakenDir ]; then rm -rf $krakenDir; fi;
    if [ -d $mergeOverlapDir ]; then rm -rf $mergeOverlapDir; fi;
    if [ -d $genomeCompareDir ]; then rm -rf $genomeCompareDir; fi;
    if [ -d $outputDir ]; then rm -rf $outputDir; fi;

    #mkdir $megahitDir
    mkdir $krakenDir
    mkdir $mergeOverlapDir
    mkdir $genomeCompareDir
    mkdir $outputDir
    """
}

process Initialize {
    input:
    val create from create

    output:
    stdout init

    executor Executor

    """
    echo -n " # Launching $WORKFLOW workflow ........................ " | tee -a $logfile; date '+%H:%M:%S %Y-%m-%d' | tee -a $logfile
    """
}

process Run_Megahit {
    input:
    val init from init

    output:
    stdout megahit

    if ( Executor == 'local' ) {
       executor "local"
    }

    script:
    """
    bash ${params.toolpath}/megahit_module/megahitRun.sh --read=${params.read} \
        		  --input1=${fastafile} \
                  --input2=${fastafile_2} \
    			  --threads=${THREADS} \
    			  --out=${megahitDir}
    """
}

process Run_Kraken {
	input:
	val megaout from megahit
	//val krakendb from databases

	output:
	stdout kraken
	
	script:
	"""
	bash ${params.toolpath}/kraken_module/kraken2Run.sh --krakendb=${databasesDir} \
				--queryfasta=${megahitDir}/*.contigs.fa \
				--report=${krakenDir}/kraken_assembled.report \
				--out=${krakenDir}/kraken1.log
	bash ${params.toolpath}/kraken_module/kraken2Run.sh --krakendb=${databasesDir} \
				--queryfasta=${fastafile} \
				--report=${krakenDir}/kraken_orig.report \
				--out=${krakenDir}/kraken2.log
	"""
}


process Run_KrakenParser {
   input: 
   val krakenout from kraken

   output:
   stdout kraken_parser

   script:
   """
   python ${params.toolpath}/kraken_module/parseKraken.py ${krakenDir}/kraken_assembled.report \
          ${krakenDir}/taxid_file.txt ${krakenDir}/parsed_kraken_phages.txt
   """
}


process Run_MergeOverlap {
	input:
    val kraken_parse from kraken_parser

	output:
	stdout merge_overlap

	script:    
	"""
	echo "Running the Merge Overlap Filter" > ${mergeOverlapDir}/mergeOverlap.log;
	python ${params.toolpath}/mergeoverlap_filter_module/mergeoverlap.py \
            --input ${krakenDir}/taxid_file.txt \
            --output_prefix ${mergeOverlapDir}/merge_overlap_out \
            --genome_directory ${genomeDir} \
            --fasta ${fastafile} \
            --threads ${THREADS} 
	"""
}

process Run_GenomeComparison {
	input:
    val filter from merge_overlap

	output:
	stdout clusters

	script:    
	"""
	echo "Running the Genome Comparison module" > ${genomeCompareDir}/genomeCompare.log;
	python ${params.toolpath}/genomeCompare_module/genome_comparison.py \
            --input ${mergeOverlapDir}/ \
            --output_dir ${genomeCompareDir}/ \
            --genome_directory ${genomeDir} \
            --kmer_length ${params.kmer_length} \
            --threshold ${params.clustering_threshold}
	"""
}

process Run_CombineOutput {
    input:
    val cluster_out from clusters

    output:
    stdout enrichseq_out

    script:
    """
    echo "Consolidating output" > ${outputDir}/genomeCompare.log;
    python ${params.toolpath}/output_module/client.py \
            --inputdir ${workingDir}/${WORKFLOW} \
            --outputdir ${outputDir}
    """
}

create.subscribe { print "$it" }
init.subscribe { print "$it" }
megahit.subscribe { print "$it" }
kraken.subscribe { print "$it" }
merge_overlap.subscribe { print "$it" }
clusters.subscribe { print "$it" } 
enrichseq_out.subscribe { print "$it"}
