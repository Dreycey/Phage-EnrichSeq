#!/usr/bin/env nextflow

Executor = 'local'

def usage() {
     log.info ''
	log.info 'Usage: nextflow run enrichseq.nf --read single --fasta /Path/to/infile.fasta --workdir /Path/to/working_directory --dbdir /Path/to/databases [--threads 4] [--log=/Path/to/run.log]'
	log.read '  --read       Ready type (single / paired / long)'
	log.info '  --fasta		Path to the input FASTA file'
	log.info '  --workdir	Path to the output working directory'
	log.info "  --dbdir		Path to the classification databases"
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

logfile = file(params.log)
fastafile = file(params.fasta)
workingDir = file(params.workdir)
databasesDir = file(params.dbdir)
// pfamDir = file("$workingDir/seqmapper/pfam_hmm")
// create_pfam = params.hmmscan

BASE = fastafile.getName()
MODULES =  "$workflow.projectDir/../modules"
SCRIPTS =  "$workflow.projectDir/../scripts"
THREADS = params.threads
WORKFLOW = "enrichseq"

megahitDir = file("$workingDir/$WORKFLOW/megahit")
krakenDir = file("$workingDir/$WORKFLOW/kraken")
brackenDir = file("$workingDir/$WORKFLOW/bracken")

translatedFile = "${workingDir}/initialize/six_frame_translation/${BASE}.translated.fasta"


process Create_Working_Directories {
    output:
    stdout create

    """
    mkdir -p $workingDir
    mkdir -p $workingDir/$WORKFLOW

    if [ -d $megahitDir ]; then rm -rf $megahitDir; fi;
    if [ -d $krakenDir ]; then rm -rf $krakenDir; fi;
    if [ -d $brackenDir ]; then rm -rf $brackenDir; fi;

    #mkdir $megahitDir
    mkdir $krakenDir
    mkdir $brackenDir

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
    //path(megahit_out) into megahit
    stdout megahit

    if ( Executor == 'local' ) {
       executor "local"
    }
  /*
    else if ( Executor == 'slurm' ) {
       clusterOptions "--ntasks-per-node $THREADS"
       executor "slurm"
    }
  */
    script:
    //megahit_out = "$megahitDir/*.contigs.fa"
    """
    bash $params.megahitpath/megahitRun.sh --read=$params.read \
        		  --input=${fastafile} \
    			  --threads=${THREADS} \
    			  --out=${megahitDir}
    """
}

/*
process Prep_Databases {
	input:
	val init from init

	output:
	stdout databases

	script:
	"""
	if [ ! -f $databasesDir/taxo.k2d ]; then
		bash $params.krakenpath/kraken2Build.sh
	else
		echo "Kraken DB already exists" | tee -a $logfile
	fi

	"""
}
*/


process Run_Kraken {
	input:
	//path assembled_fasta from megahit
	val megaout from megahit

	output:
	stdout kraken

	script:
	"""
	bash $params.krakenpath/kraken2Run.sh --krakendb=${databasesDir} \
				--queryfasta=$megahitDir/*.contigs.fa \
				--report=$krakenDir/kraken.report \
				--out=$krakenDir/kraken.log
	"""
}


process Run_Bracken {
	input:
	val krakenout from kraken

	output:
	stdout bracken

	script:
	"""
	bash $params.brackenpath/brackenBuild.sh
	bash $params.brackenpath/brackenRun.sh --krakendb=${databasesDir} \
				--input=$krakenDir/krakenout.report \
				--out=bracken_run1 \
				--read=$params.readlength
	"""
}

create.subscribe { print "$it" }
init.subscribe { print "$it" }
megahit.subscribe { print "$it" }
//databases.subscribe { print "$it"}
kraken.subscribe { print "$it" }
bracken.subscribe { print "$it" }
