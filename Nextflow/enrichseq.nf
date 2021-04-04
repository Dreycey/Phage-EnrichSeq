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

BASE = fastafile.getName()
THREADS = params.threads
WORKFLOW = "enrichseq"

megahitDir = file("$workingDir/$WORKFLOW/megahit")
krakenDir = file("$workingDir/$WORKFLOW/kraken")
brackenDir = file("$workingDir/$WORKFLOW/bracken")
blastWorkingDir = file("$workingDir/$WORKFLOW/blast")


process Create_Working_Directories {
    output:
    stdout create

    """
    mkdir -p $workingDir
    mkdir -p $workingDir/$WORKFLOW

    if [ -d $megahitDir ]; then rm -rf $megahitDir; fi;
    if [ -d $krakenDir ]; then rm -rf $krakenDir; fi;
    if [ -d $brackenDir ]; then rm -rf $brackenDir; fi;
    if [ -d $blastWorkingDir ]; then rm -rf $blastWorkingDir; fi;

    #mkdir $megahitDir
    mkdir $krakenDir
    mkdir $brackenDir
    mkdir $blastWorkingDir
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
        		  --input=${fastafile} \
    			  --threads=${THREADS} \
    			  --out=${megahitDir}
    """
}


process Prep_Databases {
	input:
	val init from init

	output:
	stdout databases

	script:
	"""
	if [[ ! -f ${databasesDir}/taxo.k2d ]]; then
		bash ${params.toolpath}/kraken_module/kraken2Build.sh
	else
		echo "Kraken DB already exists" > ${krakenDir}/kraken.log
	fi

	"""
}


process Run_Kraken {
	input:
	val megaout from megahit
	val krakendb from databases

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


process Run_Bracken {
	input:
	val krakenout from kraken

	output:
	stdout bracken

	script:
	"""
	echo "Running Run_Bracken" > ${brackenDir}/bracken.log
	bash ${params.toolpath}/bracken_module/brackenBuild.sh
	bash ${params.toolpath}/bracken_module/brackenRun.sh --krakendb=${databasesDir} \
				--input=${krakenDir}/kraken_orig.report \
				--out=${brackenDir}/bracken_run_orig \
				--read=${params.readlength}
	"""
}

process Run_BLAST {
        input:
        val blastdb from databases
	val megaout from megahit         

        output:
        stdout blast

        script:
        """
        if [[ ! -f ${params.toolpath}/blast_module/blastdb/outputMulti3.fa ]]; then
                bash ${params.toolpath}/blast_module/blastBuild.sh
        else
                echo "Blast DB already exists" > ${blastWorkingDir}/blast.log
        fi
        bash ${params.toolpath}/blast_module/blastRun.sh --blastdb=${params.toolpath}/blast_module/blastdb/outputMulti3.fa \
				 --queryfasta=${fastafile} \
				 --out=${blastWorkingDir}/blastout.txt
        """
}

create.subscribe { print "$it" }
init.subscribe { print "$it" }
megahit.subscribe { print "$it" }
//databases.subscribe { print "$it"}
kraken.subscribe { print "$it" }

bracken.subscribe { print "$it" }
blast.subscribe { print "$it" }
