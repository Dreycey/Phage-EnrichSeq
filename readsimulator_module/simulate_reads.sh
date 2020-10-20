
# print usage if not enough input!
usage() {
    echo; echo "Usage: bash simulate_reads.sh <reference genome PATH> <Outfilepath> <Threads>";
    echo; echo "Example: bash simulate_reads.sh ./sarscov2_reference.fa test_output_prefix outfiles_directory/ 4";
    exit 1;
}

# getting the global variables set correctly
ref_gen=$1;
out_name=$2;
out_directory=$3;
THREADS=$4;

if [[ ! -d $out_directory ]]
then
    mkdir $out_directory;
fi

# function for simulating Illumina reads
function simulate_illumina {
    echo "Simulating short illumina reads using ART..";
    if [ ! -f ${out_directory}/${out_name}_illumina.fq ] 
    then
        ART_CMD_ARGS="-ss MSv3 -sam -i ${ref_gen} -l 100 -f 30 -o ${out_name}_illumina";
        command ./art_bin_MountRainier/art_illumina $ART_CMD_ARGS;
        mv ${out_name}_illumina* ${out_directory};
    fi
}

function simulate_nanopore {
    echo "Simulating NanoPore reads..";
    if [ ! -f ${out_directory}/${out_name}_aligned_reads.fasta ]
    then
        TRAINING="human_NA12878_DNA_FAB49712_albacore/training"
        NANO_CMD_ARGS="-rg sarscov2_reference.fa -c ${TRAINING} -o ${out_name} -t ${THREADS}"
        python3 NanoSim/src/simulator.py genome $NANO_CMD_ARGS;
        mv ${out_name}* ${out_directory};
    fi
}

function simulate_pacbio {
    echo "Simulating Pacbio reads..";
    if [ ! -f ${out_directory}/${out_name}_pacbio.fq ]
    then
        # create the index
        perl PaSS/pacbio_mkindex.pl ${ref_gen} ${out_directory};
        # use these files to simulate the pacbio reads
        ./PaSS/PaSS -list ${out_directory}/percentage.txt -index ${out_directory}/index -m pacbio_RS  -c PaSS/sim.config -r 2000 -t ${THREADS} -o ${out_name}_pacbio -d; 
        mv ${out_name}_pacbio* ${out_directory};
    fi
}

# main function for controlling the workflow
function main {
    simulate_illumina;
    simulate_nanopore;
    simulate_pacbio;
}
main;
