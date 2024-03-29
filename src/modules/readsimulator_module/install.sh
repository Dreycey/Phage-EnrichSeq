
# setting up global vars
os_type=$1;

if  [[ $os_type == "" ]]; then
    echo; echo;
    echo "USAGE:";
    echo "bash install mac";
    echo "or";
    echo "bash install linux"; echo; echo;
    exit;
fi

# downlaod ART - for Illumina
function download_art {
    echo "Downloading ART for short read simulation..";

    if [ ! -d  "art_bin_MountRainier" ]
    then
        # commands for linux 
        if [[ $os_type == "linux" ]]
        then
            wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
            tar -zxvf artbinmountrainier2016.06.05linux64.tgz;
            rm artbinmountrainier2016.06.05linux64.tgz;    
        fi
        # commands for mac
        if [[ $os_type == "mac" ]]
        then
            wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05macos64.tgz;
            tar -zxvf artbinmountrainier2016.06.05macos64.tgz;
            rm artbinmountrainier2016.06.05macos64.tgz
        fi
    fi
}

# download NanoSim - for Nanopore
function download_nanosim {
    echo "Downloading NanoSim for Nanopore read simulation..";

    if [ ! -d "NanoSim" ]
    then
        git clone https://github.com/bcgsc/NanoSim.git;
        tar -zxvf NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_albacore.tar.gz;
    fi
}

# download PaSS - for PacBio
function download_pass {                                                     
    echo "Downloading PaSS for PacBio read simulation..";                  
    if [ ! -d "PaSS" ]
    then
        wget http://cgm.sjtu.edu.cn/PaSS/src/PaSS.tar.gz;
        tar -zxvf PaSS.tar.gz;
        rm PaSS.tar.gz;
        cd PaSS; gcc -lm -lpthread PaSS.c -o PaSS; cd ../;
    fi
} 

# download minimap2 for read mapping
function download_minimap2 {
    echo "Downloading Minimap2 read simulation testing..";
    if [ ! -d "minimap2" ]
    then
        git clone https://github.com/lh3/minimap2;
        cd minimap2;
        python setup.py install;
        cd ../;
    fi
}

function main {
    download_art;
    download_nanosim;
    download_pass;
    download_minimap2;
}
main;
