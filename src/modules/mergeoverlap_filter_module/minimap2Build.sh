#!/bin/bash
#
# Perform the build for the Blast database
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



#######################################
# downloadMappy() 
# creates a multifasta from a directory of fasta files.
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Clones the github repo for the minimap2 wrapper
#######################################
function downloadMappy() {
    echo "Downloading Minimap2 read simulation testing..";
    if [ ! -d "${DIR}/minimap2" ]
    then
        git clone https://github.com/lh3/minimap2 ${DIR}/minimap2;
        cd ${DIR}/minimap2;
        python setup.py install;
    fi    
}


#######################################
# Main function runs all of the other methods for the build
# Globals:
#   None
# Arguments:
#   None
# Outputs:
#   Downloads the mappy wrapper for python
#######################################
function main(){
    # download mappy
    downloadMappy;
}

echo "Running the minimap2 build script";
main;
