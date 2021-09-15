# ReadSimulator Module
Read simulation will take in a reference genome as input and simulate reads for 
Illumina, PacBio and Nanopore. 

> The EnrichSeq Readsimulator aims to simulate both short and long reads at different amounts based on the input genomes and their relative abundance levels. 

## Quick Usage

* install simulators
```
bash install.sh mac
```

* usage
```
usage: simulate_reads.py [-h] [-v | -q] -i CONFIG -c COVERAGE -o OUTPUT -t THREADS
```

* example
```
python simulate_reads.py -i simulate_genomes.config -rn 50000 -c 30 -o simulatedgenomes -t 4
```

***
***

## Detailed Usage

### Testing the read simulator
Pytest is used as the testing library. Both unit testin and regression testing is set up for the module. 

* To run the suite of tests, use the following command:
```
pytest --cov
```

### install dependencies

python3: Make sure to install the HTSeq python library:
```
pip3 install HTSeq
```

python3: the joblib module:
```
pip3 install joblib
```

python3: Scipy
```
pip3 install scipy
```

python3: scikit-learn
```
pip3 install -U scikit-learn
```

### install script
* Download software and create the directory structure

Mac
```
bash install.sh mac
```

Linux
```                                                                             
bash install.sh linux                                                          
```

Windows Linux Subsystem or Cygwin
```
sudo apt-get install dos2unix;
dos2unix install.sh;
dos2unix simulate_reads.sh;
bash install.sh linux;
```

### simulating reads from input reference genome
**NOTE**
When specifying different number of reads in the options, make sure to delete all \*\_simulatedreads directories:
```
rm -rf *simulatedreads
```
USAGE
```
bash simulate_reads.sh <reference genome PATH> <Outfile name> <Outfile path> <Threads>
```

EXAMPLE (using SARS-CoV-2)
```
bash simulate_reads.sh ./sarscov2_reference.fa test_output_prefix outfiles_directory/ 4
```

EXAMPLE (using an M. Smeg Phage, D29 )
```
bash simulate_reads.sh ./NC_001900.fna D29reads d29simulatedReads_directory/ 4
```
