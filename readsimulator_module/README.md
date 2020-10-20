# ReadSimulator Module
Read simulation will take in a reference genome as input and simulate reads for 
Illumina, PacBio and Nanopore. 


## Usage

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
