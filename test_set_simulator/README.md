# Test Set Simulator
This script is used to make test datasets for the goal at hand: 
using sequencing to simulated genome amplification by engrichement. 
This allows for validating a computational pipeline orientated at the task
before actually testing the task in the laboratory.

# USAGE
command:

```
python testsetmaker.py <path to config> <outfile path/name>
```

example:

```
python testsetmaker.py ./config_example.config output_dir/simulated_data.fa
```

## General Outline
The goal of this script is take input as a config file that
tells the script several things:

1. Input Reads Name
2. The percentage used for output file
   -- add stochastic option, sampling from a dist?
3. Path to the simulated reads.

## config file

* general set up
```
<name>, <percentage>, <path>
```

* example
```
phage_1 80 ~/dreycey_albin/Desktop/simulated_reads_directory/phage_1_illumina.fa
phage_2 10 ~/dreycey_albin/Desktop/simulated_reads_directory/phage_2_illumina.fa
phage_3 5 ~/dreycey_albin/Desktop/simulated_reads_directory/phage_3_illumina.fa
phage_4 5 ~/dreycey_albin/Desktop/simulated_reads_directory/phage_4_illumina.fa
```

## future updates

(1) Make it so certain regions of the genome are sampled
    -- or add this as an option
    -- Does this happen in reality (check)
(2) all for 
