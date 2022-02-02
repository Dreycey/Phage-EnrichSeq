# Output Consolidator Module
This module is used to collect and organize the important information from the pipeline; namely, phage classifcations and relative abundance.

## Plotting in ```plotter.py```
Currently contains a method to plot a pie chart given labels and abundances. Labels can be taxids or cluster ids. 

Users can define their own plotter functions for other types of graphs.

## File reorganization in ```consolidator.py```
Contains methods for file reorganization. When running ```enrichseq.nf```, each module produces its own output file(s) in their respective working directories. This script provides functions to move or copy important results to one location.

## Wrapper script ```client.py```
This script is the one called by ```enrichseq.nf``` in the ```Run_CombineOutput``` process. This script's input and output requirements should not change. 