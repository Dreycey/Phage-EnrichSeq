# Building the minimap2 module
## Installing Mappy

* Conda
```
conda install -c bioconda mappy
```
 *NOTE*: Should be included in the yml file

## Manually
* Run the following
```
bash minimap2Build.sh
```

# Running the tool
Usage:
```
python minimap2_module/minimap2wrapper.py <multi fasta with genomes> <genomes to look at> <input sequence reads> <out path name>
```

Example:
```
python minimap2_module/minimap2wrapper.py phageMulti.fa genomes_to_view.txt simulatedgenomes_illumina.fasta outname

```

where genomes_to_view.txt:
```
Ryadel  
Blessica     
D29      
Paphu     
Perseus  
```
