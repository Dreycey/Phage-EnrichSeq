##Creating BlastDB

### Creating multifasta input for blastdb creation (called PHAGEGENOMES.fa)
```
for file in phageDBgenomes/*; do echo $
file; cat $file >> PHAGEGENOMES.fa; done; for file1 in phage_genomes/*; do echo $file1; cat $file1 >> PHAGEGENOMES.fa; done
```
2. Counting number of genomes. (6589)
```
cat PHAGEGENOMES.fa | grep ">" | wc -l
```

3. Making the blast database.
```
makeblastdb -in PHAGEGENOMES.fa -dbtype nucl 
```
OUTPUT
```
Adding sequences from FASTA; added 6588 sequences in 3.8471 seconds.
```

4. Using BLASTN on simulated reads
```
blastn -db PHAGEGENOMES.fa -query ../readsimulator_module/sim_illumina.fa -out blast_out -outfmt "6 qseqid sseqid sblastnames score evalue pident length mismatch gap open gaps sseq"
```

OUTPUT example:
```
NC_023720.1-5578    Mycobacterium   N/A 88  1.81e-39    96.000  100 4   0   0   AGCGGTCTGGCCGCCGAACCCGGCGGCACCGGTGAGGATGTACAGCCCCTGCGGAGCCTTCGTGACGCGGCCGTACTTCTCCGGGTCGAAGACACCGGAC
NC_023720.1-5578    Mycobacterium   N/A 88  1.81e-39    96.000  100 4   0   0   AGCGGTCTGGCCGCCGAACCCGGCGGCACCGGTGAGGATGTACAGCCCCTGCGGAGCCTTCGTGACGCGGCCGTACTTCTCCGGGTCGAAGACACCGGAC 
```

## SCRIPTS
Python script used for making sure names work:

```
#! usr/bin/python3                                                              
import sys                                                                      
                                                                                
"""                                                                             
This file converts the names for the phage genomes                              
into identifiers that work better for the blast DB.                             
                                                                                
USAGE:                                                                          
python convertnames.py PHAGEGENOMES.fa PHAGEGENOMES_renamed.fa                  
"""                                                                             
# parameters in                                                                 
infile = open(sys.argv[1]).readlines()                                          
outfile = open(sys.argv[2], "w")                                                
                                                                                
# loop through file                                                             
index = 0                                                                       
for line in infile:                                                             
    if line.startswith(">"):                                                    
        if line.startswith(">NC"):                                              
            genename = line.strip(">").strip(",").split(" ")[1:4]               
            genename = ''.join(genename)                                        
        else:                                                                   
            genename = line.strip(">").strip(",").split(" ")[0:3]               
            genename = '_'.join(genename)                                       
        outfile.write(">" + genename)                                           
        outfile.write("\n")                                                     
        index += 1                                                              
    else:                                                                       
        outfile.write(line)                                                     
outfile.close() 
New BlastDB creation:
makeblastdb -in PHAGEGENOMES_renamed.fa
 -dbtype nucl
BLASTN:
blastn -db PHAGEGENOMES_renamed.fa -query ../readsimulator_module/sim_illumina.fa -out blast_out -outfmt "6 qseqid sseqid sblastnames score evalue pident length mismatch gap open gaps sseq"
Output:
```
