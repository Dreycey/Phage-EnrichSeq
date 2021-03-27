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
