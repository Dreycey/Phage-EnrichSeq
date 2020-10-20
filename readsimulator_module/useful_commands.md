## Individual Tool Commands                                                     
                                                                                
### ART - Illumina simulator                                                    
                                                                                
#### running art with MiSeq version 3 simulator                                 
```                                                                             
./art_bin_MountRainier/art_illumina  -ss MSv3 -sam -i ./sarscov2_reference.fa -l
100 -f 30 -o art_outfile
```                                                                             
                                                                                
#### output files                                                               
```                                                                             
art_outfile.aln art_outfile.fq art_outfile.sam                                  
```                                                                             
                                                                                
### PaSS - PacBio simulator                                                     
#### preparing index using input genome                                         
```                                                                             
perl PaSS/pacbio_mkindex.pl ./sarscov2_reference.fa ./                          
```                                                                             
                                                                                
#### using PaSS                                                                 
```                                                                             
./PaSS/PaSS -list percentage.txt -index index -m pacbio_RS -c PaSS/sim.config -r
2000 -t 3 -o output_name -d
```                                                                             
                                                                                
#### output files                                                               
```                                                                             
output_name.fq                                                                  
```                                                                             
                                                                                
### NANOSIM - Nanopore simulator                                                
#### unpacking the training model                                               
```                                                                             
tar -zxvf NanoSim/pre-trained_models/human_NA12878_DNA_FAB49712_albacore.tar.gz 
```                                                                             
                                                                                
#### running genome simulator                                                   
```                                                                             
python NanoSim/src/simulator.py genome -rg sarscov2_reference.fa -c
human_NA12878_DNA_FAB49712_albacore/training
```                                                                             
                                                                                
#### find number of reads                                                       
```                                                                             
cat simulated_aligned_reads.fasta | grep ">" | wc -l                            
```                                                                             
                                                                                
#### output files                                                               
```                                                                             
simulated_aligned_error_profile                                                 
simulated_aligned_reads.fasta                                                   
simulated_unaligned_error_profile                                               
simulated_unaligned_reads.fasta                                                 
```                                                                             
