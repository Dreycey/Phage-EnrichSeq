import sys
infile = open(sys.argv[1])
outfile = []

index = 0
for line in infile:
    if line.startswith(">"):
        if (outfile != []): outfile.close()
        genename = line
        filename = line.split(" ")[2] + ".fa" #str(index) + ".fa" #genename.replace(" ", "") + ".fa"
        outfile = open(filename,'w')
        outfile.write(line)
        index += 1
    else:
        outfile.write(line)
outfile.close()
