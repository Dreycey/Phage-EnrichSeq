import sys
config_file = open(sys.argv[1]).readlines()
max_count = int(sys.argv[2])
for count, line in enumerate(config_file):
    if count < max_count:
        path = '/Users/dreyceyalbin/Desktop/Phage-EnrichSeq/database/ref_genomes/'
        path += line.split(" ")[0].split("/")[-1].replace(",", "")
        abundance = 1/max_count
        print(f'{path}, {abundance}')
    else:
        break
