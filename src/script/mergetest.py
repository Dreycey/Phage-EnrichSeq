from src.modules.PathOrganizer_module.PathOrganizer import PathOrganizer as PathO, DuplicateGenomeError, InvalidQueryError
import sys

genome_db = "/Users/dreyceyalbin/Desktop/Phage-EnrichSeq/database/ref_genomes/"
genome_to_grab = sys.argv[1]
path_bj = PathO(genome_db)

print(path_bj.genome(genome_to_grab))
