"""
DESCRIPTION:
    This script contains the unit tests for the PathOrganizer class. The
    major functions are tested here, and these rely on the testing/directory.
"""
# std packages
# non-std packages
import pytest
from pathlib import Path
# in-house package
from PathOrganizer import PathOrganizer as PathO, DuplicateGenomeError



# SETTING UP FIXTURES
@pytest.fixture
def genome_testing_path():
    """ contains major paths """
    #return "../kraken_module/ref_genomes"
    return "testing/"

@pytest.fixture
def object_PathOrganizer(genome_testing_path):
    """ contains major paths """
    return PathO(genome_testing_path,"exampleoutdir")

# TESTING
@pytest.mark.parametrize("genome_to_grab, expected", [
    ("a", None),
    ("non_existant", None),
    ("GCF_002997835.1_ASM299783v1_genomic.fna",Path("testing/GCF_002997835.1_ASM299783v1_genomic.fna")),
    ("Streptomyces_phage_Yara|kraken:taxid|1235691 complete sequence, 68671 bp, circularly permuted, Cluster BN.fna",
     Path("testing/Streptomyces_phage_Yara|kraken:taxid|1235691 complete sequence, 68671 bp, circularly permuted, Cluster BN.fna"))
     ])
def test_genome(object_PathOrganizer, genome_to_grab, expected):
    """
    Tests the function for grabbing genomes.
    """
    genome_output_path = object_PathOrganizer.genome(genome_to_grab)
    assert genome_output_path == expected

