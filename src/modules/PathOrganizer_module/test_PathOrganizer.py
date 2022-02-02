"""
DESCRIPTION:
    This script contains the unit tests for the PathOrganizer class. The
    major functions are tested here, and these rely on the testing/directory.
"""
# std packages
# non-std packages
import pytest
import os
from pathlib import Path
# in-house package
from PathOrganizer import PathOrganizer as PathO, DuplicateGenomeError, InvalidQueryError


CURR_PATH = os.path.realpath(__file__)
CURR_PATH = os.path.dirname(CURR_PATH)

@pytest.fixture
def genome_testing_path():
    """ contains major paths """
    return f"{CURR_PATH}/testing_files/"

@pytest.fixture
def object_PathOrganizer(genome_testing_path):
    """ contains major paths """
    return PathO(genome_testing_path)

# TESTING
@pytest.mark.parametrize("genome_to_grab, expected", [
                                                        ("1636581", "GCF_002593425.1_ASM259342v1_genomic.fna"),
                                                        ("2099652", "GCF_002997835.1_ASM299783v1_genomic.fna")
                                                     ]
                        )
def test_genome(object_PathOrganizer, genome_to_grab, expected):
    """
    Tests the function for grabbing genomes.
    """
    genome_output_path = object_PathOrganizer.genome(genome_to_grab)
    assert genome_output_path.name == expected

@pytest.mark.parametrize("genome_to_grab, expected_error", [
                                                        ("a",  InvalidQueryError),
                                                        ("a1k2nmd",  InvalidQueryError),
                                                        ("KMFKM",  InvalidQueryError),
                                                        ("22O",  InvalidQueryError),
                                                        ("55884", DuplicateGenomeError)
                                                     ]
                        )
def test_errorHandeling(object_PathOrganizer, genome_to_grab, expected_error):
    """
    Tests the object for correct error handeling
    """
    with pytest.raises(expected_error) as error:
        object_PathOrganizer.genome(genome_to_grab)
