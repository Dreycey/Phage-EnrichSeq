
"""
DESCRIPTION:
    This module contains objects and methods for path disambiguation. 
"""
# std packages
import os
from pathlib import Path
from typing import Optional
import re
# non-std packages




# error handeling
class PathErrors(Exception):
    pass

class DuplicateGenomeError(PathErrors):
    """
    Exception raised if there are more than one genome 
    with a particular name.
    """

    def __init__(self, genome_name, message):
        self.genome_name = genome_name
        self.message = message

class HeaderError(PathErrors):
    """
    Exception raised if the file header is 
    not correct.
    """

    def __init__(self, genome_name, message):
        self.genome_name = genome_name
        self.message = message

# path storing class
class PathOrganizer:
    """ This data structure holds paths and retrieves information """

    def __init__(self, database):
        self.database_path: str = str(database)

    # @property
    def genome(self, genome_taxid) -> Optional[Path]:
        """ This getter grabs a genome file based on the file name prefix """
        genome_path = None
        for genome_in_db in os.listdir(self.database_path):
            full_path = Path(self.database_path) / Path(genome_in_db)
            if (int(self.get_fasta_taxid(full_path)) == int(genome_taxid)):
                if (genome_path == None):
                    genome_path = full_path
                else:
                    raise DuplicateGenomeError(genome_path, f"Duplicate genome for {full_path}, rerun the DB Build script!!")
        return genome_path

    def get_fasta_taxid(self, genome_in_db):
        """
        parse an input fasta for the taxid.
        """
        with open(genome_in_db) as fasta_file:
            try:
                taxid = fasta_file.readline().split("|kraken:taxid|")[1].split(" ")[0]
            except:
                taxid = 0 
        return taxid