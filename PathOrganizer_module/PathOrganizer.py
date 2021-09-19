
"""
DESCRIPTION:
    This module contains objects and methods for path disambiguation. 
"""
# std packages
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

# path storing class
class PathOrganizer:
    """ This data structure holds paths and retrieves information """

    def __init__(self, database, output_directory):
        self.database_path: Path = Path(database)
        self.output_path: Path = Path(output_directory)

    @property
    def genome(self, genome_file_prefix) -> Optional[Path]:
        """ This getter grabs a genome file based on the file name prefix """
        genome_path = None
        regex_expression = re.compile(f'({genome_file_prefix}.*fa$)|({genome_file_prefix}.*fna$)|({genome_file_prefix}.*fasta$)')
        for genome_path in self.database_path:
            if regex_expression.match(genome_path):
                if (genome_path == None):
                    genome_path = Path(genome_path)
                else:
                    raise DuplicateGenomeError(genome_path, "Duplicate genome for {genome_path}")
        return genome_path