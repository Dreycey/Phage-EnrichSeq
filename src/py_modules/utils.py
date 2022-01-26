"""
This module consists of utility methods 
that can be accessed by the rest of the pipeline
"""
import subprocess
from typing import List


def subproc_call(CMD: List, shell=False, output=False):
    """
    This method uses os subsystem calls to 
    start processes for input commands and arguments.

    INPUT:
        1. command line argument top start
    OUTPUT:
        1. Output from the command
    """
    output = subprocess.run(CMD, shell=shell, capture_output=output)
    return output