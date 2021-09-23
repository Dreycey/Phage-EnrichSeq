import pytest
import sys
import os
import parseKraken as parser
from typing import Dict
from pathlib import Path

''' Fixture setup '''
@pytest.fixture
def expected_phage_dict():
    return dict({ 2502430: "Mycobacterium phage Fushigi", 
            2599873: "Mycobacterium phage Paphu", 
            1913110: "Mycobacterium phage StarStuff"})

@pytest.fixture
def kraken_report_path():
    return Path(os.getcwd()) / Path("kraken_SAMPLE.report")

# @pytest.mark.parametrize("parsed_phages, expected", [
#     ({}, expected_phage_dict), # normal case
#     ("non_existant", None),
#     ("GCF_002997835.1_ASM299783v1_genomic.fna",Path("testing/GCF_002997835.1_ASM299783v1_genomic.fna")),
#     ("Streptomyces_phage_Yara|kraken:taxid|1235691 complete sequence, 68671 bp, circularly permuted, Cluster BN.fna",
#      Path("testing/Streptomyces_phage_Yara|kraken:taxid|1235691 complete sequence, 68671 bp, circularly permuted, Cluster BN.fna"))
#      ])
# def test_parseKrakenFile(parsed_phages, expected):
#     assert parser.parseKrakenFile(kraken_report_path) == expected

@pytest.mark.parametrize("phage_name, expected,", [
    ("Mycobacterium phage Paphu", True),
    ("mycobacterium phage paphu", True),
    ("Paphu", True),
    ("Kleenex", False)
    ])
def test_checkPhageExists(expected_phage_dict, phage_name, expected):
    assert parser.checkPhageExists(phage_name, expected_phage_dict) == expected


def test_saveTaxidToFile():
    pass


def test_saveNameToFile():
    pass
