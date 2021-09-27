import pytest
import sys
import os
import parseKraken as parser
from typing import Dict
from pathlib import Path

''' Fixture setup '''
@pytest.fixture
def expected_phage_dict():
    return dict({ "2502430": "Mycobacterium phage Fushigi", 
            "2599873": "Mycobacterium phage Paphu", 
            "1913110": "Mycobacterium phage StarStuff"})

@pytest.fixture
def kraken_report_path():
    return Path(os.getcwd()) / Path("kraken_testing") 


@pytest.mark.parametrize("kraken_file, expected", [
    ("kraken_SAMPLE.report", { "2502430": "Mycobacterium phage Fushigi", "2599873": "Mycobacterium phage Paphu", "1913110": "Mycobacterium phage StarStuff"}),
    ("kraken_NO_SPECIES.report", {}),
    ("kraken_SHORT_NAMES.report", { "2502430": "Fushigi", "2599873": "Paphu", "1913110": "StarStuff"})
     ])
def test_parseKrakenFile(kraken_report_path, kraken_file, expected):
    assert parser.parseKrakenFile(kraken_report_path / Path(kraken_file)) == expected


def test_parseKrakenFile_nonexistent(kraken_report_path):
    ''' FileNotFound test '''
    with pytest.raises(FileNotFoundError):
        parser.parseKrakenFile(kraken_report_path/Path("kraken_FAIL.report"))


@pytest.mark.parametrize("phage_name, expected", [
    ("Mycobacterium phage Paphu", True),
    ("mycobacterium phage paphu", True),
    ("Paphu", True),
    ("Kleenex", False) 
    ])
def test_checkPhageExists(expected_phage_dict, phage_name, expected):
    assert parser.checkPhageExists(phage_name, expected_phage_dict) == expected


@pytest.mark.parametrize("phages_dict, expected", [
    ({"2502430": "Mycobacterium phage Fushigi", "2599873": "Mycobacterium phage Paphu", "1913110": "Mycobacterium phage StarStuff"}, Path("test_out.txt")),
    ({}, None)
])
def test_saveTaxidToFile(phages_dict, expected):
    actual_outfile = parser.saveTaxidToFile("test_out.txt", phages_dict)
    # first test: return value
    assert actual_outfile == expected 
    # second test: file contents
    if actual_outfile != None:
        contents = open(actual_outfile).readlines()
        for file_phage, dict_phage in zip(contents, phages_dict.keys()):
            assert file_phage.strip() == dict_phage.strip()


@pytest.mark.parametrize("phages_dict, expected", [
    ({"2502430": "Mycobacterium phage Fushigi", "2599873": "Mycobacterium phage Paphu", "1913110": "Mycobacterium phage StarStuff"}, Path("test_out.txt")),
    ({}, None)
])
def test_saveNameToFile(phages_dict, expected):
    actual_outfile = parser.saveNameToFile("test_out.txt", phages_dict)
    # first test: return value
    assert actual_outfile == expected 
    # second test: file contents
    if actual_outfile != None:
        contents = open(actual_outfile).readlines()
        for file_phage, dict_phage in zip(contents, phages_dict.values()):
            assert file_phage.strip() == dict_phage.strip()


@pytest.mark.skip
def test_main():
    assert parser.main()