from typing import List, TextIO
import pandas as pd
import os
from collections import defaultdict

'''

#### THIS WILL PROBABLY BE USEFUL IF THIS WAS NOT STANDALONE
import shutil
from apps.binning.unitem.unitem import UniteM
from apps.read_mapping.coverm.coverm import CoverM

PATH_PKGS = ['metabat', 'run_MaxBin.pl', "UniteM",
             "coverm"]
def check_dependencies(ref_pkgs: List[str]) -> None:
    """
    Check if all dependencies are on system path and are executable
    """
    acc = [shutil.which(pkg) for pkg in ref_pkgs if shutil.which(pkg) is not None]
    for missing_pkg in acc:
        print(f'{missing_pkg} is not installed')
    return None
'''

### FUNCTIONS TO IMPORT PF FORMAT FILES
def convert_pathway_row_list_format(pathway_row:str) -> List[str]:
    """
    Convert a row in a pathway entry to a list of two strings: one of the type of information, and the information itself
    """
    end_idx = 1
    while not pathway_row[end_idx] == '\t':
        end_idx += 1
    return [pathway_row[:end_idx], pathway_row]


def convert_pl_input_to_pathways_df(pl: str) -> pd.DataFrame:
    """
    Use the name from pl to open a pathway tools file of the same name
    Read the specific file format of pathway tools output from a Metapathways output.
    """

    # Convert pl into a list with each row as an item
    with open(pl) as f:
        pl_lines = [line.rstrip() for line in f]

    # Used a two pointer style algo to determine where each pathway starts and ends,
    # then input into a dict accumulator
    start_idx = 0  # Starting index of rows in a single pathway
    pathway_acc = [] # Accumulator that we input each single pathway into
    test_acc = defaultdict(lambda:0)
    while start_idx < len(pl_lines):
        end_idx = start_idx  # Ending index of rows in a single pathway
        temp_pathway = {'ec': [], 'metacyc': []}  # Default values for ec and ko, # since these do not appear in every pathway

        while pl_lines[end_idx] != "//":  # Check every line to see if it's an end of a pathway entry
            entry_list_format = convert_pathway_row_list_format(pl_lines[end_idx])
            test_acc[entry_list_format[0].lower()] += 1

            # Special case for rows signifying ec or metacyc accessions,
            # since there can be [0,inf] of these in a pathway
            if entry_list_format[0] == "EC":
                temp_pathway['ec'].append(entry_list_format[1])
            elif entry_list_format[0] == "METACYC":
                temp_pathway['metacyc'].append(entry_list_format[1])
            else:
                temp_pathway[entry_list_format[0]] = entry_list_format[1]
            end_idx += 1
        pathway_acc.append(temp_pathway)
        start_idx += (end_idx+1-start_idx)
    return pd.DataFrame(pathway_acc)


### FUNCTIONS TO IMPORT ORF CONTIG ANNOTATION MAP
def convert_orf_contig_map_to_df(map:str) -> pd.DataFrame:
    """
    Read the ORF annotation table and change the name of ORF ID to "ID" for easier joins
    :param map:
    :return:
    """

    df = pd.read_csv(map, sep = '\t')
    df = df.rename(columns={'# ORF_ID':'ID'})
    df = df.assign(ID='ID\tO'+df['ID'])
    return df
    # Prepend required syntax to join to pathway input
###



def main():
    pass



if __name__ == '__main__':
    pass"""