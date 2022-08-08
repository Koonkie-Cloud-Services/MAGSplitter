import argparse
from collections import defaultdict
from typing import List
import pandas as pd

'''
#####
Inputs

0test.pf: pathway tools initial input from metapathways outputs
(found in /ptools/)

<sample_name>.ORF_annotation_tabel.txt: map of ORFS and contigs for an environment metagenome, outputted by metapathways 
(found in /results/annotation_table/)

orf_map.txt: map of duplicated ORFS within a metagenome, with the first ORF being the intact ORF, 
and the rest being the duplicates
(found in /results/annotation_table/)

config_info.tsv: map of contigs to MAGs created through the WGS pipeline binning process.  
(found in /binning/results/greedy)
#####
'''


# FUNCTIONS TO IMPORT PF FORMAT FILES

def convert_rxn_row_list_format(rxn_row: str) -> List[str]:
    """
    Convert a row in a reaction database to a list of two strings: the type of information, and the information itself
    :param rxn_row: reaction entry from rxn tools output
    :return: list of two strings: one of the type of information, and the information itself
    """
    end_idx = 1
    while not rxn_row[end_idx] == '\t':
        end_idx += 1
    return [rxn_row[:end_idx], rxn_row]


def convert_pl_input_to_rxn_df(pl: str) -> pd.DataFrame:
    """
    Use the name from pl to open a pathway tools file of the same name
    Read the specific file format of pathway tools output from a Metapathways output.
    :param pl: pathway tools input file name
    :return: rxn dataframe
    """

    # Convert pl into a list with each row as an item
    with open(pl) as f:
        pl_lines = [line.rstrip() for line in f]

    # Used a two pointer style algo to determine where each rxn starts and ends,
    # then input into a dict accumulator
    start_idx = 0  # Starting index of rows in a single rxn
    rxn_acc = []  # Accumulator that we input each single rxn into
    test_acc = defaultdict(lambda: 0)
    while start_idx < len(pl_lines):
        end_idx = start_idx  # Ending index of rows in a single rxn
        temp_rxn = {'ec': [], 'metacyc': []}  # Default values for ec and ko, # since these do not appear in every rxn
        while pl_lines[end_idx] != "//":  # Check every line to see if it's an end of a rxn entry
            entry_list_format = convert_rxn_row_list_format(pl_lines[end_idx])
            test_acc[entry_list_format[0].lower()] += 1
            # Special case for rows signifying ec or metacyc accessions,
            # since there can be [0,inf] of these in a rxn
            if entry_list_format[0] == "EC":
                temp_rxn['ec'].append(entry_list_format[1])
            elif entry_list_format[0] == "METACYC":
                temp_rxn['metacyc'].append(entry_list_format[1])
            else:
                temp_rxn[entry_list_format[0]] = entry_list_format[1]
            end_idx += 1
        rxn_acc.append(temp_rxn)
        start_idx += (end_idx + 1 - start_idx)
    return pd.DataFrame(rxn_acc)


# FUNCTION TO IMPORT ORF CONTIG ANNOTATION MAP
def convert_orf_contig_map_to_df(orf_map: str) -> pd.DataFrame:
    """
    Read the ORF annotation table and change the name of ORF ID to "ID" for easier joins
    :param orf_map: ORF annotation table
    :return: ORF annotation table as a dataframe
    """

    df = pd.read_csv(orf_map, sep='\t')
    df = df.rename(columns={'# ORF_ID': 'ID'})
    df = df.assign(ID='ID\tO_' + df['ID'])
    return df[["ID", "CONTIG_ID"]]  # Return only the two columns we need


def convert_contig_mag_map_to_df(contig_mag_map: str) -> pd.DataFrame:
    """
    Read the contig mag map, and extract the UniteM Bin ID and Contig ID
    :param contig_mag_map: contig mag map
    :return: contig mag map as a dataframe
    """

    df = pd.read_csv(contig_mag_map, sep='\t')
    df = df.rename(columns={'UniteM Bin ID': 'Bin_ID', 'Contig ID': 'Contig_ID'})
    return df[['Bin_ID', 'Contig_ID']]


# Within the metapathways output, duplicated rxns were removed as frequency does not contribute to ePGDB readings,
# and only slows down the runtime of pathway tools.
# Therefore, we are doing this to add back the deleted pathways using the map of deletions
def convert_duplicate_orf_map_to_list(duplicate_orf_map: str) -> List[List[str]]:
    """
    Read the duplicate ORF map, and create a list of lists,
    with each item after the first in a nested list showing deleted duplicates of first item
    :param duplicate_orf_map: duplicate orf map
    :return: duplicate orf map as a dataframe
    """
    with open(duplicate_orf_map) as f:
        orf_map = [(line.rstrip().split("\t"))
                   for line in f
                   if ('\t' in line)]
    # Prepend each orf with the required rxn tools syntax
    for rxn_idx in range(len(orf_map)):
        for orf_idx in range(len(orf_map[rxn_idx])):
            orf_map[rxn_idx][orf_idx] = 'ID\t' + orf_map[rxn_idx][orf_idx]
    return orf_map


# Since we now have the names of deleted ORFs in a clean format, we can add back into the rxn dataframe
def undo_orf_removal(df_rxn: pd.DataFrame, duplicate_orf_map: List[List[str]]) -> pd.DataFrame:
    """
    Add in extra rows that were deleted from metapathways output, done previously to optimize usage of pathway tools
    NOT RUNTIME OPTIMIZED
    :param df_rxn: rxn dataframe
    :param duplicate_orf_map: duplicate orf list of lists
    :return: rxn dataframe with extra rows
    """
    # Open the orf map file, and only choose the maps that have at least two ORFs in each row
    # The ones with only one ORF indicate that there are no duplicates of that rxn
    for rxn in duplicate_orf_map:
        intact_orf = rxn[0]
        for deleted_orf in range(1, len(rxn)):
            # Create a new dataframe with only the intact ORF, change the ID to a deleted ORF, then add it back into
            # our big dataframe
            df_rxn = pd.concat([df_rxn,
                                (df_rxn.loc[df_rxn['ID'] == intact_orf].
                                 assign(ID=rxn[deleted_orf],
                                        NAME=rxn[deleted_orf]))], ignore_index=True)
    return df_rxn


def combine_rxn_contig_map(df_rxn: pd.DataFrame, df_contig_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig map dataframe
    :param df_rxn: rxn dataframe
    :param df_contig_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """
    return df_rxn.merge(df_contig_map, on='ID', how='left')


def combine_rxn_mag_map(df_rxn: pd.DataFrame, df_mag_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig MAG map dataframe
    :param df_rxn: rxn dataframe
    :param df_mag_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """
    return df_rxn.merge(df_mag_map, on='Bin_ID', how='left')


def main():
    df_pf = convert_pl_input_to_rxn_df('example/0test.pf')
    df_duplicate_orf_map = convert_duplicate_orf_map_to_list('example/orf_map.txt')
    orf_contig_map = convert_orf_contig_map_to_df(
        "example/GAPP-5498e568-6918-4000-b27e-dbeff35eeee7.ORF_annotation_table.txt")
    contig_mag_map = convert_contig_mag_map_to_df("example/""contig_info.tsv")
    df_pf = undo_orf_removal(df_pf, df_duplicate_orf_map)
    df_pf = combine_rxn_contig_map(df_pf, contig_mag_map)


    # Parse arguments
    """
    parser = argparse.ArgumentParser(description='Convert metapathways output to ePGDB readable format')
    parser.add_argument('-i', '--input', help='metapathways output name (normally 0.pf)', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-m', '--metacyc', help='metacyc accession map', required=True)
    parser.add_argument('-c', '--contig_mag_map', help='contig mag map', required=True)
    parser.add_argument('-e', '--orf_map', help='orf map', required=True)
    args = parser.parse_args()

    import cProfile
    import pstats
    with cProfile.Profile() as pr:
        test1 = convert_pl_input_to_rxn_df('0test.pf')
        undo_orf_removal(test1, 'orf_map.txt')
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.print_stats()"""


if __name__ == '__main__':
    main()


'''
### Testing arguments

def test_convert_pl_input_to_rxn_df():

# DEPENDENCY CHECK WILL PROBABLY BE USEFUL IF THIS WAS NOT STANDALONE
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
    return None'''


"""
TESTCASE

)"""