from typing import List
import pandas as pd
import os
from collections import defaultdict

'''

#### DEPENDENCY CHECK WILL PROBABLY BE USEFUL IF THIS WAS NOT STANDALONE
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


# Within the metapathways output, duplicated pathways were removed as frequency does not contribute to ePGDB readings,
# and only slows down the runtime of pathway tools.
# Therefore, we are doing this to add back the deleted pathways using the map of deletions
def undo_orf_removal(df_pathway: pd.DataFrame, orf_map_name:str) -> pd.DataFrame:
    """
    Add in extra rows that were deleted from metapathways output, done previously to optimize usage of pathway tools
    NOT RUNTIME OPTIMIZED

    :param orf_map_name:
    :param df_pathway:
    :param orf_map:
    :return:
    """

    # Open the orf map file, and only choose the maps that have at least two ORFs in each row
    # The ones with only one ORF indicate that there are no duplicates of that pathway
    # Creates a list lists, with each item after the first in a nested list showing deleted duplicates of the first item
    with open(orf_map_name) as f:
        orf_map = [(line.rstrip().split("\t"))
                   for line in f
                   if ('\t' in line)]

    #Prepend each orf with the required pathway tools syntax
    for pathway_idx in range(len(orf_map)):
        for orf_idx in range(len(orf_map[pathway_idx])):
            orf_map[pathway_idx][orf_idx] = 'ID\t' + orf_map[pathway_idx][orf_idx]

    # Since we now have the names of deleted ORFs in a clean format, we can add back into the pathways dataframe
    for pathway in orf_map:
        intact_orf = pathway[0]
        for deleted_orf in range(1, len(pathway)):
            # Create a new dataframe with only the intact ORF, change the ID to a deleted ORF, then add it back into
            # our big dataframe
            df_pathway = pd.concat([df_pathway,
                                    (df_pathway.loc[df_pathway['ID'] == intact_orf].
                                     assign(ID=pathway[deleted_orf],
                                            NAME = pathway[deleted_orf]))])

    return df_pathway



### FUNCTION TO IMPORT ORF CONTIG ANNOTATION MAP
def convert_orf_contig_map_to_df(map:str) -> pd.DataFrame:
    """
    Read the ORF annotation table and change the name of ORF ID to "ID" for easier joins
    :param map:
    :return:
    """
    df = pd.read_csv(map, sep = '\t')
    df = df.rename(columns={'# ORF_ID':'ID'})
    df = df.assign(ID='ID\tO_'+df['ID'])
    return df

def combine_pathways_contig_map(df_pathways: pd.DataFrame, df_contig_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the pathways dataframe with the contig map dataframe
    :param df_pathways:
    :param df_contig_map:
    :return:
    """
    return df_pathways.merge(df_contig_map, on='ID', how='left')

def main():
    import cProfile
    import pstats
    with cProfile.Profile() as pr:
        test1 = convert_pl_input_to_pathways_df('0.pf')
        undo_orf_removal(test1, 'orf_map.txt')
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.print_stats()



if __name__ == '__main__':
    main()