import pandas as pd
import dask.dataframe as dd
from typing import List
import time
from multiprocessing import Pool
import itertools

# These functions are to combine each of the mapping and ePGDB dataframes into a list of dataframes split by MAGs.


# Since we now have the names of deleted ORFs in a clean format, we can add
# deleted ORF reactions back into the rxn dataframe
def undo_orf_removal_test(df_rxn: pd.DataFrame, duplicate_orf_map: List[List[str]]) -> pd.DataFrame:
    """
    Add in extra rows that were deleted from metapathways output, done previously to optimize usage of pathway tools
    NOT RUNTIME OPTIMIZED
    :param df_rxn: rxn dataframe
    :param duplicate_orf_map: duplicate orf list of lists
    :return: rxn dataframe with extra rows
    """
    # Open the orf map file, and only choose the maps that have at least two ORFs in each row
    # The ones with only one ORF indicate that there are no duplicates of that rxn
    start = time.perf_counter()
    # We use rxn_acc to accumulate all the new rows that we are going to add to the df_rxn dataframe
    rxn_acc = []
    for rxn in duplicate_orf_map:
        intact_orf = rxn[0]
        found_rxn_row = None  # This will change when the row in df_rxn is found for first time in set of matching rxns
        for deleted_orf in range(1, len(rxn)):
            # Create a new dataframe with only the intact ORF, change the ID to a deleted ORF, and accumulate them
            if found_rxn_row is None:
                found_rxn_row = df_rxn[df_rxn['ORF_ID'] == intact_orf]
                # Case where no matching ORF is found
                if found_rxn_row.empty:
                    raise Exception("No matching ORF found for ORF_ID: " + intact_orf)
                    break
            rxn_acc.append(
                found_rxn_row.
                    assign(ORF_ID=rxn[deleted_orf],
                           NAME="NAME" + rxn[deleted_orf][2:]))
    # Combine rxn dataframe with rxn accumulator dataframe
    rxn_acc = pd.concat(rxn_acc, ignore_index=True)
    df_rxn = pd.concat([df_rxn, rxn_acc], ignore_index=True)
    end = time.perf_counter()
    print("Time to undo orf removal: " + str(end - start))
    return df_rxn

# Combine the rxn dataframe with the orf:contig and contig:mag dataframes
def combine_rxn_contig_map(df_rxn: pd.DataFrame, df_contig_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig map dataframe

    :param df_rxn: rxn dataframe
    :param df_contig_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """
    return df_rxn.merge(df_contig_map, on='ORF_ID', how='left')


def combine_rxn_mag_map(df_rxn: pd.DataFrame, df_mag_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig MAG map dataframe
    Replace non-binned ORFs with BIN ID of 'non_binned_metagenome'
    :param df_rxn: rxn dataframe
    :param df_mag_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """

    df = df_rxn.merge(df_mag_map, on='CONTIG_ID', how='left')
    df['BIN_ID'] = df["BIN_ID"].fillna(value='non_binned_metagenome')
    return df


def split_full_rxn_df_by_mag(df_rxn: pd.DataFrame) -> dict:
    """
    Split the full rxn dataframe into a dict of rxn dataframes, one for each MAG {BIN_ID: df}
    For all rxns that are not in a MAG, put them in a separate dataframe (metagenome not in MAG)
    :param df_rxn: full rxn dataframe
    :return: list of rxn dataframes, one for each MAG
    """
    mags = df_rxn['BIN_ID'].unique()
    df_rxn_list = {}
    for mag in mags:
        df_rxn_list[mag] = df_rxn[df_rxn['BIN_ID'] == mag]

    return df_rxn_list


# Misc functions
def sample_name_grabber(rxn_df: pd.DataFrame) -> pd.DataFrame:
    """
    Grab the sample name from the contig dataframe.  This is under the assumption
    that the sample name is the same for all contigs in a MAG.
    :param rxn_df: 
    :return: 
    """
    # end index is one character after the last character in the string
    full_str = rxn_df['CONTIG_ID'].iloc[0]
    end_idx = len(full_str)
    # If we only have the sample name without contig ID, turn <is_only_sample> to True
    is_only_sample = False
    while end_idx > 0 and is_only_sample is False:
        if full_str[end_idx - 1] == '_':
            is_only_sample = True
            end_idx -= 1
        else:
            end_idx -= 1
    return full_str[:end_idx]


if __name__ == '__main__':
    pass
