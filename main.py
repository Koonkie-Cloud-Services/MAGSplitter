from input_to_df import *

'''
#####
Inputs

0.pf: pathway tools initial input from metapathways outputs
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

# Since we now have the names of deleted ORFs in a clean format, we can add
# deleted ORF reactions back into the rxn dataframe
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
                                (df_rxn.loc[df_rxn['ORF_ID'] == intact_orf].
                                 assign(ORF_ID=rxn[deleted_orf],
                                        NAME=rxn[deleted_orf]))], ignore_index=True)
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
    :param df_rxn: rxn dataframe
    :param df_mag_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """
    return df_rxn.merge(df_mag_map, on='CONTIG_ID', how='left')


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

def main():
    # Import MP/WGS pipeline outputs
    df_pf = convert_pl_input_to_rxn_df('example/0.pf')
    df_duplicate_orf_map = convert_duplicate_orf_map_to_list('example/orf_map.txt')
    orf_contig_map = convert_orf_contig_map_to_df("example/GAPP-5498e568-6918-4000-b27e-dbeff35eeee7.ORF_annotation_table.txt")
    contig_mag_map = convert_contig_mag_map_to_df("example/""contig_info.tsv")
    df_pf = undo_orf_removal(df_pf, df_duplicate_orf_map)
    df_pf_with_contig = combine_rxn_contig_map(df_pf, orf_contig_map)
    df_pf_with_mag = combine_rxn_mag_map(df_pf_with_contig, contig_mag_map)

    # Parse arguments
    """
    parser = argparse.ArgumentParser(description='Convert metapathways output to ePGDB readable format')
    parser.add_argument('-i', '--input', help='metapathways output name (normally 0.pf)', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-m', '--metacyc', help='metacyc accession map', required=True)
    parser.add_argument('-c', '--contig_mag_map', help='contig mag map', required=True)
    parser.add_argument('-e', '--orf_map', help='orf map', required=True)
    args = parser.parse_args()

if __name__ == '__main__':
    main()


'''
### Testing arguments



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