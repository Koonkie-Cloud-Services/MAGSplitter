import pandas as pd
import os
from typing import List
from shutil import rmtree


# File creator functions are to create all the pf files required for pathway tools in a pf folder
def pf_file_creator(pf: pd.DataFrame, file_path: str) -> None:
    """
    Create a pf file from a ORF rxn dataframe.
    Each rxn has the following information written sequentially:
    1. ORF_ID
    2. NAME
    3. STARTBASE
    4. ENDBASE
    5. FUNCTION
    6. METACYC: optional with len [0,inf]
    6. EC: optional with len [0, inf]
    6. PRODUCT TYPE
    7. "//"
    :param file_path: String representing where the file will be written
    :param pf: single rxn dataframe
    :return: None
    """
    f = open(f"{file_path}/0.pf", "a")
    # Iterate through every ORF in the rxn dataframe
    for row in pf.itertuples():
        f.write(row.ORF_ID + "\n")  # Write the ORF ID to the pf file
        f.write(row.NAME + "\n")  # Write the ORF name to the pf file
        f.write(row.STARTBASE + "\n")  # Write the ORF start base to the pf file
        f.write(row.ENDBASE + "\n")  # Write the ORF end base to the pf file
        f.write(row.FUNCTION + "\n")  # Write the ORF function to the pf file
        # If the ORF has a metacyc accession, write it to the pf file
        for metacyc in row.metacyc:
            f.write(metacyc + "\n")
        # If the ORF has an EC number, write it to the pf file
        for ec in row.ec:
            f.write(ec + "\n")
        f.write(row.PRODUCT_TYPE + "\n")  # Write the ORF product type to the pf file
        f.write("//\n")  # Signify the end of the ORF to the pf file
    f.close()
    return None


def genetic_elements_creator(file_path: str) -> None:
    """
    Create a genetic-elements.dat file for pathway tools.
    :param file_path: String representing where the file will be written
    :return: None
    """
    f = open(f"{file_path}/genetic-elements.dat", "a")
    f.write("ID\t0\n")
    f.write("NAME\t0\n")
    f.write("TYPE\t:CONTIG\n")
    f.write("ANNOT-FILE\t0.pf\n")
    f.write("//\n")
    f.close()
    return None


def organism_params_creator(file_path: str, sample_name: str) -> None:
    """
    Create an organism-params.dat file for pathway tools.
    Just copied and pasted from metapathways output, not sure if this will be right. Ask Ryan later
    :param file_path: String representing where the file will be written
    :param sample_name: String representing the sample name
    :return: None
    """
    f = open(f"{file_path}/organism-params.dat", "a")
    f.write(f"ID\t{sample_name}\n")
    f.write("STORAGE FILE\n")
    f.write(f"ABBREV-NAME\t{sample_name}\n")
    f.write("STRAIN\t1\n")
    f.write("RANK\t|species|\n")
    f.write("NCBI-TAXON-ID\t12908\n")
    f.close()
    return None


def dummy_file_creator(file_path: str, sample_name: str) -> None:
    """
    Write a dummy file in the metapathways output folder for pathway tools
    :param sample_name: name of sample inputted into metapathways
    :param file_path: String representing where the file will be written
    :return: None
    """
    f = open(f"{file_path}/{sample_name}.dummy.txt", "a")
    f.close()
    return None


def ptools_folder_creator(target_folder, sample_name:str, rxn_processed: dict[pd.DataFrame]) -> None:
    """
    Main program to create all folders of pf files, to be inputted into pathway tools,
    from a list of rxn dataframes.
    The pf folder includes "0.pf" (the main reaction ORF file), written with pf_file_creator,
    "genetic-elements.dat" (pre-written file, which I should ask Ryan about), written with genetic_elements_creator
    "organism-params.dat" (which I will also ask Ryan about), written with organism_params_creator
    "pathologic.log" (pre-written file, which I should ask Ryan about), written with pathologic_log_creator
    "<sample_name>.dummy.txt" (dummy file),
    :param target_folder: Program main folder (will be user defined later)
    :param rxn_processed: dict of mag rxn dataframes
    :return: None
    """
    # Create the results folder
    results_path = os.path.join(target_folder, "results")
    if os.path.exists(results_path) and os.path.isdir(results_path):
        rmtree(results_path)
    os.mkdir(results_path)

    # Iterate through every single mag df in the dict of mag rxn dataframes
    for mag_df in rxn_processed:
        # Create the folder for the current mag df
        mag_path = os.path.join(results_path, mag_df)
        os.mkdir(mag_path)
        # Create the pf file for the current mag df
        pf_file_creator(rxn_processed[mag_df], mag_path)
        # Create the genetic elements file
        genetic_elements_creator(mag_path)
        # Create the organism params file
        organism_params_creator(mag_path, sample_name)
        # Create the dummy file
        dummy_file_creator(mag_path, sample_name)
    return None


if __name__ == "__main__":
    pass
