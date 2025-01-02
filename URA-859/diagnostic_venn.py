#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 16:32:52 2024

@author: arun
"""

import os
import re
import subprocess
import pandas as pd
from plots_bcftools import venn_diagramify


def bcftools_isec(
    working_dir,
    sample_name_a,
    sample_name_b,
    version,
    folder_a,
    folder_b,
    output_folder,
):
    """
    Perform bcftools isec and return output vcfs as dataframes

    Parameters
    ----------
    working_dir : str
        Directory where the folder_a and folder_b is located
    sample_name_a : str
        Sample name for .vcf filename regex pattern recognition
    sample_name_b : str
        Sample name for the other .vcf filename regex pattern recognition
    version : str
        version of the vcf file for regex pattern recognition
    folder_a : str
        folder name where vcf_a files are located
    folder_b : str
        folder name where vcf_b files are located

    Raises
    ------
    FileNotFoundError
        DESCRIPTION.

    Returns
    -------
    vcfa_df_uniq : pd.Dataframe
        Variants unique to vcf_a
    vcfb_df_uniq : pd.Dataframe
        Variants unique to vcf_b
    vcfa_df_common : pd.Dataframe
        Variants from vcf_a shared with vcf_b
    vcfb_df_common : pd.Dataframe
        Variants from vcf_b shared with vcf_a

    """
    os.chdir(working_dir)

    # Construct the regular expression pattern to match the VCF file
    pattern_a = f".*{sample_name_a}.+_{version}\.vcf(\.gz)?$"
    pattern_b = f".*{sample_name_b}.+_{version}\.vcf(\.gz)?$"

    # Find the vcf files of interest
    file_a = next(
        (
            filename
            for filename in os.listdir(f"{working_dir}/{folder_a}")
            if re.search(pattern_a, filename)
        ),
        None,
    )

    file_b = next(
        (
            filename
            for filename in os.listdir(f"{working_dir}/{folder_b}")
            if re.search(pattern_b, filename)
        ),
        None,
    )

    if file_a is None or file_b is None:
        raise FileNotFoundError(
            f"No file was found with \
                                {sample_name_a} or {sample_name_b}"
        )

    # Set up relevant directories
    path_a = f"{folder_a}/{file_a}"
    path_b = f"{folder_b}/{file_b}"
    output_dir = os.path.join(output_folder)

    # Run bcftools isec
    isec_command = f"bcftools isec {path_a} {path_b} -p {output_dir}"
    print(isec_command)
    subprocess.run(isec_command, shell=True, check=True)

    # Define the columns for the VCF file DataFrame
    cols = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    ]

    # assign dataframes to return variables
    vcfa_df_uniq = pd.read_csv(
        os.path.join(output_dir, "0000.vcf"), sep="\t", comment="#", names=cols
    )
    vcfb_df_uniq = pd.read_csv(
        os.path.join(output_dir, "0001.vcf"), sep="\t", comment="#", names=cols
    )
    vcfa_df_common = pd.read_csv(
        os.path.join(output_dir, "0002.vcf"), sep="\t", comment="#", names=cols
    )
    vcfb_df_common = pd.read_csv(
        os.path.join(output_dir, "0003.vcf"), sep="\t", comment="#", names=cols
    )

    return vcfa_df_uniq, vcfb_df_uniq, vcfa_df_common, vcfb_df_common


def save_dataframe_as_tsv(df, file_path, index=False):
    """
    Saves a Pandas DataFrame as a TSV file.

    Parameters:
        df (pd.DataFrame): The DataFrame to save.
        file_path (str): The path where the TSV file will be saved.
        index (bool): Whether to include the index in the saved file.
    """
    try:
        df.to_csv(file_path, sep="\t", index=index)
        print(f"DataFrame saved as a TSV file at {file_path}")
    except Exception as e:
        print(f"An error occurred while saving the DataFrame: {e}")


def main():
    working_dir = (
        "/home/arun/Codes/HaemOnc_panel_2024-Twist_capture/"
        "Validation/plots_to_present/vcf_files"
    )
    version = "allgenesvep_filtered"

    samples_dict = [
        {
            "diagnostic": "22129Z0052",
            "validation": "22129Z0052A",
            "diagnostic_folder": "Diagnostic_runs/116495910-22129Z0052A",
            "validation_folder": "S1/new_probes_filtered",
        },
        {
            "diagnostic": "22129Z0052",
            "validation": "22129Z0052B",
            "diagnostic_folder": "Diagnostic_runs/116495910-22129Z0052B",
            "validation_folder": "S1/new_probes_filtered",
        },
        {
            "diagnostic": "23061Z0028",
            "validation": "23061Z0028",
            "diagnostic_folder": "Diagnostic_runs/121722703-23061Z0028",
            "validation_folder": "S1/new_probes_filtered",
        },
        {
            "diagnostic": "23273K0002",
            "validation": "23273K0002A",
            "diagnostic_folder": "Diagnostic_runs/125589267-23273K0002A",
            "validation_folder": "SP/new_probes_filtered",
        },
        {
            "diagnostic": "23273K0002",
            "validation": "23273K0002B",
            "diagnostic_folder": "Diagnostic_runs/125589267-23273K0002B",
            "validation_folder": "SP/new_probes_filtered",
        },
        {
            "diagnostic": "24032K0038",
            "validation": "24032K0038",
            "diagnostic_folder": "Diagnostic_runs/127898716-24032K0038",
            "validation_folder": "S1/new_probes_filtered",
        },
    ]

    for sample in samples_dict:

        output_folder = sample["diagnostic_folder"] + "/bcftools_test"

        vcfa_df_uniq, vcfb_df_uniq, vcfa_df_common, vcfb_df_common = bcftools_isec(
            working_dir,
            sample["diagnostic"],
            sample["validation"],
            version,
            sample["diagnostic_folder"],
            sample["validation_folder"],
            output_folder,
        )

        venn_diagramify(
            vcfa_df_uniq,
            vcfb_df_uniq,
            vcfa_df_common,
            "diagnostic run variants",
            "new probes variants",
            f"Sample {sample['validation']}",
        )

        output_dir = sample["diagnostic_folder"]
        output_filename = f"diagnostic_{sample['validation']}_uniq_var.tsv"
        tsv_file = f"{output_dir}/{output_filename}"
        save_dataframe_as_tsv(vcfa_df_uniq, tsv_file)


if __name__ == "__main__":
    main()
