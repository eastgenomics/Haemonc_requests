#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:42:37 2024

@author: arun
"""

import os
import re
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from plots_bcftools import list_samplenames


def give_tsv_df(sample_name, version, cwd, index_col):
    """
    Search for a VCF file matching the given sample name and version in the
    specified directory, and load it into a pandas DataFrame.

    Parameters
    ----------
    sample_name : str
        The name of the sample to search for in the tsv file.
    version : str
        The version identifier to match in the tsv file name.
    cwd : str
        The directory where the tsv files are located.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the contents of the VCF file, with
        predefined columns.

    Raises
    ------
    FileNotFoundError
        If no VCF file matching the given sample name and version is found in
        the specified directory.

    Examples
    --------
    >>> df = give_tsv_df('sample123', 'exon_stats', '/path/to/vcfs')

    """
    try:
        # Construct the regular expression pattern to match the VCF file
        pattern = f".*{sample_name}.+_{version}\\.tsv?$"
        # print(pattern)

        # Search for the first file in the directory matching the pattern.
        # If no file matches, return None instead of raising StopIteration.
        file = next(
            (filename for filename in os.listdir(cwd)
             if re.search(pattern, filename)), None
        )

        # Raise an exception if no matching file is found
        if file is None:
            raise FileNotFoundError(
                f"No file was found with {sample_name} and {version}"
            )

        # Load the VCF file into a pandas DataFrame
        dataframe = pd.read_csv(
            os.path.join(cwd, file), sep="\t", header=0, index_col=index_col
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return dataframe


def generate_genes_heatmap_df(directory_path):
    """
    Generate a heatmap using the gene_stats tsv files

    Parameters
    ----------
    directory_path : str
        Directory where the genes_stat.tsv files are stored

    Returns
    -------
    merged_df : pd.Dataframe
        Dataframe to use to produce heatmap.

    """

    # get the list of samples from directory
    list_of_samples = list_samplenames(directory_path)

    # Generate a df friendly heatmap
    first_df = give_tsv_df(list_of_samples[0], "gene_stats", directory_path, 0)
    first_col = first_df["250x"]
    merged_df = first_col.to_frame(name=list_of_samples[0])
    for sample in list_of_samples[1:]:
        print(sample)
        df = give_tsv_df(sample, "gene_stats", directory_path, 0)
        merged_df[sample] = df["250x"]

    return merged_df


def generate_exons_heatmap_df(directory_path):
    """
    Generate a heatmap using the gene_stats tsv files

    Parameters
    ----------
    directory_path : str
        Directory where the genes_stat.tsv files are stored

    Returns
    -------
    merged_df : pd.Dataframe
        Dataframe to use to produce heatmap.

    """

    # get the list of samples from directory
    list_of_samples = list_samplenames(directory_path)

    # Generate a df friendly heatmap
    first_df = give_tsv_df(list_of_samples[0], "exon_stats",
                           directory_path, [3, 4, 5])
    first_col = first_df["250x"]
    merged_df = first_col.to_frame(name=list_of_samples[0])
    for sample in list_of_samples[1:]:
        print(sample)
        df = give_tsv_df(sample, "exon_stats", directory_path, [3, 4, 5])
        merged_df[sample] = df["250x"]

    return merged_df


def plot_heatmap(dataframe, plot_width, plot_height, title, y_label):
    """
    Plot the heatmap given the dataframe and plot parameters

    Parameters
    ----------
    dataframe : pd.Dataframe
        dataframe from generate_genes_heatmap_df() or
        generate_exons_heatmapdf()
    plot_width : int
        Width of the heatmap plot
    plot_height : int
        Height of the heatmap plot
    title : str
        Heatmap title
    y_label : str
        Title on the Y axis of the heatmap plot

    Returns
    -------
    None.
        Prints out a plot

    """

    # Plot the heatmap
    plt.figure(figsize=(plot_width, plot_height))
    sns.heatmap(
        dataframe,
        annot=False,
        cmap="RdYlGn",
        fmt=".1f",
        linewidths=0.5,
        vmax=100,
        annot_kws={"size": 6},
        cbar_kws={"label": "Coverage Percentage"},
    )
    plt.title(title)
    plt.xlabel("Samples")
    plt.ylabel(y_label)
    plt.show()


def main():
    """Main function of the code"""
    # Specify the working directory
    working_dir = (
        "/home/arun/Codes/HaemOnc_panel_2024-Twist_capture"
        + "/Validation/plots_to_present/coverage_files"
    )

    # Heatmap header
    header = "Percentage coverage Heatmap for samples at 250x"

    # Generate the plots with all probes
    local_dir_tuple = (
        "/SP/old_probes",
        "/S1/old_probes",
        "/SP/new_probes",
        "/S1/new_probes",
    )

    for local_dir in local_dir_tuple:
        directory_path = working_dir + local_dir
        merged_df = generate_genes_heatmap_df(directory_path)
        plot_heatmap(merged_df, 12, 25, header, "Genes")

    # set a list of new genes added onto the new probe design
    genes_of_interest = [
        "ARAF",
        "ARID2",
        "ASXL2",
        "ATRX",
        "CD79B",
        "CSF1R",
        "CSNK1A1",
        "CTCF",
        "EIF6",
        "ERBB3",
        "FYN",
        "GNB1",
        "LUC7L2",
        "MAP3K1",
        "NTRK1",
        "PIGA",
        "PIK3CA",
        "PIK3CD",
        "PLCG1",
        "POT1",
        "PRPF8",
        "RPL11",
        "RPL5",
        "SF1",
        "SF3A1",
        "SMC1A",
        "SMC3",
        "SRCAP",
        "SUZ12",
        "U2AF2",
        "VAV1",
        "YLPM1",
        "ZBTB33",
        "ZNF318",
    ]

    local_dir_tuple = ("/SP/new_probes_with_new_genes",
                       "/S1/new_probes_with_new_genes")

    # Generate the SP and S1 plot with new probes and only new genes
    for local_dir in local_dir_tuple:
        directory_path = working_dir + local_dir
        merged_df = generate_genes_heatmap_df(directory_path)
        merged_df = merged_df.loc[genes_of_interest]
        plot_heatmap(merged_df, 10, 15, header, "Genes")

    # generate the SP and S1 plot with new probes, old genes
    # but with bigger flanking regions
    for local_dir in local_dir_tuple:
        directory_path = working_dir + local_dir
        merged_df = generate_genes_heatmap_df(directory_path)
        merged_df = merged_df.loc[~merged_df.index.isin(genes_of_interest)]
        plot_heatmap(merged_df, 12, 25, header, "Genes")

    # Generate the SP and S1 plots with the new probes,
    # but for NOTCH1 transcripts
    local_dir_tuple = ("/SP/new_probes", "/S1/new_probes")

    for local_dir in local_dir_tuple:
        directory_path = working_dir + local_dir
        merged_df = generate_exons_heatmap_df(directory_path)
        gene_of_interest = ["NOTCH1"]
        new_merged_df = merged_df[
            merged_df.index.get_level_values("gene").isin(gene_of_interest)
        ]
        plot_heatmap(new_merged_df, 12, 15, header, "Exons")


if __name__ == "__main__":
    main()
