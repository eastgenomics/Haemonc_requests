#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 12:09:14 2024

@author: arun
"""
import os
import pandas as pd
import re
import seaborn as sns
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def give_vcf_df(sample_name, version, cwd):
    """
    Search for a VCF file matching the given sample name and version in the specified directory,
    and load it into a pandas DataFrame.

    Parameters
    ----------
    sample_name : str
        The name of the sample to search for in the VCF file.
    version : str
        The version identifier to match in the VCF file name.
    cwd : str
        The directory where the VCF files are located.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the contents of the VCF file, with predefined columns.

    Examples
    --------
    >>> df = give_vcf_df('sample123', 'v1.2.3', '/path/to/vcfs')

    """
    filename = str(cwd) + str(f'/sentieon-tnbam_{version}/') + str(f'{sample_name}_markdup_recalibrated_tnhaplotyper2_normalised.vcf.gz')

    # Define the columns for the VCF file DataFrame
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            "SAMPLE"]
    # Load the VCF file into a pandas DataFrame
    df = pd.read_csv(filename, sep= "\t", comment='#',
                            names=cols)

    return df


def transformed_vcf(vcf_dataframe):
    """
   Transforms a VCF DataFrame by adding new columns for site, allele frequency (AF), and depth (DP).

    This function performs the following transformations:
    1. Adds a "site" column combining "CHROM", "POS", "REF", and "ALT" columns in the format: 'CHROM:POS_REF/ALT'.
    2. Extracts the allele frequency (AF) from the "SAMPLE" column and adds it as a new "AF" column.
    3. Extracts the depth (DP) from the "INFO" column and adds it as a new "DP" column.

    Parameters
    ----------
    vcf_dataframe : pd.DataFrame
        The input VCF DataFrame containing columns: "CHROM", "POS", "REF", "ALT", "SAMPLE", and "INFO".

    Returns
    -------
    pd.DataFrame
        A transformed VCF DataFrame with added "site", "AF", and "DP" columns.
    """

    # Create an identical dataframe
    vcf_df = vcf_dataframe

    # Create the Site column
    vcf_df['site'] = vcf_df['CHROM'].astype(str) + ':' \
                    + vcf_df['POS'].astype(str)  + '_' \
                    + vcf_df['REF'].astype(str) + '/' \
                    + vcf_df['ALT'].astype(str)

    # Create the AF
    vcf_df['AF'] = vcf_df['SAMPLE'].str.split(':', expand=True)[2]

    # Create the DP
    vcf_df['DP'] = vcf_df['INFO'].str.split(';').apply(
                   lambda x: ''.join((y for y in x if y.startswith('DP='))))
    vcf_df['DP'] = vcf_df['DP'].str.replace('DP=', '')

    return vcf_df


def merge_vcf_common_sites(vcfa_df, vcfb_df, version_a, version_b, sample):
    """
    Merge two vcf files based on their common sites

    Parameters
    ----------
    vcfa_df : pd.DataFrame
        vcf dataframe which is output of transformed_vcf(vcf_dataframe)
    vcfb_df : pd.Dataframe
        second vcf dataframe which is output of transformed_vcf(vcf_dataframe)
    version_a : str
        Version name of the vcfa_df file (eg. 'v3.2.0')
    version_b : str
        Version name of the vcfb_df file (e.g. 'v5.0.1')

    Returns
    -------
    pd.DataFrame
        Merged dataframe with the given 2 vcf files

    """

    # Define function where it returns list of common elements of two vcfs in
    # the site column
    def common_elements(list1, list2):
        return list(set(list1) & set(list2))

    common_sites = common_elements(vcfa_df.site,vcfb_df.site)
    print("Number of Chrom:Pos sites shared between two vcfs: "
          + str(len(common_sites)))

    # create dataframes containing common site rows
    vcfa_df_common = vcfa_df.loc[vcfa_df['site'].isin(common_sites)]
    vcfb_df_common = vcfb_df.loc[vcfb_df['site'].isin(common_sites)]

    # Rename the DP column to include their corresponding version
    vcfa_df_common = vcfa_df_common.rename({'DP': f'DP_{version_a}'}, axis=1)
    vcfb_df_common = vcfb_df_common.rename({'DP': f'DP_{version_b}'}, axis=1)

    vcfa_df_common = vcfa_df_common.rename({'AF': f'AF_{version_a}'}, axis=1)
    vcfb_df_common = vcfb_df_common.rename({'AF': f'AF_{version_b}'}, axis=1)

    # Merge both vcf files
    vcf_merged = pd.merge(vcfa_df_common, vcfb_df_common, on='site')

    # Store DP values as integers
    vcf_merged[f'DP_{version_a}'] = vcf_merged[f'DP_{version_a}'].astype(int)
    vcf_merged[f'DP_{version_b}'] = vcf_merged[f'DP_{version_b}'].astype(int)

    # rename the columns more clearly
    vcf_merged[f'AF_{version_a}'] = vcf_merged[f'AF_{version_a}'].astype(float)
    vcf_merged[f'AF_{version_b}'] = vcf_merged[f'AF_{version_b}'].astype(float)


    vcf_merged['AF_diff']  = abs(vcf_merged[f'AF_{version_a}'] - vcf_merged[f'AF_{version_b}'])
    vcf_merged2 = vcf_merged[(vcf_merged[['AF_diff']]>0.001).all(axis=1)]
    vcf_merged2.to_csv(f'plots/{sample}_v3_v5_AF_diff.tsv', index=False, sep = "\t")
    return vcf_merged


def make_correlation_plot_from_merged_vcf(vcf_merged, x_col, y_col, title='',
                                          xlabel=None, ylabel=None):
    """
    Generate a a correlation plot given the dataframe and 2 numerical columns

    Parameters
    ----------
    vcf_merged : pd.Dataframe
        Technicall any dataframe where you want to generate the plot from
    x_col : str
        The column name which will be the x of the plot
    y_col : str
        The column name which will be the y of the plot
    title : str, optional
        Title to include on the plot, sample should be specified here.
        The default is ''.
    xlabel : str, optional
        The x heading of the plot. The default is equal to x_col
    ylabel : str, optional
        The y heading of the plot. The default is equal to y_col

    Returns
    -------
    plot : matplotlib.axes._axes.Axes
        returns a callable pot

    """

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col

    # Calculate the R and p values
    r_value, p_value = pearsonr(vcf_merged[x_col], vcf_merged[y_col])
    print(r_value, p_value)
    n_value = len(vcf_merged.index)
    # Plot the correlation plot
    plot = sns.regplot(x=x_col, y=y_col, data=vcf_merged)
    plot.set_title(title)
    plot.set_xlabel(xlabel)
    plot.set_ylabel(ylabel)
    plot.text(x=min(vcf_merged[x_col]),
              y=max(vcf_merged[y_col]),
              s=f'R = {r_value:.2f}, p = {p_value:.2e}, n = {n_value}',
              ha='left', va='top')
    title_file = "plots/" + title + "_" + str(x_col) + "_" + str(y_col) + ".jpg"
    plt.savefig(title_file, dpi=300)
    plt.close()


samples = ["129325254-24102K0083-24NGSHO17-8128-F-96527893_S11_L001",
"129404211-24107K0073-24NGSHO18-8128-M-96527893_S46_L001",
"129450697-24109K0025-24NGSHO18-8128-M-96527893_S35_L001",
"129459591-24109K0055-24NGSHO18-8128-F-96527893_S44_L001",
"TMv2OCT20-HD734-1st_S35_L001",
"TMv2OCT20-HD734-2nd_S36_L001",
"TMv2OCT20-HD829-1st_S33_L001",
"TMv2OCT20-HD829-2nd_S34_L001",
"TMv2OCT20-oncospan-cell-line-1st_S39_L001",
"TMv2OCT20-oncospan-cell-line-2nd_S40_L001"]

for sample in samples:
    cwd = os.getcwd()
    vcf3_df = give_vcf_df(sample, "v3.2.0", cwd)
    vcf5_df = give_vcf_df(sample, "v5.0.1", cwd)

    vcf3_df = transformed_vcf(vcf3_df)
    vcf5_df = transformed_vcf(vcf5_df)

    vcf_merged = merge_vcf_common_sites(vcf3_df, vcf5_df, 'v3', 'v5', sample)
    vcf_merged.to_csv(f'plots/{sample}_v3_v5.tsv', index=False, sep = "\t")

    make_correlation_plot_from_merged_vcf(vcf_merged, "DP_v3", "DP_v5", sample)
    make_correlation_plot_from_merged_vcf(vcf_merged, "AF_v3", "AF_v5", sample)

