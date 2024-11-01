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

os.chdir("/home/arun/Codes/sentieon_update/sentieon-tnbam_comparison")
cwd = os.getcwd()
print(cwd)

# De file the general pattern to match
pattern = '.+oncospan-cell-line-1st.+_v3.2.0\.vcf.gz$'

        
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
    
    Raises
    ------
    FileNotFoundError
        If no VCF file matching the given sample name and version is found in the specified directory.

    Examples
    --------
    >>> df = give_vcf_df('sample123', 'v1.2.3', '/path/to/vcfs')

    """
    # Construct the regular expression pattern to match the VCF file
    pattern = f'.+{sample_name}.+_{version}\.vcf.gz$'
    
    # Search for the first file in the directory matching the pattern.
    # If no file matches, return None instead of raising StopIteration.
    file = next((filename for filename in os.listdir(cwd) 
                 if re.search(pattern, filename)), None)
    
    # Raise an exception if no matching file is found
    if file is None:
        raise FileNotFoundError(f"No file was found with {sample_name} and {version}")

    # Define the columns for the VCF file DataFrame
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
            "SAMPLE"]
    # Load the VCF file into a pandas DataFrame
    dataframe = pd.read_csv(os.path.join(cwd, file), sep= "\t", comment='#',
                            names=cols)
    return dataframe

vcf3_df = give_vcf_df("oncospan-cell-line-1st", "v3.2.0", cwd)
vcf5_df = give_vcf_df("oncospan-cell-line-1st", "v5.0.1", cwd)

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

vcf3_df = transformed_vcf(vcf3_df)
vcf5_df = transformed_vcf(vcf5_df)


# This code bit is to figure out if there are multiallelic entries
#vcf_test_df['site'] = vcf_test_df['CHROM'].astype(str) + ':' + vcf_test_df['POS'].astype(str)
#multiallelic_vcf = vcf_test_df[vcf_test_df['ALT'].str.contains(",")]

def merge_vcf_common_sites(vcfa_df, vcfb_df, version_a, version_b):
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
    
    # Merge both vcf files
    vcf_merged = pd.merge(vcfa_df_common, vcfb_df_common, on='site')
    
    # Store DP values as integers
    vcf_merged[f'DP_{version_a}'] = vcf_merged[f'DP_{version_a}'].astype(int)
    vcf_merged[f'DP_{version_b}'] = vcf_merged[f'DP_{version_b}'].astype(int)
    
    return vcf_merged

vcf_merged = merge_vcf_common_sites(vcf3_df, vcf5_df, 'v3', 'v5')

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
    # Plot the correlation plot
    plot = sns.regplot(x=x_col, y=y_col, data=vcf_merged)
    plot.set_title(title)
    plot.set_xlabel(xlabel)
    plot.set_ylabel(ylabel)
    plot.text(x=min(vcf_merged[x_col]),
              y=max(vcf_merged[y_col]),
              s=f'R = {r_value:.2f}, p = {p_value:.2e}', 
              ha='left', va='top')
    return plot

plot = make_correlation_plot_from_merged_vcf(vcf_merged, "DP_v3", "DP_v5", "Oncospan")



#print("Number of Chrom:Pos sites shared between two vcfs: " + str(len(common_sites)))

#print("number of rows for vcf v2: " + str(len(vcf_test_df)))

