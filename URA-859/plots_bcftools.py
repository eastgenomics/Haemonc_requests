"""
Created on Mon Nov  4 12:48:05 2024

@author: arun
"""

import os
import re
import subprocess
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr


def list_samplenames(directory):
    """
    List sample names extracted from file names in the specified directory.

    This function scans the given directory for files and extracts a specific
    pattern from each file name. The pattern consists of the first two
    segments of the file name, which are expected to be separated by hyphens
    ('-').

    Args:
        directory (str):
            The path to the directory from which to list sample
            names.

    Returns:
        patterns (list):
            A list of sample name patterns extracted from the file names.
            If the directory does not exist or cannot be accessed, an
            appropriate error message is returned as a string.

    Raises:
        FileNotFoundError: If the specified directory does not exist.
        PermissionError: If permission is denied to access the specified
        directory.

    Example:
        directory_path = '/your/directory'
        print(list_samplenames(directory_path))
    """
    try:
        # Get a list of files and directories in the specified directory
        files = os.listdir(directory)

        # Filter out directories, only keep files
        files = [f for f in files if os.path.isfile(os.path.join(directory, f))]

        # Extract the desired pattern from each file name
        patterns = []
        for file in files:
            # Split the filename at the hyphens
            segments = file.split("-")
            # Check if there are at least two segments to avoid index errors
            if len(segments) >= 2:
                # Join the first two segments and append to the patterns list
                pattern = f"{segments[0]}-{segments[1]}"
                patterns.append(pattern)

        return patterns
    except FileNotFoundError:
        return f"The directory '{directory}' does not exist."
    except PermissionError:
        return f"Permission denied to access '{directory}'."


def common_elements(list1, list2):
    return list(set(list1) & set(list2))


def bcftools_isec(working_dir, sample_name, version, folder_a, folder_b, output_folder):
    """
    Perform bcftools isec and return output vcfs as dataframes

    Parameters
    ----------
    working_dir : str
        Directory where the folder_a and folder_b is located
    sample_name : str
        Sample name for .vcf filename regex pattern recognition
    version : str
        version of the vcf file for regex pattern recognition
    folder_a : str
        folder name where vcf_a files are located
    folder_b : TYPE
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
    pattern = f".*{sample_name}.+_{version}\.vcf(\.gz)?$"

    # Find the vcf files of interest
    file_a = next(
        (
            filename
            for filename in os.listdir(f"{working_dir}/{folder_a}")
            if re.search(pattern, filename)
        ),
        None,
    )
    file_b = next(
        (
            filename
            for filename in os.listdir(f"{working_dir}/{folder_b}")
            if re.search(pattern, filename)
        ),
        None,
    )

    if file_a is None or file_b is None:
        raise FileNotFoundError(f"No file was found with {sample_name}")

    # Set up relevant directories
    path_a = f"{folder_a}/{file_a}"
    path_b = f"{folder_b}/{file_b}"
    output_dir = os.path.join(output_folder, sample_name)

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


def transformed_vcf(vcf_dataframe):
    """
    Transforms a VCF DataFrame by adding new columns for site,
    allele frequency (AF), and depth (DP).

    This function performs the following transformations:

    1. Adds a "site" column combining "CHROM", "POS", "REF", and "ALT"
       columns in the format: 'CHROM:POS_REF/ALT'.

    2. Extracts the allele frequency (AF) from the "SAMPLE" column and adds
       it as a new "AF" column.

    3. Extracts the depth (DP) from the "INFO" column and adds it as a new
       "DP" column.

    Parameters
    ----------
    vcf_dataframe : pd.DataFrame
        The input VCF DataFrame containing columns: "CHROM", "POS", "REF",
        "ALT", "SAMPLE", and "INFO".

    Returns
    -------
    pd.DataFrame
        A transformed VCF DataFrame with added "site", "AF", and "DP" columns.
    """

    # Create an identical dataframe
    vcf_df = vcf_dataframe
    if not vcf_df.empty:
        # Create the Site column
        vcf_df["site"] = (
            vcf_df["CHROM"].astype(str)
            + ":"
            + vcf_df["POS"].astype(str)
            + "_"
            + vcf_df["REF"].astype(str)
            + "/"
            + vcf_df["ALT"].astype(str)
        )

        # Create the AF
        vcf_df["AF"] = vcf_df["SAMPLE"].str.split(":", expand=True)[2]
        vcf_df["AF"] = vcf_df["AF"].astype(float)

        # Create the DP
        vcf_df["DP"] = (
            vcf_df["INFO"]
            .str.split(";")
            .apply(lambda x: "".join((y for y in x if y.startswith("DP="))))
        )
        vcf_df["DP"] = vcf_df["DP"].str.replace("DP=", "")
        vcf_df["DP"] = vcf_df["DP"].astype(int)
    else:
        pass

    return vcf_df


def remove_duplicates_vcf_df(vcf_df):
    """

    Parameters
    ----------
    vcf_df : pd.Dataframe
        Dataframe output from the transformed_vcf() function

    Returns
    -------
    vcf_df : pd.Dataframe
        Dataframe with unique variants

    """
    # drop duplicate sites
    if not vcf_df.empty:
        vcf_df = vcf_df.drop_duplicates(subset=["site", "AF", "DP"])
    else:
        pass

    return vcf_df


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

    common_sites = common_elements(vcfa_df.site, vcfb_df.site)

    print(
        "Number of Chrom:Pos sites shared between two vcfs: " + str(len(common_sites))
    )

    # create dataframes containing common site rows
    vcfa_df_common = vcfa_df.loc[vcfa_df["site"].isin(common_sites)].copy()
    vcfb_df_common = vcfb_df.loc[vcfb_df["site"].isin(common_sites)].copy()

    # Rename the DP column to include their corresponding version
    vcfa_df_common.rename({"DP": f"DP_{version_a}"}, inplace=True)
    vcfb_df_common.rename({"DP": f"DP_{version_b}"}, inplace=True)

    # Rename the AF column to include their corresponding version
    vcfa_df_common.rename({"AF": f"AF_{version_a}"}, inplace=True)
    vcfb_df_common.rename({"AF": f"AF_{version_b}"}, inplace=True)

    # Merge both vcf files
    vcf_merged = pd.merge(
        vcfa_df_common,
        vcfb_df_common,
        on="site",
        suffixes=(f"_{version_a}", f"_{version_b}"),
    )

    return vcf_merged


def make_correlation_plot_from_merged_vcf(
    vcf_merged, x_col, y_col, title="", xlabel=None, ylabel=None
):
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

    Example
    -------
    make_correlation_plot_from_merged_vcf(vcf_merged, "AF_old_probes",
                                          "AF_new_probes", "JAK5A-JAK5A")

    """

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col

    # Calculate the R and p values
    r_value, p_value = pearsonr(vcf_merged[x_col], vcf_merged[y_col])
    # print(r_value, p_value)
    # Plot the correlation plot
    plot = sns.regplot(x=x_col, y=y_col, data=vcf_merged)
    plot.set_title(title)
    plot.set_xlabel(xlabel)
    plot.set_ylabel(ylabel)
    plot.text(
        x=min(vcf_merged[x_col]),
        y=max(vcf_merged[y_col]),
        s=f"R = {r_value:.2f}, p = {p_value:.2e}",
        ha="left",
        va="top",
    )
    plt.show()

    return plot


def venn_diagramify(vcfa_df, vcfb_df, vcf_merged, version_a, version_b, samplename):
    """
    Generate a venn diagram using the 2 transformed vcf dfs, return variants
    which are unique to to each vcf df.

    Parameters
    ----------
    vcfa_df : pd.Dataframe
        vcf dataframe which is output of transformed_vcf(vcf_dataframe)
    vcfb_df : pd. Dataframe
        another vcf dataframe which is output of transformed_vcf(vcf_dataframe)
    version_a : str
        distinct characteristic of vcfa_df
    version_b : str
        distinct characteristic of vcfb_df
    samplename : str
        Name of the sample where the vcf files originated

    Returns
    -------
    pd.Dataframe
        dataframe with unique variants for vcfb_df
    pd.Dataframe
        dataframe with unique variants for vcfb_df
    matplotlib_venn._common.VennDiagram
        Venn diagram plot

    """

    plot = venn2(
        subsets=(len(vcfa_df), len(vcfb_df), len(vcf_merged)),
        set_labels=(version_a, version_b),
    )

    plt.title(samplename, fontweight="bold", fontsize=20, pad=30)
    plt.show()

    return plot


def main():

    path = (
        "/home/arun/Codes/HaemOnc_panel_2024-Twist_capture/"
        + "Validation/plots_to_present/vcf_files"
    )

    os.chdir(path)
    cwd = os.getcwd()
    print(cwd)

    # set up the directories
    WORKING_DIRECTORY = path
    sp_directory = WORKING_DIRECTORY + "/SP"
    sp_old_probes_directory = sp_directory + "/old_probes_filtered"
    sp_new_probes_directory = sp_directory + "/new_probes_filtered"

    s1_directory = WORKING_DIRECTORY + "/S1"

    os.chdir(WORKING_DIRECTORY)
    cwd = os.getcwd()

    # Get the list of samples across different probes
    old_sample_list = list_samplenames(sp_old_probes_directory)
    new_sample_list = list_samplenames(sp_new_probes_directory)

    # Get the list of sample names that are common between old run and new run
    common_sample_list = common_elements(old_sample_list, new_sample_list)
    # print(len(common_sample_list))

    # prepare unique_variants_df
    unique_variants_df = pd.DataFrame()
    unique_variants_df_PASS = pd.DataFrame()

    for sample in common_sample_list:

        # Bcftools isec SP samples
        result = bcftools_isec(
            sp_directory,
            sample,
            "allgenesvep_filtered",
            "old_probes_filtered",
            "new_probes_filtered",
            "bcftools_isec",
        )

        sp_vcfold_df_uniq = result[0]
        sp_vcfnew_df_uniq = result[1]
        sp_vcfold_df_common = result[2]
        sp_vcfnew_df_common = result[3]

        # Bcftools isec S1 samples
        s1_result = bcftools_isec(
            s1_directory,
            sample,
            "allgenesvep_filtered",
            "old_probes_filtered",
            "new_probes_filtered",
            "bcftools_isec",
        )

        s1_vcfold_df_uniq = s1_result[0]
        s1_vcfnew_df_uniq = s1_result[1]
        s1_vcfold_df_common = s1_result[2]
        s1_vcfnew_df_common = s1_result[3]

        # transform dataframes for VAF and AF
        sp_vcfold_df_uniq = remove_duplicates_vcf_df(transformed_vcf(sp_vcfold_df_uniq))
        sp_vcfnew_df_uniq = remove_duplicates_vcf_df(transformed_vcf(sp_vcfnew_df_uniq))
        sp_vcfold_df_common = remove_duplicates_vcf_df(
            transformed_vcf(sp_vcfold_df_common)
        )
        sp_vcfnew_df_common = remove_duplicates_vcf_df(
            transformed_vcf(sp_vcfnew_df_common)
        )

        s1_vcfold_df_uniq = remove_duplicates_vcf_df(transformed_vcf(s1_vcfold_df_uniq))
        s1_vcfnew_df_uniq = remove_duplicates_vcf_df(transformed_vcf(s1_vcfnew_df_uniq))
        s1_vcfold_df_common = remove_duplicates_vcf_df(
            transformed_vcf(s1_vcfold_df_common)
        )
        s1_vcfnew_df_common = remove_duplicates_vcf_df(
            transformed_vcf(s1_vcfnew_df_common)
        )

        # Generate merged df
        sp_vcf_merged = merge_vcf_common_sites(
            sp_vcfold_df_common, sp_vcfnew_df_common, "old_probes", "new_probes"
        )
        s1_vcf_merged = merge_vcf_common_sites(
            s1_vcfold_df_common, s1_vcfnew_df_common, "old_probes", "new_probes"
        )

        # Print out the sample name
        print(sample)

        # Make AF correlation plot
        make_correlation_plot_from_merged_vcf(
            sp_vcf_merged, "AF_old_probes", "AF_new_probes", f"{sample} SP"
        )

        make_correlation_plot_from_merged_vcf(
            s1_vcf_merged, "AF_old_probes", "AF_new_probes", f"{sample} S1"
        )

        # Make DP correlation plots
        make_correlation_plot_from_merged_vcf(
            sp_vcf_merged, "DP_old_probes", "DP_new_probes", f"{sample} SP"
        )
        make_correlation_plot_from_merged_vcf(
            s1_vcf_merged, "DP_old_probes", "DP_new_probes", f"{sample} S1"
        )

        # Venn Diagrams
        venn_diagramify(
            sp_vcfold_df_uniq,
            sp_vcfnew_df_uniq,
            sp_vcf_merged,
            "old probes",
            "new probes",
            f"{sample} SP",
        )
        venn_diagramify(
            s1_vcfold_df_uniq,
            s1_vcfnew_df_uniq,
            s1_vcf_merged,
            "old probes",
            "new probes",
            f"{sample} S1",
        )

        # Bcftools isec with PASS variants
        result = bcftools_isec(
            sp_directory,
            sample,
            "allgenesvep_filtered",
            "old_probes_filtered_PASS",
            "new_probes_filtered_PASS",
            "bcftools_isec_PASS",
        )
        sp_oldpass_df_uniq = result[0]
        sp_newpass_df_uniq = result[1]
        sp_oldpass_df_common = result[2]
        sp_newpass_df_common = result[3]

        result = bcftools_isec(
            s1_directory,
            sample,
            "allgenesvep_filtered",
            "old_probes_filtered_PASS",
            "new_probes_filtered_PASS",
            "bcftools_isec_PASS",
        )
        s1_oldpass_df_uniq = result[0]
        s1_newpass_df_uniq = result[1]
        s1_oldpass_df_common = result[2]
        s1_newpass_df_common = result[3]

        # Venn Diagrams with only PASS variants
        venn_diagramify(
            sp_oldpass_df_uniq,
            sp_newpass_df_uniq,
            sp_newpass_df_common,
            "old probes",
            "new probes",
            f"{sample} PASS variants SP",
        )

        venn_diagramify(
            s1_oldpass_df_uniq,
            s1_newpass_df_uniq,
            s1_newpass_df_common,
            "old probes",
            "new probes",
            f"{sample} PASS variants S1",
        )

        # prepare unique variants dataframes
        if not sp_vcfold_df_uniq.empty:
            sp_vcfold_df_uniq.loc[:, "Flowcell"] = "SP"
            sp_vcfold_df_uniq.loc[:, "Probe set"] = "old probes"
            sp_vcfold_df_uniq.loc[:, "Samplename"] = sample

        if not sp_vcfnew_df_uniq.empty:
            sp_vcfnew_df_uniq.loc[:, "Flowcell"] = "SP"
            sp_vcfnew_df_uniq.loc[:, "Probe set"] = "new probes"
            sp_vcfnew_df_uniq.loc[:, "Samplename"] = sample

        if not s1_vcfold_df_uniq.empty:
            s1_vcfold_df_uniq.loc[:, "Flowcell"] = "S1"
            s1_vcfold_df_uniq.loc[:, "Probe set"] = "old probes"
            s1_vcfold_df_uniq.loc[:, "Samplename"] = sample

        if not s1_vcfnew_df_uniq.empty:
            s1_vcfnew_df_uniq.loc[:, "Flowcell"] = "S1"
            s1_vcfnew_df_uniq.loc[:, "Probe set"] = "new probes"
            s1_vcfnew_df_uniq.loc[:, "Samplename"] = sample

        # Concatenate all dataframes into 1
        dataframe_tuple = (
            sp_vcfold_df_uniq,
            sp_vcfnew_df_uniq,
            s1_vcfold_df_uniq,
            s1_vcfnew_df_uniq,
        )

        for dataframe in dataframe_tuple:
            unique_variants_df = pd.concat(
                [unique_variants_df, dataframe], axis=0, ignore_index=True
            )

        # prepare unique variants dataframes for PASS variants
        if not sp_oldpass_df_uniq.empty:
            sp_oldpass_df_uniq.loc[:, "Flowcell"] = "SP"
            sp_oldpass_df_uniq.loc[:, "Probe set"] = "old probes"
            sp_oldpass_df_uniq.loc[:, "Samplename"] = sample

        if not sp_newpass_df_uniq.empty:
            sp_newpass_df_uniq.loc[:, "Flowcell"] = "SP"
            sp_newpass_df_uniq.loc[:, "Probe set"] = "new probes"
            sp_newpass_df_uniq.loc[:, "Samplename"] = sample

        if not s1_oldpass_df_uniq.empty:
            s1_oldpass_df_uniq.loc[:, "Flowcell"] = "S1"
            s1_oldpass_df_uniq.loc[:, "Probe set"] = "old probes"
            s1_oldpass_df_uniq.loc[:, "Samplename"] = sample

        if not s1_newpass_df_uniq.empty:
            s1_newpass_df_uniq.loc[:, "Flowcell"] = "S1"
            s1_newpass_df_uniq.loc[:, "Probe set"] = "new probes"
            s1_newpass_df_uniq.loc[:, "Samplename"] = sample

        # Concatenate all dataframes into 1
        PASS_dataframe_tuple = (
            sp_oldpass_df_uniq,
            sp_newpass_df_uniq,
            s1_oldpass_df_uniq,
            s1_newpass_df_uniq,
        )

        for dataframe in PASS_dataframe_tuple:
            unique_variants_df_PASS = pd.concat(
                [unique_variants_df_PASS, dataframe], axis=0, ignore_index=True
            )

    unique_variants_df.to_csv("unique_variants_list.tsv",
                              sep="\t", index=False)
    unique_variants_df_PASS.to_csv("unique_variant_list_PASS.tsv",
                                   sep="\t", index=False)


if __name__ == "__main__":
    main()
