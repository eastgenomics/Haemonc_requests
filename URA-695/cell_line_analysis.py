from cProfile import label
import pandas as pd
from scipy.stats import pearsonr
from correlation_plots import transformed_vcf, merge_vcf_common_sites
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from plotnine import *


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
                            names=cols, low_memory=False)

    return df

def read_cell_line_truths(sample):
    """
    Read in cell lines that are bed files


    Returns
    -------
    pd.DataFrame
        Merged dataframe with the given 2 vcf files

    """

    cell_line = pd.read_csv(f'cell_lines/{sample}_variants.bed', sep= "\t")
    cell_line = cell_line[["GRCh38 co-ordinates", "NGS allelic frequency (%)"]]
    cell_line.columns =['site', 'expected_af']
    cell_line = cell_line[cell_line.expected_af.str.contains("%")]
    cell_line['expected_af'] = cell_line['expected_af'].str.split('%', expand=True)[0]
    cell_line["expected_af"] = pd.to_numeric(cell_line["expected_af"])
    cell_line["expected_af"] = cell_line["expected_af"]

    return cell_line

def correlation_plot(query_vcf, truth_df,query_col, truth_col, title='',
                                          xlabel=None, ylabel=None):

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = truth_col
    if ylabel is None:
        ylabel = query_col


    # merge the query and truth on shared sites
    merged_dat = pd.merge(query_vcf, truth_df, on='site')
    # Calculate the R and p values
    r_value, p_value = pearsonr(merged_dat[query_col], merged_dat[truth_col])
    print(r_value, p_value)
    n_value = len(merged_dat.index)
    # Plot the correlation plot
    plot = sns.regplot(x=truth_col, y=query_col, data=merged_dat)
    plot.set_title(title)
    plot.set_xlabel(xlabel)
    plot.set_ylabel(ylabel)
    # plot.axis([0, 100, 0, 100])
    plot.text(x=min(merged_dat[truth_col]),
              y=max(merged_dat[query_col]),
              s=f'R = {r_value:.2f}, p = {p_value:.2e}, n = {n_value}',
              ha='left', va='top')
    title_file = "plots/" + title + "_" + str(query_col) + "_" + str(truth_col) + ".jpg"
    plt.savefig(title_file, dpi=300)
    plt.close()

def all_boxplot(query_vcf, truth_df, query_col, query_col2, truth_col, title='',
                                          xlabel=None, ylabel=None):

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = query_col
    if ylabel is None:
        ylabel = truth_col

    # merge the query and truth on shared sites
    merged_dat = pd.merge(query_vcf, truth_df, on='site')
    data = [merged_dat[query_col].tolist(), merged_dat[query_col2].tolist(), merged_dat[truth_col].tolist()]
    title_file = "plots/" + title + "_boxplots_" + str(query_col) + "_" + str(query_col2) + "_" + str(truth_col) + ".jpg"

    # Create a figure and axis
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

    for ax, group_data, title in zip(axes, data, [query_col, query_col2, truth_col]):
        sns.boxplot(data=group_data, ax=ax, color='lightblue')
        sns.stripplot(data=group_data, ax=ax, color='black', alpha=0.5, jitter=True)
        ax.set_title(title)
        ax.set_ylabel('Allele Frequency')
        ax.set_xticks([]) 
        ax.text(0.1, 0.95, 'n: {}'.format(len(group_data)),
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=12, bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()
    plt.savefig(title_file, dpi=300)
    plt.close()

def density_plot(query_vcf, truth_df, query_col, truth_col, title='',
                                          xlabel=None, ylabel=None):

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = query_col
    if ylabel is None:
        ylabel = truth_col

    # merge the query and truth on shared sites
    merged_dat = pd.merge(query_vcf, truth_df, on='site')
    plot_title =  title + "_density"

    # Create a figure and axis
    p = (
    ggplot(merged_dat)
    + geom_density(aes(x=query_col), color="blue")
    + geom_density(aes(x=truth_col), color="red")
    + ggtitle(plot_title)
    )

    p.save("plots/" + title + "_density_" + str(query_col) +  "_" + str(truth_col) + ".jpg")

def differences_barplot(query_vcf, truth_df, query_col, query_col2, truth_col, title='',
                                          xlabel=None, ylabel=None):

    # Set x and y labels if not provided
    if xlabel is None:
        xlabel = query_col
    if ylabel is None:
        ylabel = truth_col

    # merge the query and truth on shared sites
    merged_dat = pd.merge(query_vcf, truth_df, on='site')
    merged_dat[f'diff_{truth_col}_{query_col}']  = abs(merged_dat[truth_col] - merged_dat[query_col])
    merged_dat[f'diff_{truth_col}_{query_col2}']  = abs(merged_dat[truth_col] - merged_dat[query_col2])

    merged_dat.to_csv(f'plots/{title}_observed_expected_merged.tsv', index=False, sep = "\t")
    # Set up the bar widths and positions
    bar_width = 0.4
    y_positions = np.arange(len(merged_dat))

    # Create the horizontal bar plot
    plt.figure(figsize=(20,8))
    plt.barh(y_positions - bar_width/2, merged_dat[f'diff_{truth_col}_{query_col}'].tolist(), height=bar_width, label=f'diff_{truth_col}_{query_col}', color='lightblue')
    plt.barh(y_positions + bar_width/2, merged_dat[f'diff_{truth_col}_{query_col2}'].tolist(), height=bar_width, label=f'diff_{truth_col}_{query_col2}', color='salmon')

    # Set the y-ticks and labels
    plt.yticks(y_positions, merged_dat['site'], fontsize = 10)

    # Add a title and labels
    plt.xlabel('Absolute AF difference (%)')
    plt.title(f'Absolute difference between expected AF and observed AF of {title}')
    plt.legend()

    # Display the plot
    plt.tight_layout()
    title_file = "plots/" + title + "_barplot_difference_AF_" + str(truth_col) + ".jpg"
    plt.savefig(title_file, dpi=300)
    plt.close()


def main():
    #########################  Read in cell lines  #########################

    ### oncospan
    oncospan_xlsx = open('cell_lines/hd827-hd832-hd833-high-confidence-variants.xlsx', 'rb')
    oncospan = pd.read_excel(oncospan_xlsx,sheet_name='Sheet1', skiprows=2)
    oncospan = oncospan.iloc[:,[2,3,10]]
    oncospan.columns = ['chr', 'pos', 'expected_af']
    oncospan['site'] = oncospan['chr'].astype(str) + ':' + oncospan['pos'].astype(str)
    oncospan = oncospan[['site', 'expected_af']]
    oncospan["expected_af"] = pd.to_numeric(oncospan["expected_af"])
    oncospan["expected_af"] = oncospan["expected_af"] * 100

    ### HD734
    hd734 = read_cell_line_truths('HD734')

    ### HD829
    hd829 = read_cell_line_truths('HD829')

    ######################  Read in cell lines VCFs  ######################

    samples = {"TMv2OCT20-HD734-1st_S35_L001" : hd734,
                "TMv2OCT20-HD734-2nd_S36_L001": hd734,
                "TMv2OCT20-HD829-1st_S33_L001": hd829,
                "TMv2OCT20-HD829-2nd_S34_L001": hd829,
                "TMv2OCT20-oncospan-cell-line-1st_S39_L001": oncospan,
                "TMv2OCT20-oncospan-cell-line-2nd_S40_L001": oncospan}
    cwd = os.getcwd()
    for sample, truth in samples.items():
        print(sample)
        vcf3_df = give_vcf_df(sample, "v3.2.0", cwd)
        vcf5_df = give_vcf_df(sample, "v5.0.1", cwd)

        vcf3_df = transformed_vcf(vcf3_df)
        vcf5_df = transformed_vcf(vcf5_df)
        vcf3_df["DP"] = pd.to_numeric(vcf3_df["DP"])
        vcf3_df = vcf3_df[vcf3_df['DP'] > 99]
        vcf5_df["DP"] = pd.to_numeric(vcf5_df["DP"])
        vcf5_df = vcf5_df[vcf5_df['DP'] > 99] 

        vcf_merged = merge_vcf_common_sites(vcf3_df, vcf5_df, 'v3', 'v5', sample)
        # drop low DP samples
        vcf_merged = vcf_merged[~vcf_merged.FILTER_x.str.contains("multiallelic")]
        vcf_merged['site'] = vcf_merged['site'].str.split('_', expand=True)[0]
        vcf_merged['AF_v3'] = vcf_merged['AF_v3'] * 100
        vcf_merged['AF_v5'] = vcf_merged['AF_v5'] * 100

        print("Plotting AF diff to expected plots")
        differences_barplot(vcf_merged, truth, 'AF_v3', 'AF_v5', 'expected_af', sample)

        print("Plotting boxplots")
        all_boxplot(vcf_merged, truth, 'AF_v3', 'AF_v5', 'expected_af', sample)

        print("Plotting density plots")
        density_plot(vcf_merged, truth, 'AF_v3', 'expected_af', sample)
        density_plot(vcf_merged, truth, 'AF_v5', 'expected_af', sample)

        print("Plotting correlation plots")
        correlation_plot(vcf_merged, truth, 'AF_v3', 'expected_af', sample)
        correlation_plot(vcf_merged, truth, 'AF_v5', 'expected_af', sample)





if __name__ == "__main__":
    main()