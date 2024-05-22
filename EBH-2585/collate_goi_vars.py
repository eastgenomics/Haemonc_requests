#!/usr/bin/env python

"""
Author: Jay Miles
Date: 22/05/2024

Purpose:
Given a list of genes of interest (specified in 'goi'), this script extracts
all variants in those genes from a defined set of annotated VCFs. It stores
variants in an intermediate pandas dataframe while parsing files, then writes
the completed dataframe to an .xlsx workbook ('output'). Variants in the
separate genes of interest are written to separate worksheets.

Prerequisites:
Uses the os, re, pandas and pysam modules. Also assumes that there exists a
subdirectory (specified in 'input_dir') containing at least one VCF file. VCF
files are expected to be named in the format
"{parent DNAnexus project}_{sample name}.vcf". An example would be:

parent project: 002_220606_A01303_0075_BH555NDRX2_MYE
sample name: 2204363-22129Z0001-1-FFPE-MCL-MYE-M-EGG2_S9_L001
VCF file name: 002_220606_A01303_0075_BH555NDRX2_MYE_204363-22129Z0001-1-FFPE-MCL-MYE-M-EGG2_S9_L001.vcf
"""


import os
import re
import pandas as pd
from pysam import VariantFile


# define output file, fields to populate it with, and genes of interest
input_dir = "all_filtered_vcfs"
output = "ebh-2585_variants.xlsx"

fields = ["project", "sample", "CHROM", "POS", "ID", "REF", "ALT", "QUAL",
    "FILTER", "INFO", "FORMAT", "format_values"]

goi = ["ANKRD26", "CEBPA", "DDX41", "ETV6", "GATA2", "RUNX1", "TP53"]

# define a dataframe to hold variants of interest
df=pd.DataFrame(columns=fields)

# iterate over downloaded VCFs
for vcf in os.listdir(input_dir):

    # get project and sample from file name
    try:
        project = re.findall(r"(002_\d{6}_(.*?)MYE)", vcf)[0][0]
        sample = re.findall(r"_MYE_((.*?))[\.]", vcf)[0][0]
    except IndexError:
        print(f"ERROR: Regex issue for {vcf}")

    # empty VCFs won't work, so skip them
    try:
        vcf_data = VariantFile(f"{input_dir}/{vcf}")
    except ValueError:
        print(f"ERROR: VariantFile parsing issue for {vcf}")
        continue

    # iterate over variant records
    for var in vcf_data.fetch():

        # some VCFs don't have an INFO/SYMBOL field
        try:
            symbol = var.info["SYMBOL"][0]
        except KeyError:
            symbol = ""

        # if variant is in a GOI, append to df
        if symbol in goi:
            new_row = [project, sample] + [f.strip() for f in str(var).split("\t")]
            df.loc[len(df)] = new_row

# write df to output xlsx - one sheet for each GOI
with pd.ExcelWriter(output) as writer:
    for gene in goi:

        gene_df = df[df['INFO'].str.contains(gene)]
        gene_df = gene_df.sort_values(by=['CHROM', 'POS', 'REF', 'ALT'], ignore_index=True)
        gene_df.to_excel(writer, sheet_name=gene, header=True, index=False)
