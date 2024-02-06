#!/usr/bin/env python

import sys
import glob
import statistics

output_path = sys.argv[1]

# define VCFs to look at
small_var_vcfs = glob.glob(f'{output_path}*_intersect_small_vars.vcf')
indel_vcfs = glob.glob(f'{output_path}*_intersect_indels.vcf')

# look at different VCF types separately
for group in ("small variant", small_var_vcfs), ("indel", indel_vcfs):

    counts = []

    # get number of non-blank lines (i.e. variants) in each VCF
    for vcf in group[1]:
        with open(vcf, 'r') as reader:
            lines = [line for line in reader.readlines() if line.strip()]

        counts.append(len(lines))

    # get descriptive statistics
    counts_min = min(counts)
    counts_max = max(counts)
    counts_mean = sum(counts) / len(counts)
    counts_median = statistics.median(counts)

    # print results
    print(f"Numbers of variants in {group[0]} VCFs:\n"
        f"Mean: {counts_mean}\n"
        f"Median: {counts_median}\n"
        f"Maximum: {counts_max}\n"
        f"Minimum: {counts_min}\n")
