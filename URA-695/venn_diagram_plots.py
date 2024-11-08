#do venn diagram with number of variants 
#reported in each output from bcftools isec
#single clinical sample mutect2 (bam to vcf)
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import io
import os
#function to read vcfs that were output from bcftools isec
# Things to do to polish the code
# make it save a file for each sample and it its appropriate directory

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
#read the vcfs to use for the venn diagram of variants
list_of_samples = os.listdir('bcftools_isec/')
list_of_samples.sort()

img_filenames = []
for sample in list_of_samples:
    prod_vcf = read_vcf(f"bcftools_isec/{sample}/0000.vcf")
    dev_vcf = read_vcf(f"bcftools_isec/{sample}/0001.vcf")
    dev_shared = read_vcf(f"bcftools_isec/{sample}/0003.vcf")
    #get the number of variants in each vcf
    #the column POS was chosen as an example as all have the same length
    set1_size = len(prod_vcf['POS'])
    set2_size = len(dev_vcf['POS'])
    overlap_size = len(dev_shared['POS'])
    # Create a Venn diagram
    plot = venn2(subsets=(set1_size, set2_size, overlap_size), set_labels=('Sentieon_v3', 'Sentieon_v5'))
    samplename = sample.split("_")[0]
    plt.title(f"{samplename}")
    plt.savefig(f"bcftools_isec/{sample}/{sample}_diagram.jpeg")
    img_filenames.append(f"bcftools_isec/{sample}/{sample}_diagram.jpeg")
    plt.close()
