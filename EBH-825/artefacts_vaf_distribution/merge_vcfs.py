import pandas as pd
import glob
import os
import gzip

# this contains novaseq samples, with less than 1% contamination
# these are the vcfs normally used in the MAF file
os.chdir('../../../../data/merged_vcf/vcfs/')
print(os.getcwd())
all_vcfs = glob.glob("*.vcf.gz")

cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
all_sample_dict = []

for sample_vcf in all_vcfs:
  # read in sample
  vcf_df = pd.read_csv(sample_vcf, sep="\t",
                    names=cols, compression='infer', comment = "#")
  # Select rows where it starts with chr (removes the headers)
  header = []
  with gzip.open(sample_vcf) as fh:
    for line in fh.readlines():
        if sample_vcf.endswith('.gz'):
            line = line.decode()
        if line.startswith('#'):
            header.append(line.rstrip('\n'))
        else:
            break

  # find what info headers shoud be from the vcf headers CSQ line
  csq = list(filter(lambda x:'CSQ' in x, header))[0]
  csq_headers = csq.split("Format: ")[1].split(">")[0]
  info_colnames = csq_headers.split("|")
  # last element in list has a " attached so need to remove it
  info_colnames[-1] = info_colnames[-1].replace('"','')
  vcf_df[info_colnames] = vcf_df['INFO'].str.split('CSQ=', 1, expand=True)[1].str.split('|',expand=True)
  # loop through the variants  in the vcf to extract relevant info
  for idx, row in vcf_df.iterrows():
    sample_dict = {}
    sample_dict["sample_name"] = sample_vcf.split('_')[0]
    sample_dict["chr"] = vcf_df.at[idx, "CHROM"]
    sample_dict["pos"] = vcf_df.at[idx, "POS"]
    sample_dict["ref"] = vcf_df.at[idx, "REF"]
    sample_dict["alt"] = vcf_df.at[idx, "ALT"]
    sample_dict["vaf"] = vcf_df.at[idx, "SAMPLE"].split(':')[2]
    # all annotations are the same per vcf so we can a  pply the info_col
    sample_dict["gene"] = vcf_df.at[idx, "SYMBOL"]
    sample_dict["HGVSc"] = vcf_df.at[idx, "HGVSc"]
    sample_dict["transcript"] = vcf_df.at[idx, "HGVSc"].split(':')[0]
    sample_dict["prev_count"] = vcf_df.at[idx, "Prev_AC"]
    sample_dict["total_count"] = vcf_df.at[idx, "Prev_NS"]
    all_sample_dict.append(sample_dict)


df = pd.DataFrame.from_dict(all_sample_dict)
df['chrom_ref'] = df['chr'] + '_' + df['ref']

df.prev_count = pd.to_numeric(df.prev_count)
df['prev_count'] = df['prev_count'].fillna(0).astype(int)

# select more than 20 to get common variants
filtered_df = df[df.prev_count > 20]
df.to_csv(os.path.join(os.path.realpath('../../../'), 'results/EBH-825/220201_allvcf_merged_2.tsv'), sep="\t", header=True, index=False)
filtered_df.to_csv(os.path.join(os.path.realpath('../../../'), 'results/EBH-825/filtered_10_perc_all_variants_2.tsv'), sep="\t", header=True, index=False)



