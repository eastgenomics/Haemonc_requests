## DI-1801 Generate counts from Genie data
This folder contains scripts to take Genie MAF data, aggregate counts of each variant,
convert each variant to a VCF-like description and output to a CSV and VCF.

### Subset Genie data
The original Genie data has millions of rows, therefore we are going to subset the data to only variants in genies that are present in the Uranus pipeline for simpler processing. This requires a haemonc BED file, an email to identify yourself with Entrez and any extra symbols to include (optional).
Example command:
```
python subset_genie_data.py \
  --input_maf data_mutations_extended.txt \
  --bed_file GCF_000001405.39_GRCh38.p13_genomic_20200815.exons_20bp_5bp_uranus_panel_v1.0.0.bed \
  --entrez_email my.email@nhs.net \
  --extra_symbols GNAS-AS1 \
  --output data_mutations_extended_haemonc_genes.txt
```

### Merge sample info
We need to merge in the clinical data (patient IDs, cancer types etc.) from Genie for each sample so that we can generate
count information.
Example command:
```
python merge_sample_info.py \
  --input_maf data_mutations_extended_haemonc_genes.txt \
  --clinical_info data_clinical_sample.txt \
  --output data_mutations_extended_haemonc_genes_clinical_info.txt
```

### Generate count data
For each unique variant, generate counts for all cancers, all haemonc cancers and each
haemonc cancer type, dependent on the input file of cancer types we're interested in.
Example command:
```
python generate_count_data.py \
  --input data_mutations_extended_haemonc_genes_clinical_info.txt \
  --haemonc_cancer_types haemonc_cancer_types.txt \
  --output genie_17_aggregated_counts_GRCh37.txt
```

### Convert counts to VCF
Because each variant is in MAF-like format, in order to represent Genie variants in the same way
as variants would be represented in a Uranus sample VCF, we want to convert them to a VCF-like
description. This also means we can perform liftover to GRCh38 properly.
A JSON file representing details of the fields to keep as INFO fields in the output VCF is required;
an example `info_fields.json` is provided, and a FASTA (sourced from Ensembl) is also required for the GRCh37 reference genome.
Example command:
```
python convert_counts_to_vcf.py \
  --input genie_17_aggregated_counts_GRCh37.txt \
  --fasta Homo_sapiens.GRCh37.dna.toplevel.fa.gz \
  --info_fields info_fields.json \
  --output_csv genie_17_aggregated_counts_GRCh37_vcf_description.csv \
  --output_vcf genie_17_aggregated_counts_GRCh37.vcf \
```

Note: the VCF should then be sorted, compressed and normalised with bcftools and lifted over to GRCh38 with Picard LiftoverVcf.

### Add GRCh38 liftover to counts CSV
To make a final CSV with GRCh38 lifted over variants and their aggregate counts, we want to read the CSV with count data in GRCh37 and the lifted over VCF back into a dataframe and merge the GRCh38 variant info with the CSV to write a final CSV file.
Example command:
```
python add_grch38_liftover.py \
  --csv genie_17_aggregated_counts_GRCh37_vcf_description.csv \
  --vcf Genie_v17_aggregated_GRCh38_v1.0.0.vcf.gz \
  --output Genie_v17_aggregated_GRCh38_v1.0.0.csv
```

### Create HTML table file
From the final CSV we want to generate a searchable HTML table with Tabulator.
A JSON file representing the columns to keep in the HTML file from the Genie count data CSV is required; an example `columns.json` file is provided. We also provide the name of the Jinja HTML template (`jinja_template.html`) in the /templates directory which the CSV data will be inserted into.
Example command:
```
python create_html.py \
    --input Genie_v17_aggregated_GRCh38_v1.0.0.csv \
    --columns columns.json \
    --template jinja_template.html \
    --output Genie_v17_aggregated_GRCh38_v1.0.0.html
```
