# germline_from_somatic
Identifying germline variants from somatic data

## File summary
### unarchive_and_download.sh

This script identifies all DNAnexus projects of the form 002_*_MYE created after 31/05/2022, and identifies all files within these projects of the form *_markdup_recalibrated_tnhaplotyper2_allgenesvep.vcf. It then iterates over these files to download them if they are live, or unarchive them if necessary.

### collate_goi_vars.py

This script iterates over the downloaded VCFs and extracts all variants in a specified set of genes. It stores variants in an intermediate pandas dataframe while parsing files, then writes the completed dataframe to an .xlsx workbook. Variants in the separate genes of interest are written to separate worksheets. The data written out for each variant consists of the sample and project of origin, plus its entire VCF record.
