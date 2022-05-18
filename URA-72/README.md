# Cosmic HaemOnc analysis README


## Introduction
The aim of this analysis is to see whether the mane transcripts retain more of the variants in cosmic than the alternative transcripts of a gene. Script COSMIC_HaemOnc_new_genes.R is used to generate the results table. The genes are provided from a text file with the gene names listed. To trigger the running of the script 
```Rscript COSMIC_HaemOnc_new_genes.R -c annotation_sources/cosmic/V95_38_MUTANT.csv -i HaemOnc_genes_v2.txt -o HaemOnc_v2_COSMs.tsv```
where the out.txt is the name of the desired output file. If this is missing it will use the input filename with results prefix.

## Method

From COSMIC, the file from "COSMIC Mutation Data" (https://cancer.sanger.ac.uk/cosmic/download) was filtered for cancer "haematopoietic_and_lymphoid_tissue" and downloaded. 

There are number of ways to track a COSMIC variants, either genomic position, Genomic mutation identifier (COSV) or Legacy mutation identifier (COSM).
For this analysis, the Genomic mutation identifier (COSV) will filtered for anything not null then the Legacy mutation identifier (COSM) will be filtered for coding ID (there are COSN which is non coding). Then cosmic variants with no protein mutation description are removed as some COSM are outside exons but labelled COSM as they are thought to impact other transcript functions such as splicing, promoter regions etc.

COSMIC is based off ENSEMBL transcripts and not all ensembl transcripts can be easily mapped to RefSeq transcripts which is what we use. So this analysis will only be looking at the ENSEMBL but their respective RefSeq transcripts will be noted to make this analysis more clear for us. The Refseq and ENSEMBL transcripts were mapped using CCDS (e.g https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=ALLFIELDS&DATA=CXCR4&ORGANISM=0&BUILDS=CURRENTBUILDS).
