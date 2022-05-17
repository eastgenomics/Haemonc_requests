###########################   Load packages   ###########################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

###########################   Read data in   ###########################

option_list <- list(
  make_option(c("-c", "--cosmic_data"), action="store", default=NA,type='character',
              help="COSMIC Mutation Data"),
  make_option(c("-i", "--input_gene"), action="store", default=NA,type='character',
              help="A text file with row of genes"),
  make_option(c("-o", "--output_name"), action="store", default="COSMIC_output.tsv",type='character',
              help="String of the output filename")
)

opt = parse_args(OptionParser(option_list=option_list))

gene_list <- read.delim(opt$input_gene, sep = "\t", stringsAsFactors = F, header = F)
gene_list <- gene_list[,1] # this is a dataframe so select the first (and only) column

dat <- read_csv(opt$cosmic_data) %>%
   separate(GENE_NAME,  c("Gene", "Transcript"), "_")

###########################   Main pipeline   ###########################

summary_tbl <- dat %>%
  # filtering
  filter(Gene %in% gene_list,
         str_detect(LEGACY_MUTATION_ID, "COSM"),
         GENOMIC_MUTATION_ID != "null",
         MUTATION_DESCRIPTION != "Unknown") %>%

  group_by(Gene) %>%
  # calculates unique cosmic IDs without getting dups
  mutate(Total = n_distinct(LEGACY_MUTATION_ID))  %>%
  select(Gene, ACCESSION_NUMBER, LEGACY_MUTATION_ID,Total) %>%
  # Calculate unique legacy_mutation_IDs (COSMs) per transcripts
  group_by(ACCESSION_NUMBER) %>%
  summarise(n_COSMs = n_distinct(LEGACY_MUTATION_ID), across()) %>%
  # order columns
  select(Gene,ACCESSION_NUMBER, n_COSMs, Total) %>%
  unique() %>%
  # order rows by gene and transcript
  arrange(Gene, ACCESSION_NUMBER) %>%
  # rename some headers to be more clear
  rename(Transcript= ACCESSION_NUMBER,
         total_COSMs= Total)

# print some on the console for quick check
summary_tbl

###########################   Save output   ###########################
write.table(summary_tbl, opt$o, sep = "\t",row.names=FALSE)
