setwd("~/Documents/Projects/myeloid/haemonc_smalltasks/bin/Haemonc_requests/EBH-929/data_vcf")

library(tidyverse)
library(ggplot2)
library(ggpubr)
mydir = "~/Documents/Projects/myeloid/haemonc_smalltasks/bin/Haemonc_requests/EBH-929/data_vcf"
total_files =  list.files(path = paste0(mydir))



for(sample in 1:length(total_files)) {
  print(sample)
  sample_folder = list.files(path = paste0(mydir))[sample]
  # vcf_1 is the old capture
  # vcf_3 is the new capture
  vcf_1 = readr::read_delim(paste(sample_folder, "/002_norm.vcf" , sep = ""), comment = '#', col_names = F, delim = "\t")
  colnames(vcf_1) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  vcf_1$chr_pos_ref <- paste(vcf_1$CHROM, vcf_1$POS, vcf_1$REF, sep = "_")
  vcf_1$VAF = as.numeric(str_split_fixed(vcf_1$SAMPLE, ":", 4)[,3])
  vcf_3 = readr::read_delim(paste(sample_folder, "/003_norm.vcf" , sep = ""), comment = '#', col_names = F, delim = "\t")
  colnames(vcf_3) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  vcf_3$chr_pos_ref <- paste(vcf_3$CHROM, vcf_3$POS, vcf_3$REF, sep = "_")
  vcf_3$VAF = as.numeric(str_split_fixed(vcf_3$SAMPLE, ":", 4)[,3])
  
  
  assign(paste(sample_folder, "_0002" , sep = ""), vcf_1)
  assign(paste(sample_folder, "_0003" , sep = ""), vcf_3)

}

rm(vcf_1, vcf_3)

pdf("~/Documents/Projects/myeloid/haemonc_smalltasks/bin/Haemonc_requests/EBH-929/Shared_VAF_distributions.pdf")
for(sample in total_files) {
  print(sample)
  vcf_1 = get(paste(sample,"_0002", sep = ""))
  vcf_3 = get(paste(sample,"_0003", sep = ""))
  print(identical(vcf_1$chr_pos_ref, vcf_3$chr_pos_ref))
  vcf_all <- as.data.frame(cbind(vcf_1$VAF, vcf_3$VAF, vcf_1$chr_pos_ref))
  colnames(vcf_all) <- c("Old_capture", "New_capture", "Site")
  vcf_all$Old_capture <- as.numeric(vcf_all$Old_capture)
  vcf_all$New_capture <- as.numeric(vcf_all$New_capture)
  
  p = ggscatter(vcf_all, x = "Old_capture", y = "New_capture",
            add = "reg.line", conf.int = TRUE, title = sample )
  p2 = p + stat_cor(method = "pearson")
  print(p2)
}
dev.off()

library(plotly)
p = plot_ly(vcf_all, x = vcf_all$Old_capture, y = vcf_all$New_capture, text = paste("Site: ", vcf_all$Site),
        type = "scatter",mode = "markers")
