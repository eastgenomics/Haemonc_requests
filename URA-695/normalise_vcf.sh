## normalise the VCFs to deal with multiallelic variants
for vcf_file in sentieon-tnbam_v3.2.0/*.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    bcftools norm -f GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18_noChr.fa.gz -m -any --keep-sum AD  "$vcf_file"  | bcftools view -i 'ALT!="*"' |  bgzip  > sentieon-tnbam_v3.2.0/"${sample}_normalised.vcf.gz"
done

for vcf_file in sentieon-tnbam_v5.0.1/*.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    bcftools norm -f GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18_noChr.fa.gz -m -any --keep-sum AD  "$vcf_file" | bcftools view -i 'ALT!="*"' | bgzip  > sentieon-tnbam_v5.0.1/"${sample}_normalised.vcf.gz"
done