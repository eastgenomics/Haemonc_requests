#!/bin/bash/

# First the file is downloaded for variant filtering before they are compared
rm -r sentieon-tnbam_comparison; mkdir sentieon-tnbam_comparison

for vcf_file in sentieon-tnbam_v3.2.0/*.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    bcftools norm -a "$vcf_file" | 
    bcftools filter -i "FORMAT/DP > 99 && FORMAT/AF > 0.03 && FILTER='PASS'" |
    bedtools intersect -header -u -a "$vcf_file" -b coding_unrestricted_GRCh38_myeloid_5bp_flank_v2.1.0.bed > sentieon-tnbam_comparison/"${sample}_filtered_v3.2.0.vcf"
    done

for vcf_file in sentieon-tnbam_v5.0.1/*.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    bcftools norm -a "$vcf_file" |
    bcftools filter -i "FORMAT/DP > 99 && FORMAT/AF > 0.03 && FILTER='PASS'" |
    bedtools intersect -header -u -a "$vcf_file" -b coding_unrestricted_GRCh38_myeloid_5bp_flank_v2.1.0.bed > sentieon-tnbam_comparison/"${sample}_filtered_v5.0.1.vcf"
    done

for i in sentieon-tnbam_comparison/*.vcf; do bgzip $i; done
for i in sentieon-tnbam_comparison/*.vcf.gz; do tabix -p vcf $i; done

rm -r bcftools_isec ; mkdir bcftools_isec

# Create stuff for venn diagrams
for vcf_file in sentieon-tnbam_v3.2.0/*.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    bcftools isec -p bcftools_isec/${sample} sentieon-tnbam_comparison/${sample}_filtered_v5.0.1.vcf.gz sentieon-tnbam_comparison/${sample}_filtered_v3.2.0.vcf.gz
    done


