
for vcf_file in sentieon-tnbam_v3.2.0/*markdup_recalibrated_tnhaplotyper2.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    echo $sample
    #icdiff -H --unified=0 <(zcat sentieon-tnbam_v3.2.0/"${sample}.vcf.gz" | grep -v ^"#" | cut -f 1,2,4,5,10 ) <(zcat sentieon-tnbam_v5.0.1/"${sample}.vcf.gz" | grep -v ^"#" | cut -f 1,2,4,5,10 ) | aha  > sentieon-tnbam_comparison/${sample}__key_cols_vcf_diff.html
    #icdiff -H --unified=0 <(zcat sentieon-tnbam_v3.2.0/"${sample}.vcf.gz") <(zcat sentieon-tnbam_v5.0.1/"${sample}.vcf.gz" ) | aha  > sentieon-tnbam_comparison/${sample}_vcf_diff.html
    done

for vcf_file in sentieon-tnbam_v3.2.0/*markdup_recalibrated_tnhaplotyper2.vcf.gz
    do sample=$(basename "$vcf_file" .vcf.gz)
    uniq_v3=$(zgrep -v ^"#" sentieon-tnbam_v3.2.0/"${sample}.vcf.gz" | cut -f 7 | sort | uniq -c)
    uniq_v5=$(zgrep -v ^"#" sentieon-tnbam_v5.0.1/"${sample}.vcf.gz" | cut -f 7 | sort | uniq -c)
    printf "\n%s\t%s\t\n%s\t%s\t%s\t%s\n" "$sample" "sentieon_v3" "$uniq_v3" "sentieon_v5" "$uniq_v5" >> filter_change.tsv
    done