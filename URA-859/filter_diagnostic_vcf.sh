#!/bin/bash

cd "/home/arun/Codes/HaemOnc_panel_2024-Twist_capture/Validation/plots_to_present/vcf_files/Diagnostic_runs" || return

FILTER_CRITERIA='FORMAT/DP > 250 && FORMAT/AF > 0.05'

for dir in */; do
    # Ensure we are working with directories
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"

        # Find all VCF files in the directory
        for vcf_file in "$dir"*allgenesvep.vcf; do
            # Check if VCF files exist
            if [ -f "$vcf_file" ]; then

                echo "  Checking file: $vcf_file"

                # Define a temporary file for chr modification
                tmp_vcf="${vcf_file%.vcf}_nochr.vcf"
                
                # Check if the file contains "chr" in the chromosome column and remove it
                if grep -q "^chr" "$vcf_file"; then
                    echo "  'chr' found in $vcf_file. Removing prefixes..."
                    # Use sed to remove 'chr' prefixes
                    sed 's/^chr//' "$vcf_file" > "$tmp_vcf"
                    mv "$tmp_vcf" "$vcf_file"
                    echo "  Prefixes removed."
                else
                    echo "  No 'chr' prefix found. Skipping modification."
                fi

                echo "  Filtering file: $vcf_file"
                
                # Define output file name
                output_file="${vcf_file%.vcf}_filtered.vcf"
                
                # Run bcftools filter
                bcftools filter -i "$FILTER_CRITERIA" -o "$output_file" -Ov "$vcf_file"
                bgzip -f "$output_file"
                tabix -p vcf "${output_file}.gz"

                echo "  Output written to: $output_file"
            else
                echo "  No VCF files found in $dir"
            fi
        done
    fi
done

#bcftools filter -i 'FORMAT/DP > 250 && FORMAT/AF > 0.05' -Oz -o ${sample}_filtered.vcf.gz ${vcf_file}
#tabix -p vcf ${sample}_filtered.vcf.gz