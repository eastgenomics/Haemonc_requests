#!/bin/bash

cd "/home/arun/Codes/HaemOnc_panel_2024-Twist_capture/Validation/plots_to_present/vcf_files" || exit

for file in SP/new_probes/*.vcf; do bgzip "$file"; done
for file in SP/old_probes/*.vcf; do bgzip "$file"; done
for file in S1/new_probes/*.vcf; do bgzip "$file"; done
for file in S1/old_probes/*.vcf; do bgzip "$file"; done

for file in SP/new_probes/*.vcf.gz; do tabix -p vcf "$file"; done
for file in SP/old_probes/*.vcf.gz; do tabix -p vcf "$file"; done
for file in S1/new_probes/*.vcf.gz; do tabix -p vcf "$file"; done
for file in S1/old_probes/*.vcf.gz; do tabix -p vcf "$file"; done

filter_PASS_and_index_vcfs() {

    # Set up the function inputs
    local input_dir="$1"
    local output_dir="$2"

    # Step 1: Remove and recreate the output directory
    rm -rf "$output_dir" && mkdir -p "$output_dir"

    # Step 2: Filter VCF files
    for vcf_file in "$input_dir"/*.vcf.gz; do
        local sample=$(basename "$vcf_file" .vcf.gz)
        bcftools filter -i 'FORMAT/DP > 250 && FORMAT/AF > 0.05 && FILTER="PASS"' \
            -Oz -o "$output_dir/${sample}_filtered.vcf.gz" "$vcf_file"
    done

    # Step 3: Index the filtered VCF files
    for file in "$output_dir"/*.vcf.gz; do
        tabix -p vcf "$file"
    done

}


filter_withoutPASS_and_index_vcfs() {

    # Set up the function inputs
    local input_dir="$1"
    local output_dir="$2"

    # Step 1: Remove and recreate the output directory
    rm -rf "$output_dir" && mkdir -p "$output_dir"

    # Step 2: Filter VCF files
    for vcf_file in "$input_dir"/*.vcf.gz; do
        local sample=$(basename "$vcf_file" .vcf.gz)
        bcftools filter -i 'FORMAT/DP > 250 && FORMAT/AF > 0.05' \
            -Oz -o "$output_dir/${sample}_filtered.vcf.gz" "$vcf_file"
    done

    # Step 3: Index the filtered VCF files
    for file in "$output_dir"/*.vcf.gz; do
        tabix -p vcf "$file"
    done

}

# Example
filter_PASS_and_index_vcfs "SP/new_probes" "SP/new_probes_filtered_PASS"
filter_PASS_and_index_vcfs "SP/old_probes" "SP/old_probes_filtered_PASS"
filter_PASS_and_index_vcfs "S1/new_probes" "S1/new_probes_filtered_PASS"
filter_PASS_and_index_vcfs "S1/old_probes" "S1/old_probes_filtered_PASS"
filter_withoutPASS_and_index_vcfs "SP/new_probes" "SP/new_probes_filtered"
filter_withoutPASS_and_index_vcfs "SP/old_probes" "SP/old_probes_filtered"
filter_withoutPASS_and_index_vcfs "S1/new_probes" "S1/new_probes_filtered"
filter_withoutPASS_and_index_vcfs "S1/old_probes" "S1/old_probes_filtered"