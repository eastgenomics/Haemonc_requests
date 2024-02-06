#!/bin/bash

project="$1"

# output of Sentieon somatic tnhaplotyper 2
small_var_vcfs=$(dx find data --name *_markdup_recalibrated_tnhaplotyper2.vcf.gz --project "$project" --brief)

# output of cgppindel
indel_vcfs=$(dx find data --name *_vs_TA2_S59_L008_tumor.flagged.vcf.gz --project "$project" --brief)

# MYE capture bed (project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9vKbv0433Gg2bP9GP72pxKJ)
mye_bed="coding_unrestricted_GRCh38_myeloid_5bp_flank_v2.0.0.bed"

# create folders

echo "Using DNAnexus project: ${project}"
echo " "
echo "Creating folders"

folders="vcfs vcfs/${project} output output/${project}"

for folder in $folders; do
    mkdir "$folder"
done

# download small variant and indel VCFs

echo "Downloading files"

for vcf in $small_var_vcfs; do

    name=$(dx describe "$vcf" --json | jq -r '.name')
    dl_path="vcfs/${project}/${name}"
    unzipped="${dl_path%%.gz}"

    if [[ ! -e "$dl_path" ]] && [[ ! -e "$unzipped" ]]; then
        echo "Downloading ${name}"
        dx download "$vcf" -o "$dl_path" -f --no-progress
    fi

    if [[ -e "$dl_path" ]]; then
        gunzip "$dl_path"
    fi
done

for vcf in $indel_vcfs; do

    name=$(dx describe "$vcf" --json | jq -r '.name')
    dl_path="vcfs/${project}/${name}"
    unzipped="${dl_path%%.gz}"

    if [[ ! -e "$dl_path" ]] && [[ ! -e "$unzipped" ]]; then
        echo "Downloading ${name}"
        dx download "$vcf" -o "$dl_path" -f --no-progress
    fi

    if [[ -e "$dl_path" ]]; then
        gunzip "$dl_path"
    fi
done

# intersect VCFs with bed file
# -u: write original A entry once if any overlaps are found in B

echo "Intersecting files"

for vcf in "vcfs/${project}/"*_markdup_recalibrated_tnhaplotyper2.vcf; do

    name=$(basename "$vcf")
    output="output/${project}/${name%%.vcf}_intersect_small_vars.vcf"

    if [[ ! -e "$output" ]]; then
        bedtools intersect -a "${vcf}" -b "$mye_bed" -u > "$output"
    fi
done

for vcf in "vcfs/${project}/"*_vs_TA2_S59_L008_tumor.flagged.vcf; do

    name=$(basename "$vcf")
    output="output/${project}/${name%%.vcf}_intersect_indels.vcf"

    if [[ ! -e "$output" ]]; then
        bedtools intersect -a "${vcf}" -b "$mye_bed" -u > "$output"
    fi
done

# get variant averages

echo "Getting statistics"
echo " "
output_path="output/${project}/"
python ebh_2402.py "$output_path"
