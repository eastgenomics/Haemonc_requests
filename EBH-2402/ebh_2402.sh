#!/bin/bash

project="$1"

# output of Sentieon somatic tnhaplotyper 2
small_var_vcfs=$(dx find data --name *_markdup_recalibrated_tnhaplotyper2.vcf.gz --project "$project" --brief)

# output of cgppindel
indel_vcfs=$(dx find data --name *_vs_TA2_S59_L008_tumor.flagged.vcf.gz --project "$project" --brief)

# MYE capture bed
mye_bed="coding_unrestricted_GRCh38_myeloid_5bp_flank_v2.0.0.bed"

# create folders
printf "\nUsing DNAnexus project: %s\n\nCreating folders\n" "$project"

folders="vcfs vcfs/${project} output output/${project}"

for folder in $folders; do
    mkdir "$folder"
done

# download and process small variant VCFs
printf "\nProcessing small variant VCFs\n\n"

for vcf in $small_var_vcfs; do

    name=$(dx describe "$vcf" --json | jq -r '.name')
    dl_path="vcfs/${project}/${name}"
    unzipped="${dl_path%%.gz}"

    # if file not already downloaded, check if archived
    if [[ ! -e "$dl_path" ]] && [[ ! -e "$unzipped" ]]; then

        state=$(dx describe "$vcf" --json | jq -r '.archivalState')

        # if live, download and unzip it
        if [[ "$state" == "live" ]]; then
            dx download "$vcf" -o "$dl_path" -f --no-progress
            gunzip "$dl_path"

        # if not live, print message
        else
            echo "${name} could not be downloaded, file is ${state}"
        fi
    fi

    # if file previously downloaded but not unzipped, unzip it
    if [[ -e "$dl_path" ]]; then
        gunzip "$dl_path"
    fi

    # if unzipped file is present, intersect it with the MYE bed
    # bedtools -u: each entry in A is output once if it occurs anywhere in B
    if [[ -e "$unzipped" ]]; then
        output="output/${project}/${name%%.vcf.gz}_intersect_small_vars.vcf"
        bedtools intersect -a "$unzipped" -b "$mye_bed" -u > "$output"
    else
        echo "${name} not processed, file not available"
    fi
done

# download and process indel VCFs
printf "\nProcessing indel VCFs\n\n"

for vcf in $indel_vcfs; do

    name=$(dx describe "$vcf" --json | jq -r '.name')
    dl_path="vcfs/${project}/${name}"
    unzipped="${dl_path%%.gz}"

    # if file not already downloaded, check if archived
    if [[ ! -e "$dl_path" ]] && [[ ! -e "$unzipped" ]]; then

        state=$(dx describe "$vcf" --json | jq -r '.archivalState')

        # if live, download and unzip it
        if [[ "$state" == "live" ]]; then
            dx download "$vcf" -o "$dl_path" -f --no-progress
            gunzip "$dl_path"

        # if not live, print message
        else
            echo "${name} could not be downloaded, file is ${state}"
        fi
    fi

    # if file previously downloaded but not unzipped, unzip it
    if [[ -e "$dl_path" ]]; then
        gunzip "$dl_path"
    fi

    # if unzipped file is present, process it
    if [[ -e "$unzipped" ]]; then
        output="output/${project}/${name%%.vcf.gz}_intersect_indels.vcf"
        bedtools intersect -a "$unzipped" -b "$mye_bed" -u > "$output"
    else
        echo "${name} not processed, file not available"
    fi
done

# if output files available, get average variant numbers

output_path="output/${project}/"

if [ "$(ls -A "$output_path")" ]; then
    printf "\nGetting statistics\n\n"
    python ebh_2402.py "$output_path"

else
    printf "\nNo intersected VCFs generated, cannot generate statistics\n\n"
fi
