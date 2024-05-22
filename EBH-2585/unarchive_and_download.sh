#!/bin/bash

output_dir="all_filtered_vcfs"
mkdir -p "$output_dir"

# get ids of all prod MYE projects since 31/05/2022
projects=$(dx find projects --name 002_*_MYE --created-after 2022-05-31 --brief)

for project in $projects; do
    project_name=$(dx describe "$project" --json | jq -r '.name')

    # get annotated VCFs for all cases in that project
	files=$(dx find data --name *_markdup_recalibrated_tnhaplotyper2_allgenesvep.vcf --project "$project" --brief)

    # initialise list of files to unarchive
    to_unarchive=""

	for file in $files; do

        # identify file's name and archival state
        name=$(dx describe "$file" --json | jq -r '.name')
        state=$(dx describe "$file" --json | jq -r '.archivalState')

        # get sample id from filename, construct name for downloaded file
        sample="${name%%_markdup_recalibrated_tnhaplotyper2_allgenesvep.vcf}"
        dl_name="${output_dir}/${project_name}_${sample}.vcf"

        # if file isn't already downloaded,
        if [[ ! -e "$dl_name" ]]; then

            # if live, download it
            if [[ "$state" == "live" ]]; then
                echo "${project_name} ${sample} - downloading"
                dx download "$file" -o "$dl_name" -f --no-progress

            # if archived, add to list of files to unarchive
            elif [[ "archived archival" == *"$state"* ]]; then
                to_unarchive="${to_unarchive}${file} "

            # otherwise just print what it's current state is
            else
                echo "${project_name} ${sample} - ${state}"
            fi
        # else
        #     echo "${project_name} ${sample} - already downloaded"
        fi
    done

    # if to_unarchive isn't an empty string, unarchive the files in it
    if [[ "${to_unarchive// /}" != "" ]]; then
        echo "${project_name} - starting unarchiving"
        dx unarchive "$to_unarchive" -y
    else
        echo "${project_name} - no files to unarchive"
    fi
done
