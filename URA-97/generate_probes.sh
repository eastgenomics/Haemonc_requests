# TWIST have sent two folders 1X and 2X of probes.
# We need to merge these two probe files to get one for the whole
# capture. This probe file is used in Picard which is required for
# the multiQC report.

# cut first three columns as the fourth col contains genes for some probes

cat Probes_merged_ok_Cambridge_NHS_Myeloid_1X_TE-96527893_hg38_211119151014.bed Probes_merged_ok_Cambridge_NHS_Myeloid_2X_TE-96527893_hg38_211119151112.bed | cut -f 1,2,3 > Probes_GRCh38_HaemOnc_v2.0.bed 
