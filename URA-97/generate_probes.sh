# cut first three columns as the fourth contains genes for some probes

cat Probes_merged_ok_Cambridge_NHS_Myeloid_1X_TE-96527893_hg38_211119151014.bed Probes_merged_ok_Cambridge_NHS_Myeloid_2X_TE-96527893_hg38_211119151112.bed | cut -f 1,2,3 > Probes_GRCh38_HaemOnc_v2.0.bed 
