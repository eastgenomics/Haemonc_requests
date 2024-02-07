# EBH-2402

Original request: "Ilenia would like to know how many variants are called per case (in general) within the target regions of the MYE assay, prior to any filtering."

The ebh_2402.sh Bash script takes a DNAnexus project ID as input. It downloads the unfiltered VCFs, intersects them with the MYE assay bed file, counts the numbers of variants in the output of intersection, and generates descriptive statistics via the ebh_2402.py Python script.

## Resources used

The small variant VCFs used are the outputs of Sentieon somatic tnhaplotyper 2 (*_markdup_recalibrated_tnhaplotyper2.vcf.gz).

The MYE bed file used is 'coding_unrestricted_GRCh38_myeloid_5bp_flank_v2.0.0.bed' (project-Fkb6Gkj433GVVvj73J7x8KbV:file-G9vKbv0433Gg2bP9GP72pxKJ).

The five most recent 002_*_MYE runs (as of 06/02/2024) are:

- 002_240202_A01295_0311_AHWJWFDRX3_MYE (project-Gfz6X9j4b57gJ1Kzqg4Qp88y)
- 002_240130_A01303_0330_AHWL32DRX3_MYE (project-Gfx83v04GkZXQXK2Q63vYj4Q)
- 002_240117_A01303_0325_AHW5H5DRX3_MYE (project-GfZZP5Q4Bp89jZ80yyGvQKFY)
- 002_240105_A01295_0292_BHTH75DRX3_MYE (project-GfGjKF04Q5Z5Z47z31QzJgXV)
- 002_240105_A01295_0291_AHVJF5DRX3_MYE (project-GfGjGZj4XZq86PPy80GFqKXK)

## Outputs

The outputs of ebh_2402.sh for each of the above projects, as well as the results across all five projects, are included in EBH-2402_results.xlsx.
