## Run this script when current location is in the data file
mkdir anonymised_outputs

# As we need the sample names to be random, we will assign a random five
# digit number
# All output samples will be elsewhere in another folder called anonymised_outputs
for j in $(find . -type f -name "*.tsv");
    do echo $j;
    random_num=${RANDOM:0:5};
    echo S$random_num;
    cut -f 2- $j > anonymised_outputs/S${random_num}_markdup_recalibrated_tnhaplotyper2_allgenes.tsv
done > anonymised_outputs/log_out.txt
