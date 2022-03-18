# tmp contains the list of files in current directory
tmp=( $(ls | sort | uniq))
# number of files to loop through, take one away as loop starts at 0
len=$(($(ls | wc -l)-1))

mkdir anonymised_outputs

# As we need the sample names to be random, we will assign a random five
# digit number
# All output samples will be else in another folder called anonymised_outputs
for j in `seq 0 1 $len`;
    do echo ${tmp[$j]};
    random_num=${RANDOM:0:5};
    echo S$random_num;
    python ../bin/anonymise_samples.py -f ${tmp[$j]}  -o S${random_num}_markdup_recalibrated_tnhaplotyper2_panels -p anonymised_outputs
done > anonymised_outputs/log_out.txt