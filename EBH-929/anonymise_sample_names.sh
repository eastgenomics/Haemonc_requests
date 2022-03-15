# 14/03/2022 AD

# Anonymise sample names
# As the order of the fastq files do not follow the samplesheet order,
# we can assign s1, s2 as we loop through the list of sample fastqs


mkdir annonymised_outputs
tmp=( $(ls | sort | uniq))
len=( $(ls | wc -l) / 4)
for j in `seq 0 4 $len`; 
    do echo $j; 
    for ((i=$j; i<$(($j + 4)); i++)); 
        do echo ${tmp[$i]}; scp ${tmp[$i]}  annonymised_outputs/$(echo ${tmp[$i]} | sed "s/.*MYE/S${j}/g");
    done;
done


