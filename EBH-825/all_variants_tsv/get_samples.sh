
projects=$(dx find projects --name "002_*MYE")

for project in $projects;  do      echo $project; file_id=$(dx find data --json --path="$project" --name="*_allgenes.tsv" | jq -r '.[].id');  for file in $file_id; do         dx download $file --no-progress; done; done


# ! The last column is wrapped so for everyline there is extra 5 rows, so to know wc -l 
# we need to divide by 6 

# count number of lines in each file and add to count
count=0
for file in *.tsv;
    do echo $file;
    let count=count+$(cat $file | wc -l | awk '{ print ($1-1)/6 }');
    echo $count;
done

# concat all samples skipping the first line for every flile except the first
awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > all_variants_220124.tsv

# check that the number of lines in all == count
echo $count
wc -l all_variants_220124.tsv | awk '{ print ($1-1)/6 }'
