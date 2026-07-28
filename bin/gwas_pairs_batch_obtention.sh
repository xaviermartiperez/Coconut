#!/bin/bash

#This script intends to create a new file that contains each possible pairwise
#comparison in the GWAS list provided

#Exit immediately if any command within it returns a non-zero exit (error) status
set -e
#Treat references to unset variables as errors and make the script exit if so
set -u

#Create a merged file with the GWAS that have >=25 independent SNPs from all arguments
#passed (the result of .collect())
touch selected_GWAS.txt

for file in "$@"; do
  # Extract the number from the filename using parameter expansion
  number=${file#*-+}
  # Check if the number is greater than or equal to 25
  if [ "$number" -ge 25 ]; then
    cat "$file" >> selected_GWAS.txt
  fi
done

SELECTED_GWAS=selected_GWAS.txt

filenames=($(cat "$SELECTED_GWAS"))

#Create list with all pairs
touch all_comparisons.txt

#Loop through all pairwise combinations
for ((i=0; i<${#filenames[@]}; i++)); do
    for ((j=i+1; j<${#filenames[@]}; j++)); do
        #Extract the corresponding rows from the data frame and save them to the new file
        echo -e "${filenames[i]}\t${filenames[j]}" >> all_comparisons.txt
    done
done

#Split the file with all comparisons 
split --lines=10 --numeric-suffixes=1 --suffix-length=3 all_comparisons.txt comparison_batch_

#Clean up indesired files
rm selected_GWAS.txt
rm all_comparisons.txt