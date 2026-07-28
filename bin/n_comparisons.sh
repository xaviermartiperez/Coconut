#!/bin/bash

#This script takes the UKBB list with the selected GWAS and stores the URL of each
#of them in a separate file. This allows Nextflow to parallelise the following 
#processes.

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

#Assign the file to a variable
SELECTED_GWAS=selected_GWAS.txt

#Get the number of pairwise comparisons (for the Bonferroni testing of the genetic
#correlation)
#Count the number of rows in the data frame (excluding the header)
NUM_ROWS=$(wc -l < $SELECTED_GWAS)

# Calculate the number of pairwise comparisons
N_COMPARISONS=$((NUM_ROWS * (NUM_ROWS - 1) / 2))

#Put the result in a file
echo $N_COMPARISONS > n_pairwise_comparisons.txt

#Remove unnecesary file
rm selected_GWAS.txt