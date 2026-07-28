#!/bin/bash

#This script not only trims SNPs by its p-value but also
#separates each SNP so that they can be processed in future steps.

#Set the input as the first argument passed to the script
INPUT=$1

#Define the max p-value
MAXP=$2

#Define PALM analysis
ANALYSIS=$3

#Create a temporary file
TMP=TMP.txt

#Remove the headers
tail -n +2 "$INPUT" | \

#Trim SNPs by maximum p-value. This is incredibly helpful to
#avoid unnecessary sampling branch lengths and inferring SNPs
#likelihoods that will not be used later on. Keep those with
#p-val <= 5e-8 for PALM and those with p-val <= 5e-8 in the first 
#GWAS OR the second GWAS for JPALM.
if [ "$ANALYSIS" = "JPALM" ]; then
    awk -F'\t' -v maxp="$MAXP" '($10 <= maxp) || ($15 <= maxp)' "$INPUT_FILE" > "$TMP"
elif [ "$ANALYSIS" = "PALM" ]; then
    awk -F'\t' -v maxp="$MAXP" '($10 <= maxp)' "$INPUT_FILE" > "$TMP"
else
    echo "Analysis not recognized. Please use either 'PALM' or 'JPALM'."
    exit 1
fi

#Initialize batch size and counter
batch_size=$4
batch_counter=0
batch_number=1

#Read the TSV file and split into separate files
while IFS= read -r line || [[ -n "$line" ]]; do
    #Increment the counter
    ((batch_counter++))

    #Generate the output file name for the current batch
    output_file="batch_${batch_number}_${INPUT}"

    #Append the line to the current batch file
    echo -e "$line" >> "$output_file"

    #Check if the batch is full, and reset counter
    if ((batch_counter % batch_size == 0)); then
        ((batch_number++))
    fi
done < "$TMP"

#Remove TMP file
rm "$TMP"