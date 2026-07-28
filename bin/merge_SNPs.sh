#!/bin/bash

#This script merged the selected SNPs, then sorts them and
#subsequently removes any duplicate entries

#Define input
INPUT=("$@")

#Sort the input by the variant column, remove duplicates and stores it in a file
sort -k5,5 -u "${INPUT[@]}" > merged_SNPs.tsv

#Capture the last line, which refers to the header
headers=$(tail -n 1 merged_SNPs.tsv)

#Remove the last line of the file by printing in a temporary file the previous line for each iteration
awk 'NR>1{print last} {last=$0}' merged_SNPs.tsv > temp_file.tsv

#Rebuild the file with the header on top
echo -e "$headers\n$(cat temp_file.tsv)" > merged_SNPs.tsv

#Remove the temporary file
rm temp_file.tsv