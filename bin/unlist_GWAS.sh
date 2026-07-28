#!/bin/bash

# This script takes the GWAS list and stores the path of each GWAS 
# in a separate file. This allows Nextflow to parallelise the following 
# processes.

# Exit immediately if any command within it returns a non-zero exit (error) status
set -e
# Treat references to unset variables as errors and make the script exit if so
set -u

# Define objects
GWAS_LIST=$1

# Create separate text files for each gwas (skips header). Also sets NA for missing
# columns 3 and/or 4.
tail -n +2 "$GWAS_LIST" | awk '{
    # Check if column 3 exists; if not, set to NA
    if ($3 == "") $3 = "NA"
    # Check if column 4 exists; if not, set to NA
    if ($4 == "") $4 = "NA"
    filename = "gwas_" NR ".txt"
    print $1 "\t" $2 "\t" $3 "\t" $4 > filename
}'
