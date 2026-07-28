#!/bin/bash

#Load Phyton module
module load Python/3.9.6-GCCcore-11.2.0

#Define inputs
METADATA=$1
MAXPVAL=$2
PALM_TOOL=$3
PATH_TO_LIK=$4

#Create directory and if it already exists, the existing directory is used.
mkdir -p $OUTDIR

#Obtain file with minimum number of independent, significant SNPs in each GWAS
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i ~ /^min_independent_significant_SNPs@/) header[i]=$i} NR==2 {for (i=1; i<=NF; i++) if (header[i]) print header[i], $i}' \
$1 > $(basename "$METADATA" _selected_SNPs.tsv)_significant_independent_SNPs.txt

#Separate the GWAS in the name of the file
INFO=($(echo "$1" | tr '-' ' '))

#Modify the elements so that they match the names in the file
GWAS_1="${INFO[0]}.tsv"
GWAS_2="${INFO[1]%_selected_SNPs.tsv}.tsv"

#Call the script
python $PALM_TOOL \
    --traitDir ${PATH_TO_LIK}/ \
    --metadata $METADATA \
    --maxp $MAXPVAL \
    --B 10000 \
    --traits $GWAS_1,$GWAS_2 \
    > $(basename "$METADATA" _selected_SNPs.tsv)_J_PALM.txt