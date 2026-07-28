#!/bin/bash

# Load Phyton module
module load Python/3.9.6-GCCcore-11.2.0

#Define inputs
GWAS=$1
MAXPVAL=$2
PALM_TOOL=$3
LIK_DIR=$4


#Call the script
python $PALM_TOOL \
    --traitDir ${LIK_DIR}/ \
    --metadata $GWAS \
    --maxp $MAXPVAL \
    --B 10000 \
    > $(basename "$GWAS" _selected_SNPs.tsv)_marginal_PALM.txt

