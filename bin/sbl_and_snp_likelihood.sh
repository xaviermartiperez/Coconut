#!/bin/bash

#Load Relate and Phyton modules
module load Relate/1-GCCcore-10.2.0
module load Python/3.9.6-GCCcore-11.2.0

#Define the variables that don't change
MU=1.25e-8
BATCH=$1
POP=$2
RELATE_TOOL=$3
LIK_TOOL=$4
COAL_FILE=$5
ANC_MUT_DIR=$6
LIK_OUTDIR=$7

#Iterate over each line (row) of the provided batch of SNPs
while IFS=$'\t' read -r BP foo_1 foo_2 CHR VARIANT LD DAF rest_of_the_columns
do
    POP_CHR_ANC_MUT=${POP}_chr${CHR}_upd
    INPUT_RELATE=${ANC_MUT_DIR}/${POP_CHR_ANC_MUT}
    OUTPUT_RELATE=${VARIANT}_out

    #Call the script and store it to obtain the SNP likelihood (max 5min to do so)
    timeout 5m $RELATE_TOOL -i $INPUT_RELATE -o $OUTPUT_RELATE \
    -m $MU --coal $COAL_FILE --format b --num_samples 5 --first_bp $BP --last_bp $BP --seed 1 || true

    #Remove the files that are not necessary
    rm *_out.anc
    rm *_out.mut
    rm *_out.dist

    #Define the input for the SNP Likelihood script
    INPUT_LIK=${OUTPUT_RELATE}.timeb

    #Check if the input of the SNP likelihood script exists and is not empty.
    #If there is no file, skip the rest of the iteration and start the next one    
    if [ ! -s "$INPUT_LIK" ]; then
        # If the file does not exist or is empty, move to the next iteration
        echo "Skipping iteration, $INPUT_LIK does not exist or is empty."
        continue
    fi

    #Create the necessary directory
    OUTDIR=${LIK_OUTDIR}/ld_${LD}
    mkdir -p "$OUTDIR"

    #Define the ouput file name "root"
    OUTPUT_LIK=bp${BP}

    #Call the script and store it to obtain the SNP likelihood (max 5min to do so)
    timeout 5m python $LIK_TOOL --times $OUTPUT_RELATE \
    --popFreq $DAF --out $OUTPUT_LIK --coal $COAL_FILE || true

    #Remove the files that are not necessary (again)
    rm *_out.timeb

    #Copy the output to the specified folder if the file exists and is not empty
    if [ -e "${OUTPUT_LIK}.quad_fit.npy" ] && [ -s "${OUTPUT_LIK}.quad_fit.npy" ]; then
        cp ${OUTPUT_LIK}.quad_fit.npy $OUTDIR/.
    fi

done < $BATCH
