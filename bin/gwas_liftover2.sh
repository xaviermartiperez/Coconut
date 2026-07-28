#!/bin/bash

#This script creates a BED-formatted file for each selected SNPs
#GWAS file. Then executes the liftover tool. Finally it subsitutes
#the position in the GWAS with the new coordinates.

#Define GWAS
GWAS=$(cat ${1})

#Define input directory
INDIR="$2"

#Define GWAS list
LIST="$3"

#Read the 4th column of the list and check if it's 'true' or 'false'
lift_over=$(awk 'NR==1 { print $4 }' "$LIST")

#Define genome
GENOME="$4"

#Define output directory for GWAS lifted
OUTDIR="/homes/users/ajimenez/scratch/gwas_tfm/lifted"

# Create directory and if it already exists, the existing directory is used.
mkdir -p $OUTDIR

#Define the final output file for GWAS lifted
GWAS_OUTPUT="${OUTDIR}/${GWAS%.tsv}_${GENOME}.tsv"

if [ "$lift_over" == "true" ]; then
    if [ -e "$GWAS_OUTPUT" ]; then
        #If the lifted file exists, emit the expected output file and exit 
        echo "${GWAS_OUTPUT}" > ${GWAS%.tsv}_${GENOME}.txt
        exit 0
    fi
    
    #Define BED file obtained from the GWAS
    TEMP_BED=${GWAS%.tsv}.bed

    #Define LiftOver output file
    OUTPUT_LIFTOVER=${GWAS%.tsv}_lifted.bed

    #Define file with unlifted SNPs
    UNLIFTED=${OUTDIR}/${GWAS%.tsv}_${GENOME}_unlifted.txt

    #Specify the path to the LiftOver tool
    LIFTOVER_TOOL="/homes/users/ajimenez/scratch/gcat/liftOver"

    #Specify the chain file for the genome conversion (needed for the liftOver tool to work)
    if [[ "$GENOME" == "hg38" ]]; then
        CHAIN_FILE="/homes/users/ajimenez/scratch/gcat/hg19ToHg38.over.chain"
    else
        CHAIN_FILE="/homes/users/ajimenez/scratch/gcat/hg38ToHg19.over.chain"
    fi

    #Create pseudo-BED file for the LiftOver tool to work
    awk 'BEGIN { OFS="\t" } 
    NR>1 {
        chr = "chr" $9;
        start = $10 - 1;
        end = $10;
        variant = $1;
        print chr, start, end, variant
    }' ${INDIR}/${GWAS} > $TEMP_BED

    #Call the LiftOver tool
    "$LIFTOVER_TOOL" "$TEMP_BED" "$CHAIN_FILE" "$OUTPUT_LIFTOVER" "$UNLIFTED"

    #Add the new variant column as a new column to the OUPUT_LIFTOVER file.
    #Put the output in the TEMP_BED file (thus recycling it)
    awk 'BEGIN { OFS="\t" } { gsub("chr", "", $1); new_col = $1 ":" $3; print "chr"$1, $2, $3, $4, new_col }' $OUTPUT_LIFTOVER > $TEMP_BED
    
    #Add headers (just in case)
    sed -i '1ichr\tpos_start\tpos_end\told_variant\tnew_variant' $TEMP_BED

    #Change the old variant with the new variant in the GWAS file
    #Create a temporary file to store the mappings from old_variant to new_variant
    temp_mapping_file=$(mktemp)

    #Read the second file and populate the temporary mapping file
    while IFS=$'\t' read -r chr pos_start pos_end old_variant new_variant; do
        # Skip the header line
        if [ "$old_variant" != "old_variant" ]; then
            echo -e "${old_variant}\t${new_variant}" >> "$temp_mapping_file"
        fi
    done < $TEMP_BED

    #Process the first file and replace old_variant with new_variant where applicable, removing unmatched rows
    awk -v OFS="\t" '
    BEGIN {
        # Read the mapping file into an associative array
        while ((getline < "'$temp_mapping_file'") > 0) {
            variant_map[$1] = $2
        }
    }
    NR==1 {
        # Print the header
        print $0
    }
    NR>1 {
        # Only print the row if a mapping exists
        if ($1 in variant_map) {
            variant = variant_map[$1]
            $1 = variant
            split(variant, pos, ":");
            $9 = pos[1];
            $10 = pos[2];
            print $0
        }
    }
    ' ${INDIR}/${GWAS} > $GWAS_OUTPUT

    #Remove the temporary file
    rm "$temp_mapping_file" "$TEMP_BED" "$OUTPUT_LIFTOVER"

    #Emit the expected output file
    echo "${GWAS_OUTPUT}" > ${GWAS%.tsv}_${GENOME}.txt
    exit 0

else
    #If the lift Over is not necessary, emit the exepected output
    ######### REVISAR INDIR FUNCIONARA?
    echo "${INDIR}/${GWAS}" > ${GWAS%.tsv}_${GENOME}.txt
    exit 0
fi