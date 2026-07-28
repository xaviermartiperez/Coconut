#!/bin/bash

## This script intends to add the reference allele, 
## chromosome and position column to GWAS summary statistics
## when needed.

# Define objects
GWAS_INFO=$1
GWAS_DIR=$2
VCF_DIR=$3

read gwas_name genome_version missing_ref_alt missing_chr_pos < "$GWAS_INFO"

GWAS_FILE_PATH=${GWAS_DIR}/${gwas_name}
FINAL_GWAS_FILE_PATH="${GWAS_DIR}/$(echo "$gwas_name" | sed 's/formatted_/fixed_/')"

if [[ "$missing_ref_alt" == FALSE && "$missing_chr_pos" == FALSE ]]; then
    # Emit the expected file name
    mv "$GWAS_FILE_PATH" "$FINAL_GWAS_FILE_PATH"
    # Emit the expected output
    echo -e "$(basename "$FINAL_GWAS_FILE_PATH")\t$genome_version" > "$(basename ${FINAL_GWAS_FILE_PATH%.tsv}).txt" 
    exit 0
else
    #Load necessary packages
    module load BCFtools/1.21-GCC-12.3.0

    # Define according dbSNP file based on the genome version
    #if [[ "$genome_version" == "hg38" ]]; then
        # Whatever
    #elif [[ "$genome_version" == "hg19" ]]; then
        # Whatever
    #else
    #    echo "Invalid genome version. Please use 'hg38' or 'hg19'."
    #    exit 1
    #fi

    #Define an intermediate file to store the raw rsIDs
    INTER=list_rsIDs.txt.tmp
        
    #Extract unique rsIDs from the GWAS
    awk 'NR>1 {print $11}' "$GWAS_FILE_PATH" | sort -u > $INTER

    #Build the extraction format for bcftools based on what is missing
    FORMAT="%ID" 
    if [[ "$missing_chr_pos" == "TRUE" ]]; then
        FORMAT="$FORMAT\t%CHROM\t%POS"
    fi
    if [[ "$missing_ref_alt" == "TRUE" ]]; then
        FORMAT="$FORMAT\t%REF"
    fi

    for CHR in {1..22}; do
        #Name the appropriate dbSNP dataset for this chromosome
        VCF_FILE="${VCF_DIR}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

        OUT="vcf_data_chr${CHR}.tsv.tmp"

        #Query dbSNP for the current chromosome
        bcftools query -f "$FORMAT\n" -i "ID=@$INTER" "$VCF_FILE" > "$OUT"

        
        #Add header to the output file
        HEADER="rsID"
        if [[ "$missing_chr_pos" == "TRUE" ]]; then
            HEADER="$HEADER\tchr\tpos"
        fi
        if [[ "$missing_ref_alt" == "TRUE" ]]; then
            HEADER="$HEADER\tref"
        fi

        sed -i "1i $HEADER" "$OUT"

    done

    #Merge all the datasets into one
    FINAL_TEMP=temp_vcf_data.tsv
    HEADER_WRITTEN=false
    > "$FINAL_TEMP"  #Create an empty file

    #Iterate over the file names and concatenate them
    for CHR in {1..22}; do
        CHR_FILE="temp_vcf_data_chr${CHR}.tsv"
        #Check if it's the first file and write the header
        if [[ "$HEADER_WRITTEN" = false ]]; then
            cat "$CHR_FILE" >> "$FINAL_TEMP"
            HEADER_WRITTEN=true
        else
            #Skip the header and append the contents
            tail -n +2 "$CHR_FILE" >> "$FINAL_TEMP"
        fi
    done
  
    TEMP_GWAS="temp_gwas.tsv"
    
    #Subtitute in the GWAS the missing information 
    awk -v missing_chr_pos="$missing_chr_pos" -v missing_ref_alt="$missing_ref_alt" '
    NR==FNR {
        if (missing_chr_pos == "TRUE" && missing_ref_alt == "TRUE") {
            data[$1] = $2 "\t" $3 "\t" $4;
            variant[$1] = $2 ":" $3;
        } else if (missing_chr_pos == "TRUE") {
            data[$1] = $2 "\t" $3;
            variant[$1] = $2 ":" $3;
        } else if (missing_ref_alt == "TRUE") {
            data[$1] = $2;
        }
        next;
    }
    FNR==1 {print; next}
    $11 in data {
        if (missing_chr_pos == "TRUE") {
            $1 = variant[$11];
            split(data[$11], vals, "\t");
            $9 = vals[1];
            $10 = vals[2];
        }
        if (missing_ref_alt == "TRUE") {
            split(data[$11], vals, "\t");
            $7 = vals[length(vals)];

            if ($7 == $12) {
                # If reference allele is the effect allele, change the beta to assign the effect to the alt
                $3 = -$3;
                # Alt allele is the non effect allele
                $8 = $13
            } else {
                # If the reference allele is not the effect allele, assign the effect allele to alt
                $8 = $12}
        }
        print;
    }' OFS='\t' "$FINAL_TEMP" "$GWAS_FILE_PATH" > "$TEMP_GWAS.tmp"

    rm "$GWAS_FILE_PATH"

    #Remove the effect_allele column if missing_ref_alt = TRUE
    if [[ "$missing_ref_alt" == "TRUE" ]]; then
        awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' "$TEMP_GWAS.tmp" > "$TEMP_GWAS"
        rm "$TEMP_GWAS.tmp"  #Remove the temporary file
    else
        mv "$TEMP_GWAS.tmp" "$TEMP_GWAS"
    fi
    
    #Move the final file
    mv "$TEMP_GWAS" "$FINAL_GWAS_FILE_PATH"
    echo -e "$(basename "$FINAL_GWAS_FILE_PATH")\t$genome_version" > "$(basename ${FINAL_GWAS_FILE_PATH%.tsv}).txt" 

    #Remove the useless files
    rm -f "$INTER" "$TEMP_GWAS" temp_vcf_data_chr*.tsv
fi