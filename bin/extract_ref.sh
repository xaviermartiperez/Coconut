#!/bin/bash
#SBATCH --job-name=ref_alleles
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH -o slurm-%j.out 
#SBATCH --partition=haswell
#SBATCH --mail-user=ainara.jimenez@upf.edu
#SBATCH --mail-type=ALL

##This script intends to add the variant information (chr:pos:ref:alt)
##to the raw reQTL data that Rotival and Quintana-Murci sent to me.
##Coordinates are for Grch37/hg19

#Load necessary packages
module load BCFtools 

##Define objects
#Define the file to which the variant informations need to be added (in .tsv format)
IN=/homes/users/ajimenez/scratch/gwas_tfm/GCST90092809_buildGRCh37.tsv.gz
#Define the column numbers where the IDs and chromosomes are stored
ID_COLUMN=1
#Define an intermediate file to store the raw rsIDs
INTER=/homes/users/ajimenez/scratch/gwas_tfm/unique_rsIDs.txt
#Define the output file (without the extension)
BASE_OUT=/homes/users/ajimenez/scratch/gwas_tfm/id_ref_GCST90092809
#Define the SNPs databases path that are inteded to be used for the 
#annotation based on the rsID
dbSNP=/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf/
#Define temp file
TEMP=/homes/users/ajimenez/gwas_tfm/gcat/temp.tsv

#Create an array of chromosomes (1 to 22)
CHROMOSOMES=( $(seq 1 22) )

##Loop to do the whole process per chromosome
for CHR in "${CHROMOSOMES[@]}"; do
  # Name the output file for this chromosome
  OUT=${BASE_OUT}_chr${CHR}.tsv
  
  # Name the appropriate dbSNP dataset for this chromosome
  dbSNP_CHR=${dbSNP}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

  # Extract unique rsIDs from the GWAS
  gunzip -c "$IN" | cut -f $ID_COLUMN | sort -u > $INTER

  # Query dbSNP for the current chromosome
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\n' -i "ID=@$INTER" $dbSNP_CHR > $OUT
  
  # Add headers to the output file
  sed -i '1i chr\tpos\trsID\tref' $OUT
  
  # Create chr:pos column
  awk -F '\t' 'NR==1 {print $0, "chr:pos"} NR>1 {print $0, $1":"$2}' OFS='\t' $OUT > ${OUT}.tmp
  mv ${OUT}.tmp $OUT
done

##Merge all the datasets into one
#Name the final output file
OUTPUT_FILE=${BASE_OUT}.tsv
#Set a condition for the header situation
HEADER_WRITTEN=false
# Create an empty output file
> "$OUTPUT_FILE"
# Iterate over the file names and concatenate them
for ((i=1; i<=22; i++))
do
    INPUT_FILE=${BASE_OUT}_chr${i}.tsv
    # Check if it's the first file and write the header
    if [[ "$HEADER_WRITTEN" = false ]]; then
        cat "$INPUT_FILE" >> "$OUTPUT_FILE"
        HEADER_WRITTEN=true
    else
        # Skip the header and append the contents
        tail -n +2 "$INPUT_FILE" >> "$OUTPUT_FILE"
    fi
done
echo "Merging completed."

##Clean-up
rm $INTER
rm ${BASE_OUT}_chr*

