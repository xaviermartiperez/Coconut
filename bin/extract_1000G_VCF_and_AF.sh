#!/bin/bash

#This script extracts the VCFs of the subpopulations of the 1000G from its "global"
#VCF file with the help of the Sample IDs file. Then it creates a TSV file with
#The chr, pos, ref, alt and af columns.

#Load modules
module load BCFtools

#Define objects
IDs=$1
#IDs=/homes/users/xmarti/scratch/gcat/igsr-1000_genomes_phase_3_release.tsv
VCFDir=$2
#VCFDir=/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf
BASE_OUTDIR_VCF=${3}/1000G_VCFs
#BASE_OUTDIR_VCF=/homes/users/xmarti/scratch/gcat/relate/1000G_VCFs
POP=$4
#POP=FIN
CHR=$5
#CHR=22
OUTDIR_AF=${3}/${POP}/daf
#OUTDIR_AF=/homes/users/xmarti/scratch/gcat/relate/${POP}/daf

##SUB-population VCF extraction
if [ -e "${BASE_OUTDIR_VCF}/${POP}/${POP}_chr${CHR}.vcf.gz" ]; then
    #Create the report
    echo "${POP}_chr${CHR}.vcf.gz file already exists." > ${POP}_chr${CHR}_VCF_report.txt
else
    #Extract the sample IDs that correspond to the specified population
    awk -F'\t' -v subpop="$POP" '$4 == subpop {print $1}' $IDs > IDs.txt
    #Create directory if it doesn't exist already
    mkdir -p ${BASE_OUTDIR_VCF}/${POPULATION}
    #Call bcftools command to extract VCF with specified samples. Add a "." when no "alt is found
    #Don't specify an output name and pipe it to the next command. 
    bcftools view -S IDs.txt --force-samples --trim-alt-alleles -Ou \
    ${VCFDir}/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
    #Recalculate AF considering the AC of the specified subpopulation
    bcftools +fill-tags - -Oz -o ${BASE_OUTDIR_VCF}/${POP}/${POP}_chr${CHR}.vcf.gz -- -t AF,AN,AC,AC_Hemi,AC_Hom,AC_Het
    #Create the report
    echo "${POP}_chr${CHR}.vcf.gz file has been created." > ${POP}_chr${CHR}_VCF_report.txt
    #Clean-up
    rm IDs.txt
fi

##Alternate (NOT derived) allele frequency extraction
if [ -e "${OUTDIR_AF}/${POP}_chr${CHR}_af.tsv" ]; then   
    #Create the report
    echo "${POP}_chr${CHR}_af.tsv file already exists." > ${POP}_chr${CHR}_af_report.txt
else
    #Create directory if it doesn't exist already
    mkdir -p ${OUTDIR_AF}
    #Query the columns (and regions) of interest
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' \
    ${BASE_OUTDIR_VCF}/${POP}/${POP}_chr${CHR}.vcf.gz > ${OUTDIR_AF}/${POP}_chr${CHR}_af.tsv
    #Split the rows with multiple "alt" alleles into separate rows with one just allele and its frequency
    awk -F'\t' '{ if (index($4, ",") > 0) { alt_count = split($4, alts, ","); split($5, afs, ","); for (i = 1; i <= alt_count; i++) print $1 "\t" $2 "\t" $3 "\t" alts[i] "\t" afs[i]; } else print }' ${OUTDIR_AF}/${POP}_chr${CHR}_af.tsv > TMP.txt 
    #Remove rows with a "." in the alt column (AF=0), rows with AF=1 and the rows that are not SNPs
    awk -F'\t' '!(/\./ && ($4 == ".")) && ($5 != 1) && (length($3) == 1) && (length($4) == 1)' TMP.txt > ${OUTDIR_AF}/${POP}_chr${CHR}_af.tsv
    #Add column headers as first line (this command does work on Linux but not on MacOS)
    sed -i '1i chr\tpos\tref\talt\taf' ${OUTDIR_AF}/${POP}_chr${CHR}_af.tsv
    #Create the report
    echo "${POP}_chr${CHR}_af.tsv file has been created." > ${POP}_chr${CHR}_af_report.txt
    #Clean-up
    rm TMP.txt
fi
