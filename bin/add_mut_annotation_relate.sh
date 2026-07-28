#!/bin/bash

#This script annotates the .mut files by overwriting them, so it is advisable that these
#files are not the original ones, but a copy. Anc and Mut raw files were previously renamed:
#apparently the good ones (after refinement) contained "prene" in the name of the file. I
#essentially renamed the "prene" files taking out the "prene", thus also removing the other
#anc/mut files. Remember that this is a copy so it's not as if I'm deleting original data.
#The script to rename the files is called "rename_relate_anc_mut_files.sh"

#Chromosome 13 doesn't have a mut file, it's in the "mut.zip". I renamed this file with the
#command line

#Load modules
module load relate/1.1.6-x86_64
module load GCC/12.3.0

#Define paths
PATH_TO_RELATE=/homes/users/xmarti/scratch/relate
#INDIR=${1}/ALL_RESULTS
INDIR=/homes/users/xmarti/scratch/gcat/relate/ALL_RESULTS
#Define chromosome
CHR=$2
#CHR=13

#Call Relate script to add population annotation to .mut.gz files
$PATH_TO_RELATE/bin/RelateFileFormats \
--mode GenerateSNPAnnotations \
--haps ${INDIR}/1000GP_Phase3_mask_chr${CHR}.haps.gz \
--sample ${INDIR}/1000GP_Phase3_mask_chr${CHR}.sample.gz \
--poplabels ${INDIR}/1000GP_Phase3_mask_chr${CHR}.poplabels \
-o ${INDIR}/1000GP_Phase3_mask_chr${CHR} \
--mut ${INDIR}/1000GP_Phase3_mask_chr${CHR}.mut.gz

#The previous step creates a .mut file instead of overwriting the .mut.gz,
#For simplicity, replace the .mut.gz for the .mut (by removing the .mut.gz and
#compressing the file)
rm ${INDIR}/1000GP_Phase3_mask_chr${CHR}.mut.gz
gzip ${INDIR}/1000GP_Phase3_mask_chr${CHR}.mut

echo "Done!"