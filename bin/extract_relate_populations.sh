#!/bin/bash

#This script extracts the anc/mut files for the specified subpopulation and
#chromosome from the joint Relate anc/mut files provided by Spiedel. It is
#advisble (according to the software) to have annotated the .mut files before
#doing this step.

#Load modules
module load relate/1.1.6-x86_64
module load GCC/12.3.0

#Define objects
PATH_TO_RELATE=/homes/users/xmarti/scratch/relate
BASE_RELATE_DIR=$3
#BASE_RELATE_DIR=/homes/users/xmarti/scratch/gcat/relate
INDIR=${BASE_RELATE_DIR}/ALL_RESULTS
POPULATION=$1
#POPULATION=IBS
CHR=$2
#CHR=13


if [ -e "${BASE_RELATE_DIR}/${POPULATION}/anc_mut/${POPULATION}_chr${CHR}_upd.anc" ]; then
    echo "${POPULATION}_chr${CHR}_upd.* files already exist." > ${POPULATION}_chr${CHR}_extract_report.txt
else
    mkdir -p ${BASE_RELATE_DIR}/${POPULATION}/anc_mut
    #Call Relate script to extract subtrees corresponding to subpopulations
    $PATH_TO_RELATE/bin/RelateExtract \
    --mode SubTreesForSubpopulation \
    --anc ${INDIR}/1000GP_Phase3_mask_chr${CHR}.anc.gz \
    --mut ${INDIR}/1000GP_Phase3_mask_chr${CHR}.mut.gz \
    --poplabels ${INDIR}/1000GP_Phase3_mask_chr${CHR}.poplabels \
    --pop_of_interest ${POPULATION} \
    -o ${BASE_RELATE_DIR}/${POPULATION}/anc_mut/${POPULATION}_chr${CHR}_upd
    echo "${POPULATION}_chr${CHR}_upd.* files have been created." > ${POPULATION}_chr${CHR}_extract_report.txt
fi