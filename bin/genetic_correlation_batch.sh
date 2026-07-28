#!/bin/bash

#This script implements LDSC (2 to 3 version) to obtain the genetic
#correlation between two GWAS traits

#Load Miniconda module to load environment
module load Miniconda3/4.9.2

#Define project and output directory
prDir=$3
PATH_TO_LOGS="$prDir/rg"

#Create directory and if it already exists, the existing directory is used.
mkdir -p $PATH_TO_LOGS


#Load ldsc environment (this automatically loads 2-to-3 branch)
source activate ${prDir}/conda_env/ldsc

#Set paths
PATH_TO_LDSC=${prDir}/ldsc
PATH_TO_SUMSTATS=${prDir}/jpalm_gwas_sumstats
PATH_TO_LD_SCORES=${prDir}/gcat/ldscore/eur_w_ld_chr/
#PATH_TO_LD_SCORES=/homes/users/ajimenez/scratch/gcat/ldscore/eur_w_ld_chr/
#PATH_TO_LOGS="/homes/users/ajimenez/scratch/rg"

#Read values from the GWAS pair batch file (it ends in .tsv) and assign them to variables
while IFS=$'\t' read -r GWAS_1 GWAS_2
do
  #Critical step: check if genetic correlation file already exists
  if [ -e "${PATH_TO_LOGS}/${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg.log" ]; then
    #Store the genetic correlation value and the p-value in a file
    awk '/Genetic Correlation:/ {genetic_corr=$3} /P:/ {print genetic_corr "\t" $2}' ${PATH_TO_LOGS}/${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg.log > pre_result.txt
    
    #Read the values from "pre_result.txt"
    read -r VALUE_1 VALUE_2 < pre_result.txt
    
    #Check if VALUE_1 is "nan" and set it to 0.
    if [ "$VALUE_1" == "nan" ]; then
      VALUE_1=0
    fi
    
    #Read the number of comparisons
    MULTIPLIER=$(cat $2)
    PRE_RESULT=$(awk -v x="$VALUE_2" -v y="$MULTIPLIER" 'BEGIN {printf "%.5e", x * y}')
    
    #Add a leading zero before the decimal point (just in case Groovy doesn't read well Xe-X numbers)
    RESULT=$(printf "%.5f" "$PRE_RESULT")
    
    #Add the GWAS to the file, in the beginning of it. Modify the extension from .sumstats.gz to .tsv
    echo -e "$GWAS_1\t$GWAS_2\t$VALUE_1\t$RESULT" >> ${GWAS_1%.tsv}-+${GWAS_2%.tsv}-+${VALUE_1}-+${RESULT}
    
    #Remove now-undesired files
    rm pre_result.txt
  else
    #Call the function. Make sure to change the extension of the GWAS and load them from the
    #sumstats directory. Set a timeout because apparently it sometimes gets stuck
    timeout 900 python $PATH_TO_LDSC/ldsc.py \
    --rg $PATH_TO_SUMSTATS/${GWAS_1%.tsv}.sumstats.gz,$PATH_TO_SUMSTATS/${GWAS_2%.tsv}.sumstats.gz \
    --ref-ld-chr $PATH_TO_LD_SCORES \
    --w-ld-chr $PATH_TO_LD_SCORES \
    --out ${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg

    #Check if the timeout occurred
    if [ $? -eq 124 ]; then
      #Timeout occurred, set rg = 0  and p-value = 1
      echo -e "0\t1" > pre_result.txt
    else
      #Timeout did not occur, continue with the rest of the script
      #Store the genetic correlation value and the p-value in a file
      awk '/Genetic Correlation:/ {genetic_corr=$3} /P:/ {print genetic_corr "\t" $2}' ${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg.log > pre_result.txt
    fi

    #Read the values from "pre_result.txt"
    read -r VALUE_1 VALUE_2 < pre_result.txt

    #Check if VALUE_1 is "nan" and set it to 0.
    if [ "$VALUE_1" == "nan" ]; then
      VALUE_1=0
    fi
    
    #Read the number of comparisons
    MULTIPLIER=$(cat $2)
    PRE_RESULT=$(awk -v x="$VALUE_2" -v y="$MULTIPLIER" 'BEGIN {printf "%.5e", x * y}')

    #Add a leading zero before the decimal point (just in case Groovy doesn't read well Xe-X numbers)
    RESULT=$(printf "%.5f" "$PRE_RESULT")

    #Add the GWAS to the file, in the beginning of it. Modify the extension from .sumstats.gz to .tsv
    echo -e "$GWAS_1\t$GWAS_2\t$VALUE_1\t$RESULT" >> ${GWAS_1%.tsv}-+${GWAS_2%.tsv}-+${VALUE_1}-+${RESULT}

    #Copy the log in the logs path
    cp ${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg.log "${PATH_TO_LOGS}/."
    
    #Remove now-undesired files
    rm pre_result.txt
    rm ${GWAS_1%.tsv}_${GWAS_2%.tsv}_rg.log
  fi
done < "$1"
