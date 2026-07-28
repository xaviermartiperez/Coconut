#!/usr/bin/env Rscript

#This script intends to intersect the filtered GWAS SNPs with the SNPs present in
#Relate. After that, it only retains the SNP with the highest p-value per LD block

#Suppress only automatic warnings from R
options(warn = -1)

# Packages
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))

# Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
relate <- readRDS(files[1])
gwas_name <- read.delim(files[2], header = FALSE)
gwasdir <- files[3]
gwas_path <- paste0(gwasdir, "/", gwas_name[1])
gwas <- fread(file = gwas_path, header = TRUE, sep = "\t")

#Discard chromosome and position columns in GWAS file (it is provided by relate, so we avoid duplicated)
gwas$chr <- NULL
gwas$pos <- NULL

#Remove useless columns of GWAS 1 from this step onwards
gwas$n_samples <- NULL
gwas$z <- NULL
gwas$rsID <- NULL

# Intersect the selection statistics with the GWAS by common "variant" column
relate_gwas <- inner_join(relate, gwas, by = "variant")

# Remove the other data
rm(relate, gwas)

# Allele validation
# Create new column "allele_info" with 1 if "anc" is "ref" and "der" is "alt",
# -1 for "anc" is "alt" and "der" is "ref", and NA for other cases. Besides a filter,
# this is also a good check to check that alleles in the Relate analysis have been
# flipped according to ancestry (that is, they are not all 1)
relate_gwas$allele_info <- ifelse(relate_gwas$anc == relate_gwas$ref &
                                       relate_gwas$der == relate_gwas$alt, 1,
                                     ifelse(relate_gwas$anc == relate_gwas$alt &
                                              relate_gwas$der == relate_gwas$ref, -1, NA))

# Remove SNPs with NA in "allele_info" column (if any)
relate_gwas <- relate_gwas[!is.na(relate_gwas$allele_info), ]

# SNP selection per LD block

## Function to select the SNP(s) with lowest p-value in every LD block
select_SNPs <- function(data) {
  #Create a vector with the bins that are unique to use for iteration
  ld_blocks <- sort(unique(data$LD_block))
  #Create an empty list for the indexes of the lowest p-vals per LD block
  snp_indexes <- vector("list", length = length(ld_blocks))
  #for loop to extract the indexes of the SNPs with lowest p-values
  for (i in seq_along(ld_blocks)) {
    #Find the lowest p-value per bin. This p_value may be the same for different SNPs
    best_pval <- min(data$pval[data$LD_block==ld_blocks[i]])
    #Use the previous value to obtain the indexes of the SNPs with that p-value
    pval_indexes <- which(data$pval==best_pval & data$LD_block==ld_blocks[i])
    #Store the indexes in the list
    snp_indexes[[i]] <- pval_indexes
  }
  #Unlist the indexes to create a single vector
  snp_indexes <- unlist(snp_indexes)
  #Create a new data frame with just the selected SNPs
  selected_snps_df <- data[snp_indexes, ]
  #Return the selected data frame
  return(selected_snps_df)
}

## Call the function
selected_snps <- select_SNPs(relate_gwas)

# Create output file
output_file_name <- paste0(gsub("\\.tsv$", "", gwas_name), "_selected_SNPs.tsv")
write_tsv(selected_snps, file = output_file_name)
