#!/usr/bin/env Rscript

#This script intends to intersect the filtered GWAS SNPs with the SNPs present in
#Relate. After that, it only retains the SNP with the highest p-value per LD block

#Load Packages
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))

#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
relate <- readRDS(files[2])

#Load the genetic correlation text file with the GWAS names
gwas_info <- read.delim(file = files[1], header = FALSE, sep = "\t")

#Define the path from where to load the GWAS_1 (with the latter included)
gwas_1_path <- paste0(files[3], "/", gwas_info[1,1])

#Load the GWAS_1 file
gwas_1 <- fread(file = gwas_1_path, header = TRUE, sep = "\t")

#Discard chromosome and position columns in GWAS 1 file (it is provided by relate, so we avoid duplicated)
gwas_1$chr <- NULL
gwas_1$pos <- NULL

#Remove useless columns of GWAS 1 from this step onwards
gwas_1$n_samples <- NULL
gwas_1$z <- NULL
gwas_1$rsID <- NULL

#Create objects with the names of the GWAS
gwas_1_name <- paste0("@",gwas_info[1,1])
gwas_2_name <- paste0("@",gwas_info[1,2])

#Intersect the selection statistics with the GWAS by common "variant" column
relate_gwas_1 <- inner_join(relate, gwas_1, by = "variant")

#Remove the other data
rm(relate, gwas_1)

#Define the path from where to load the GWAS_1 (with the latter included)
gwas_2_path <- paste0(files[3], "/", gwas_info[1,2])

#Load the GWAS_2 file
gwas_2 <- fread(file = gwas_2_path, header = TRUE, sep = "\t")

#Discard chromosome and position columns in GWAS 2 file (it is provided by relate, so we avoid duplicated)
gwas_2$chr <- NULL
gwas_2$pos <- NULL

#Remove useless columns of GWAS 2 from this step onwards
gwas_2$n_samples <- NULL
gwas_2$z <- NULL
gwas_2$rsID <- NULL

#Repeat intersection with the already intersected data frame
relate_gwas_1_2 <- inner_join(relate_gwas_1, gwas_2, by = "variant", suffix = c(gwas_1_name, gwas_2_name))

#Remove the other data
rm(relate_gwas_1, gwas_2)

##Allele validation

#Check that the ref and alt alleles of both GWAS match (or, at least, that they
#are flipped and they are not other alleles). Create new column "allele_info_1" 
#with 1 reference and alternate alleles match and -1 if "ref" is the other's "alt"
#and vice versa. NA is for any other result
relate_gwas_1_2$allele_info_1 <- ifelse(relate_gwas_1_2[[paste0("ref",gwas_1_name)]] == relate_gwas_1_2[[paste0("ref",gwas_2_name)]] &
                                        relate_gwas_1_2[[paste0("alt",gwas_1_name)]] == relate_gwas_1_2[[paste0("alt",gwas_2_name)]], 1,
                                      ifelse(relate_gwas_1_2[[paste0("ref",gwas_1_name)]] == relate_gwas_1_2[[paste0("alt",gwas_2_name)]] &
                                               relate_gwas_1_2[[paste0("alt",gwas_1_name)]] == relate_gwas_1_2[[paste0("ref",gwas_2_name)]], -1, NA))

# Remove SNPs with NA in "allele_info_1" column (if any)
relate_gwas_1_2 <- relate_gwas_1_2[!is.na(relate_gwas_1_2$allele_info_1), ]

# Flip the beta sign of the second GWAS
relate_gwas_1_2[[paste0("beta",gwas_2_name)]] <- relate_gwas_1_2[[paste0("beta",gwas_2_name)]] * relate_gwas_1_2$allele_info_1


#Create new column "allele_info_2" with 1 if "anc" is "ref" from GWAS 1 and "der"
#is "alt" for GWAS 1, -1 for "anc" is "alt" for GWAS 1 and "der" is "ref" for 
#GWAS 1, and NA for other cases. GWAS 2 doesn't matter because it has already been
#accounted in the previous "allele_info_1"
relate_gwas_1_2$allele_info_2 <- ifelse(relate_gwas_1_2$anc == relate_gwas_1_2[[paste0("ref",gwas_1_name)]] &
                                       relate_gwas_1_2$der == relate_gwas_1_2[[paste0("alt",gwas_1_name)]], 1,
                                     ifelse(relate_gwas_1_2$anc == relate_gwas_1_2[[paste0("alt",gwas_1_name)]] &
                                              relate_gwas_1_2$der == relate_gwas_1_2[[paste0("ref",gwas_1_name)]], -1, NA))

#Remove SNPs with NA in "allele_info_2" column (if any)
relate_gwas_1_2 <- relate_gwas_1_2[!is.na(relate_gwas_1_2$allele_info_2),]

#In this case, we don't flip the alleles as this task is handled by the PALM scripts.

#Remove "ref" and "alt" alleles from the second gwas
relate_gwas_1_2[[paste0("ref",gwas_2_name)]] <- NULL
relate_gwas_1_2[[paste0("alt",gwas_2_name)]] <- NULL

#Rename the ones in the first one as just "ref" and "alt" (no ...@foo.). This is
#needed so that palm.py can work well
colnames(relate_gwas_1_2)[colnames(relate_gwas_1_2) == paste0("ref",gwas_1_name)] <- "ref"
colnames(relate_gwas_1_2)[colnames(relate_gwas_1_2) == paste0("alt",gwas_1_name)] <- "alt"

#Remove allele_info columns (no longer used)
relate_gwas_1_2$allele_info_1 <- NULL
relate_gwas_1_2$allele_info_2 <- NULL

#SNP selection per LD block

##Create function to select the SNP(s) with lowest p-value in every LD block
select_SNPs_jpalm <- function(data, pval_gwas_1, pval_gwas_2) {
  #Create a vector with the bins that are unique to use for iteration
  ld_blocks <- sort(unique(data$LD_block))
  #Create an empty list for the indexes of the lowest p-vals per LD block
  snp_indexes_gwas_1 <- vector("list", length = length(ld_blocks))
  snp_indexes_gwas_2 <- vector("list", length = length(ld_blocks))
  #for loop to extract the indexes of the SNPs with lowest p-values
  for (i in seq_along(ld_blocks)) {
    #Find the lowest p-value per bin. This p_value may be the same for different SNPs
    best_pval_gwas_1 <- min(data[[pval_gwas_1]][data$LD_block==ld_blocks[i]])
    best_pval_gwas_2 <- min(data[[pval_gwas_2]][data$LD_block==ld_blocks[i]])
    #Use the previous value to obtain the indexes of the SNPs with that p-value
    pval_indexes_gwas_1 <- which(data[[pval_gwas_1]]==best_pval_gwas_1 & data$LD_block==ld_blocks[i])
    pval_indexes_gwas_2 <- which(data[[pval_gwas_2]]==best_pval_gwas_2 & data$LD_block==ld_blocks[i])
    #Store the indexes in the list
    snp_indexes_gwas_1[[i]] <- pval_indexes_gwas_1
    snp_indexes_gwas_2[[i]] <- pval_indexes_gwas_2
  }
  #Unlist the indexes to create a single vector
  snp_indexes_gwas_1 <- unlist(snp_indexes_gwas_1)
  snp_indexes_gwas_2 <- unlist(snp_indexes_gwas_2)
  #Create a new data frame with just the selected SNPs
  selected_snps_gwas_1_df <- data[snp_indexes_gwas_1,]
  selected_snps_gwas_2_df <- data[snp_indexes_gwas_2,]
  #Join data frames
  selected_snps_df <- rbind(selected_snps_gwas_1_df, selected_snps_gwas_2_df)
  #Remove potential duplicated SNPs
  selected_snps_df <- distinct(selected_snps_df)
  #Return the selected data frame
  return(selected_snps_df)
}

##Create necessary objects to call the function
pvals_gwas_1 <- paste0("pval",gwas_1_name)
pvals_gwas_2 <- paste0("pval",gwas_2_name)

##Call the function
selected_snps <- select_SNPs_jpalm(relate_gwas_1_2, pvals_gwas_1, pvals_gwas_2)

#Number of significant SNPs in independent LD blocks

##Create new set for significant SNPs for GWAS 1
significant <- selected_snps[selected_snps[[paste0("pval", gwas_1_name)]]<as.numeric(files[4]),]

##Check the numner of significant independent LD Blocks are in this selection
num_unique_ld_blocks <- length(unique(significant$LD_block))

##Add this info as a column for GWAS 1
selected_snps[[paste0("min_independent_significant_SNPs",gwas_1_name)]] <- num_unique_ld_blocks

##Repeat for GWAS 2
significant <- selected_snps[selected_snps[[paste0("pval", gwas_2_name)]]<as.numeric(files[4]),]
num_unique_ld_blocks <- length(unique(significant$LD_block))
selected_snps[[paste0("min_independent_significant_SNPs",gwas_2_name)]] <- num_unique_ld_blocks

#Create output file
output_file_name_1 <- paste0(gsub("\\.tsv$", "", gwas_info[1,1]), "-", gsub("\\.tsv$", "", gwas_info[1,2]), "_selected_SNPs.tsv")
write_tsv(selected_snps, file = output_file_name_1)