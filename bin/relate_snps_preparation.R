#!/usr/bin/env Rscript

#This script intends to process the Relate .mut files so that the SNPs that it
#contains (which are already biallelic and have the ancestral/derived information)
#can be intersected with the GWAS SNPs . It also has the purpose of adding the LD
#block information to each SNP (as inferred in Berisa and Pickrell, 2015).


#####Preparation#####

#Packages
require(dplyr)
require(stringr)
require(tidyverse)

#Working directory
#setwd("/Users/xavi_marti/Documents/academia/phd/data/gcat/ukbb/relate")
#setwd("/homes/users/xmarti/scratch/gcat")

##Load objects
#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
#Read the population
population <- as.character(files[1])
#population <- "FIN"
#Read the Relate directory (base directory path, needs to be refined later on)
relate_dir <- as.character(files[2])
#relate_dir <- "/homes/users/xmarti/scratch/gcat/relate"
#Read the output directory (the one where you want the files to be saved to)
outdir <- as.character(files[3])
#outdir <- "/homes/users/xmarti/scratch/gcat"
#Read LD blocks
ld <- read.delim(file = files[4], header = TRUE, sep = "\t")
#ld <- read.delim(file = "ld_blocks.bed", header = TRUE, sep = "\t")

#####Data processing#####

##Create a single data frame

#Create a function to do the unification of the files into one data frame
unify_relate_chromosomes <- function(path_to_files) {
  #Define the chromosome list
  chromosomes <- 1:22
  #Create empty vector to store the individual chromosome files
  df_list <- vector("list", length = length(chromosomes))
  #Loop to store the chromosome files into the list
  for (i in seq_along(chromosomes)) {
    #Load the individual file 
    chr_df <- read.delim(paste0(path_to_files, "/", population, "/anc_mut/",
     population, "_chr", as.character(chromosomes[i]), "_upd.mut"), sep = ";", header = TRUE)
    #Add chromosome column
    chr_df$chr <- as.numeric(chromosomes[i])
    #Store the file in the list
    df_list[[i]] <- chr_df
  }
  #Merge all the chromosomes into the unified data frame
  merged_df <- do.call(rbind, df_list)
  #Return the unified data frame
  return(merged_df)
}

#Call the function
relate_df <- unify_relate_chromosomes(relate_dir)

##Create necessary columns and remove useless ones

#Rename the position column for standardisation purposes
relate_df <- relate_df %>% rename(pos = pos_of_snp)
#Create variant column, necessary for the intersection
relate_df$variant <- paste0(relate_df$chr, ":", relate_df$pos)
#Create "anc" and "der" columns from the one that has both alleles
relate_df <- separate(relate_df, ancestral_allele.alternative_allele, into = c("anc", "der"), sep = "/")

#Remove useless columns
relate_df$snp <- NULL
relate_df$dist <- NULL
relate_df$`rs.id` <- NULL
relate_df$tree_index <- NULL
relate_df$branch_indices <- NULL
relate_df$is_not_mapping <- NULL
relate_df$is_flipped <- NULL
relate_df$age_begin <- NULL
relate_df$age_end <- NULL
relate_df$upstream_allele <- NULL
relate_df$downstream_allele <- NULL
relate_df$X <- NULL
relate_df[[population]] <- NULL

##Addition of the LD block information

#Create a function to add the LD block (it can be thought of as a bin) to the
#Relate data frame
add_LD_blocks <- function(data) {
  #Remove the "chr" part from the chr column
  ld$chr <- as.numeric(gsub("chr", "", ld$chr))
  #Create a column with the number of the LD block (the row number)
  ld$LD_block <- row_number(ld$chr)
  #Define chromosomes
  chromosomes <- 1:22
  #Create for loop to iterate over chromosomes
  for (i in chromosomes) {
    start_positions <- ld$start[ld$chr==i]
    data$LD_start[data$chr==i] <-
      start_positions[(findInterval(data$pos[data$chr==i],
                                   ld$stop[ld$chr==i] - 1e-10) + 1)] #The +1 is added because otherwise it selects the previous interval
    #Match the positions to obtain the LD block
    data$LD_block[data$chr==i] <-
      ld$LD_block[match(data$LD_start[data$chr==i], ld$start)]
  }
  #Remove the now useless LD_start column
  data$LD_start <- NULL
  #Return the data frame with the added column of interest
  return(data)
}
#Call the function
relate_df <- add_LD_blocks(relate_df)

##Add DAF of GCAT SNPs. NOTE: This step required a previous step done in the
##script "extract_relate_daf.sh", in the cluster

#Create function to check that the alleles match and to add the DAF
add_relate_DAF <- function(data, path_to_files) {
  #Define chromosomes
  chromosomes <- 1:22
  #Create for loop to add the DAF per chromosome
  for (i in chromosomes) {
    #Load the file with the DAF and the alleles per chromosome
    daf_ch <- read.delim(paste0(path_to_files, "/", population, "/daf/",
                                population, "_chr", i, "_af.tsv"), 
                          header = TRUE, sep = "\t")
    #Create "variant" column for intersection by posotion
    daf_ch$variant <- paste0(daf_ch$chr, ":", daf_ch$pos)
    #Add the DAF to the data frame of interest
    data$AF[data$chr==i] <- daf_ch$af[match(data$variant[data$chr==i], daf_ch$variant)]
    #Add the reference allele to the data frame of interest
    data$ref[data$chr==i] <- daf_ch$ref[match(data$variant[data$chr==i], daf_ch$variant)]
    #Add the alternate allele to the data frame of interest
    data$alt[data$chr==i] <- daf_ch$alt[match(data$variant[data$chr==i], daf_ch$variant)]
  }
  #Set 0 when "ref" is "anc" and "alt" is "der"
  data$af_modifier[data$ref == data$anc & data$alt == data$der] <- 0
  #Set 1 when "ref" is "der" and "alt" is "anc"
  data$af_modifier[data$ref == data$der & data$alt == data$anc] <- 1
  #Remove rows with NA values
  data <- data[!is.na(data$AF),]
  #Create DAF column
  data$AF <- as.numeric(data$AF)
  data$DAF <- abs(data$af_modifier - data$AF)
  #Remove useless columns
  data$AF <- NULL
  data$ref <- NULL
  data$alt <- NULL
  data$af_modifier <- NULL
  #Return the data frame with the new information
  return(data)
}
#Call the function
relate_df <- add_relate_DAF(relate_df, relate_dir)

#####Save the result#####
output_path <- paste0(outdir, "/relate_snps_", population, ".rds")
saveRDS(relate_df, file = output_path)

#Save just the name of the GWAS into this file
report_file_name <- paste0("relate_snps_", population, "_report.txt")
output_name <- paste0("relate_snps_", population, ".rds")
writeLines(paste0(output_name, " file has been created."), report_file_name)


