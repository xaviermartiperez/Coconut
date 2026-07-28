#!/usr/bin/env Rscript

# Load Packages
require(dplyr)
require(stringr)
require(tidyverse)
require(data.table)

# Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)

# Load GWAS and LD blocks
gwas_info <- read.delim(files[1], header = FALSE)
gwasdir <- files[2]
gwas_path <- paste0(gwasdir, "/", gwas_info[1])
gwas <- fread(gwas_path, header = TRUE, sep = "\t")
ld_hg38 <- files[3]
ld_hg19 <- files[4]
if (gwas_info[2] == "hg38") {
  ld <- fread(ld_hg38)
} else if (gwas_info[2] == "hg19") {
  ld <- fread(ld_hg19)
} else {
  stop("Invalid genome version. Please use 'hg38' or 'hg19'.")
}
max_pval <- as.numeric(files[5])

# Remove SNPs above max. p-value
gwas <- gwas[gwas[["pval"]] < max_pval, ]

#Function to add LD blocks
add_LD_blocks <- function(data) {
  #Remove the "chr" part from the chr column
  ld$chr <- as.numeric(gsub("chr", "", ld$chr))
  
  #Create a column with the number of the LD block (the row number)
  ld$LD_block <- seq_len(nrow(ld))
  
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
gwas <- add_LD_blocks(gwas)

#Find number of LD blocks that have at least 1 significant SNP
result <- as.character(length(unique(gwas$LD_block)))

#Create output file and save just the name of the GWAS into this file
output_name <- paste0(sub("\\..*", "", gwas_info[1]), "-+", result)
writeLines(paste(as.character(gwas_info[1]), as.character(gwas_info[2]), sep = "\t"), output_name)

#Remove unnecessary data
rm(gwas)
rm(ld)