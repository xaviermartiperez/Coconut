#!/usr/bin/env Rscript

#Suppress only automatic warnings from R
#options(warn = -1)

#Load Packages
require(dplyr)
require(stringr)
require(tidyverse)
require(stringi)
require(data.table)

#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
gwas_info <- read.delim(files[1], header = FALSE, sep = "\t")
gwas_dir <- files[2]
gwas_path <- paste0(gwas_dir, "/", gwas_info[1])
genome_version <- gwas_info[2]
effect_allele_is_ref <- gwas_info[3]
n_sample <- as.numeric(gwas_info[4])
hapmap3_hg38 <- files[3]
hapmap3_hg19 <- files[4]
if (genome_version == "hg38") {
  hapmap3 <- hapmap3_hg38
} else if (genome_version == "hg19") {
  hapmap3 <- hapmap3_hg19
} else {
  stop("Invalid genome version. Please use 'hg38' or 'hg19'.")
}
# Define name of the report output file
output_name <- paste0("formatted_", sub("\\..*", "", gwas_info[1]), ".tsv")

# Define final output name for formatted GWAS
final_output_name <- paste0("formatted_", sub("\\..*", "", gwas_info[1]), ".tsv")
final_output_path <- paste0(gwas_dir, "/", final_output_name)
other_final_output_path <- paste0(gwas_dir, "/", "fixed_", sub("\\..*", "", gwas_info[1]), ".tsv")

if (file.exists(final_output_path) || file.exists(other_final_output_path)) {
  #If the file already exists, emit the expected output
  final_output_content <- paste(output_name, genome_version, FALSE, FALSE, sep = "\t")
  writeLines(final_output_content, sub("\\.tsv$", ".txt", output_name))
} else {
  #Function to load GWAS according extension
  load_GWAS <- function(data) {
    extension <- sub(".*\\.(.+)\\.(gz|bgz)$|.*\\.(.+)$", "\\1\\3", basename(data))
    switch(extension,
           "tsv" = fread(data, header = TRUE, sep = "\t"),
           "csv" = fread(data, header = TRUE, sep = ","),
           "txt" = fread(data, header = TRUE),
           "regenie" = fread(data, header = TRUE),
           "vcf" = {
             vcf <- read.vcfR(data)
             #We extract only the information of interest
             data <- as.data.frame(vcf@fix)
             format <- as.data.frame(extract_gt_tidy(vcf, format_fields = c("ES", "SE", "LP","SS"), 
                                                     gt_column_prepend = "", alleles = FALSE))
             cbind(data, format)[, -c(6:10)]
           }
    )
  }
  
  #Load GWAS
  gwas <- load_GWAS(gwas_path)
  print("Finished loading GWAS")
  #Function to change header
  change_headers <- function(data, gwas_name) {
    #Check if it is harmonized
    harm <- grepl("\\.h\\.tsv\\.gz$", gwas_name)
    hm_var <- ifelse(harm, grepl("^[0-9]+-GCST[0-9]+-.*\\.h\\.tsv\\.gz$", gwas_name), FALSE)
    
    #If the file name presents a determined structure, select the harmonized variables together with the standard error and p-value
    if (harm == TRUE) {
      if (hm_var == TRUE) {
        selected_col <- grepl("^hm_", colnames(data)) | colnames(data) %in% c("standard_error", "p_value")
        data <- data[, ..selected_col]
        #Remove the "hm_" prefix of the columns' name
        colnames(data) <- sub("^hm_", "", colnames(data)) 
      } 
    }
    
    #Possible names for the header
    variant <- c("^(variant|uniqID|variant_ID|markername|SNPID|ID|SNP|cptid)$")
    beta <- c("^(beta|effect|effects|b|effect_size|es)$")
    se <- c("^(se|standard_error|StandardError|StdErr|sebeta)$")
    pval <- c("^(pval|p_val|p-val||p.val|pvalue|p_value|p-value|p.value|gc_pvalue|p)$")
    LP <- c("^(lp|neg_log_10_p_value|pvalue_mlog|mlogp|log10p)$")
    chr <- c("^(chr|chromosome|chrom|snp-chr|chr_ID|#chrom)$")
    pos <- c("^(pos|genpos|base_pair_location|snppos|snp_pos|bp|position|chr_pos|pos\\(b37\\))$")
    z <- c("^(z|zscore|z-score|GC_zscore|tstat|t_stat|t-statistic)$")
    n_samples <- c("^(N|n_complete_samples|sample_size|n_total|TotalSampleSize|SS)$")
    ref <- c("^(ref|reference_allele|ref_allele)$")
    alt <- c("^(alt|alternate_allele|alt_allele)$")
    effect_allele <- c("^(effect_allele|A1|allele1|allele_1|inc_allele|EA)$")
    non_effect_allele <- c("^(other_allele|A2|allele2|allele_2|non_effect_allele|dec_allele|NEA|allele0)$")
    odd_ratio <- c("^(OR|odd_ratio)$")
    
    #If a column presents rsID, change the column name for rsID
    colnames(data) <- ifelse(apply(head(data, 50), 2, function(col) { 
      any(stri_detect_regex(col, "^rs", case_insensitive = TRUE))
    }), "rsID", colnames(data))
    
    #Change column names to the stablished format
    colnames(data) <- stri_replace_all_regex(colnames(data),
                                             pattern = c(variant, beta, se, z, ref, alt,
                                                         pval, LP, chr, pos, n_samples,
                                                         effect_allele, non_effect_allele, odd_ratio),
                                             replacement = c("variant", "beta", "se", "z", "ref", "alt",
                                                             "pval", "LP", "chr", "pos", "n_samples",
                                                             "effect_allele", "non_effect_allele", "odd_ratio"),
                                             vectorize = FALSE, case_insensitive = TRUE)
    return(data)
  }
  
  #Call the function to change headers
  gwas <- change_headers(gwas, gwas_info[1])
  
  #If it is necessary transform odd ratio to beta
  if (!("beta" %in% colnames(gwas)) && "odd_ratio" %in% colnames(gwas)) {
    gwas$odd_ratio <- as.numeric(gwas$odd_ratio)
    gwas$beta <- log(gwas$odds_ratio)
  }
  
  #Determine beta as a numeric variable
  gwas$beta <- as.numeric(gwas$beta)
  
  #If it is necessary, transform -log10(p-value) to p-value
  if ("LP" %in% colnames(gwas) && !("pval" %in% colnames(gwas))) {
    gwas$LP <- as.numeric(gwas$LP)
    gwas$pval <- 10^(-gwas$LP)
  }
  print(head(gwas$LP, 10))
  print(head(gwas$pval, 10))
  #If it is necessary, calculate SE
  if (!("se" %in% colnames(gwas))) {
    if ("z" %in% colnames(gwas)) {
      gwas$z <- as.numeric(gwas$z)
      gwas$se <- gwas$beta / gwas$z
    } else {
      gwas$se <- 0
    }
  }
  
  #If it is necessary, calculate Z
  if (!("z" %in% colnames(gwas)) && "se" %in% colnames(gwas) && "se" != 0) {
    gwas$se <- as.numeric(gwas$se)
    gwas$z <- gwas$beta / gwas$se
  }
  
  #If it is necessary, add number of samples
  if (!("n_samples" %in% colnames(gwas))) {
    gwas$n_samples = n_sample
  }
  
  # Function to change the header names in allele columns
  alleles_f <- function(data, effect_allele_is_ref) {
    # If NA or "no"/"n", return data as-is
    if (is.na(effect_allele_is_ref) || tolower(effect_allele_is_ref) %in% c("no", "n")) {
      if ("effect_allele" %in% colnames(data)) {
        colnames(data)[match("effect_allele", colnames(data))] <- "alt"
        colnames(data)[match("non_effect_allele", colnames(data))] <- "ref"
      }
      return(data)
      # If it "yes" or "y", flip the beta sign
    } else if (tolower(effect_allele_is_ref) %in% c("yes", "y")) {
      data$beta <- -data$beta
      if ("effect_allele" %in% colnames(data)) {
        colnames(data)[match("effect_allele", colnames(data))] <- "ref"
        colnames(data)[match("non_effect_allele", colnames(data))] <- "alt"
      }
      return(data)
    } else {
      stop("Invalid value for effect_allele_is_ref. Please leave blank, 'yes', 'no', 'y', or 'n'.")
    }
  }
  
  #Determine Reference and Alternative alleles
  ##There are GWAS  files where there is effect/non_effect columns and a ref allele column, therefore based on
  ##these columns we extract ref and alt allele columns.
  if ("effect_allele" %in% colnames(gwas) & ("ref" %in% colnames(gwas))) {
    gwas$ref <- ifelse(
      gwas$ref == "EA", gwas$effect_allele,
      ifelse(gwas$ref == "OA", gwas$non_effect_allele,
             ifelse(gwas$ref %in% c("A", "C", "G", "T"), gwas$ref, NA)))
    
    gwas$alt <- ifelse(
      gwas$ref == "EA", gwas$non_effect_allele,
      ifelse(gwas$ref == "OA", gwas$effect_allele,
             ifelse(gwas$ref %in% c("A", "C", "G", "T"),
                    ifelse(gwas$effect_allele == gwas$ref, gwas$non_effect_allele, gwas$effect_allele), NA)))
    
    #Assign beta to the ALT allele
    gwas$beta <- ifelse(gwas$ref == gwas$effect_allele, -gwas$beta, gwas$beta)
    
    ##For the rest of cases we call the function to change names.
  } else {
    gwas <- alleles_f(gwas, effect_allele_is_ref)
  }
  
  #Create function to add missing columns from variant column or from the function previously created for rsID
  #It also discards SNPs from chr X, Y, MT and variants that are not 1 substitution base, and creates a variant column. 
  variant_function <- function(gwas) {
    #Determine if there is missing columns about chromosome, position or alleles.
    columns_variants <- c("chr", "pos", "ref", "alt")
    missing_columns <- setdiff(columns_variants, colnames(gwas))
    
    # Initialize variants (it is TRUE when the information needs to be extracted through VCF files)
    missing_ref_alt <- FALSE
    missing_chr_pos <- FALSE
    
    #Replace any separator (-, /, or _) in the variant variable with ":"
    if("variant" %in% colnames(gwas)) {
      gwas$variant <- gsub("[-/_]", ":", gwas$variant)
    }
    
    #Add missing columns
    if (length(missing_columns) != 0) {
      #From variant column structured as chr:pos:ref:alt
      if ("variant" %in% colnames(gwas) && all(grepl("^([0-9]+|X|Y|MT|M):[0-9]+:[A-Za-z]+:[A-Za-z]+$", head(gwas$variant, 10)))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 4)
        if (!"chr" %in% colnames(gwas)) gwas$chr <- split_variant[, 1]
        if (!"pos" %in% colnames(gwas)) gwas$pos <- split_variant[, 2]
        if (!"ref" %in% colnames(gwas)) gwas$ref <- split_variant[, 3]
        if (!"alt" %in% colnames(gwas)) gwas$alt <- split_variant[, 4]
        
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
        
        if ("effect_allele" %in% colnames(gwas) && is.null(effect_allele_is_ref)) {
          gwas$beta <- ifelse(gwas$ref == gwas$effect_allele, -gwas$beta, gwas$beta)
        }
        
        #From variant column structured as chr:pos
      } else if ("variant" %in% colnames(gwas) && any(grepl("^[0-9]+:[0-9]+$", head(gwas$variant, 10))) && identical(missing_columns, c("chr", "pos"))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 2)
        if (!"chr" %in% colnames(gwas)) gwas$chr <- split_variant[, 1]
        if (!"pos" %in% colnames(gwas)) gwas$pos <- split_variant[, 2]
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
        #From rsID
      } else if ("rsID" %in% colnames(gwas)) {
        missing_ref_alt <- any(c("ref", "alt") %in% missing_columns)
        missing_chr_pos <- any(c("chr", "pos") %in% missing_columns)
      }
    }
    
    if(missing_chr_pos == FALSE) {
      #Determine position and chromosome as numeric variables
      gwas$pos <- as.numeric(gwas$pos)
      #Discard SNPs from chr X, Y, MT and variants that are not 1 substitution base
      gwas$chr <- as.numeric(gwas$chr)
      gwas <- gwas[!is.na(gwas$chr), ]
      
      #If the variant still maintains the structure "chr:pos:ref:alt", remove the allele information (:ref:alt) to retain only the "chr:pos" structure.
      if ("variant" %in% colnames(gwas)) {
        # Remove "chr" prefix if present
        gwas$variant <- sub("^chr", "", gwas$variant)
        # Remove allele information if present
        gwas$variant <- ifelse(
          grepl("^(1[0-9]|2[0-2]|[1-9]):[0-9]+:[A-Za-z]+:[A-Za-z]+$", gwas$variant), 
          sub(":([A-Za-z]+:[A-Za-z]+)$", "", gwas$variant), gwas$variant)
      } else {
        #Add column variant if it is necessary
        gwas$variant <- paste0(gwas$chr, ":", gwas$pos)
      }
      
      #Determine variant as a string
      gwas$variant <- as.character(gwas$variant)
      
    } else {
      gwas$chr <- NA
      gwas$pos <- NA
      gwas$variant <- NA
    }
    
    if(missing_ref_alt==FALSE) {
      #Determine ref and alt columns as string and select variants with 1 substitution base.
      gwas$ref <- as.character(gwas$ref)
      gwas$alt <- as.character(gwas$alt)
      gwas <- gwas[nchar(gwas$ref) == 1 & nchar(gwas$alt) == 1, ]
    } else {
      gwas$effect_allele <- as.character(gwas$effect_allele)
      gwas$non_effect_allele <- as.character(gwas$non_effect_allele)
      gwas <- gwas[nchar(gwas$effect_allele) == 1 & nchar(gwas$non_effect_allele) == 1, ]
      gwas$ref <- NA
      gwas$alt <- NA
    }
    return(list(gwas = gwas, missing_ref_alt = missing_ref_alt, missing_chr_pos = missing_chr_pos))
  }
  
  #Call the function
  results <- variant_function(gwas)
  gwas <- results$gwas
  missing_ref_alt <- results$missing_ref_alt
  missing_chr_pos <- results$missing_chr_pos
  
  #Function to add rsID from chromosome and position column based on HapMap3
  add_rsID = function(data){
    if (!("rsID" %in% colnames(data))) {
      HM3 <- fread(hapmap3)
      data <- data %>%
        left_join(HM3, by = c("variant"))
    }
    return(data)
  }
  
  #Call the function
  gwas <- add_rsID(gwas)
  
  #Verify that the GWAS includes all the columns required for the downstream analysis.
  columns_req <- c("rsID", "variant", "n_samples", "beta", "se", 
                  "pval", "z", "ref", "alt", "chr", "pos")
  missing_columns <- setdiff(columns_req, colnames(gwas))
  
  if (length(missing_columns) > 0) {
    message("There is one or more columns missing, revise ", gwas_info[1], " file to verify it. In case there is no missing column, please adjust the header to match the script requirements.")
    stop("The following columns are missing: ", paste(missing_columns, collapse = ", "))
  } else {
    #Create GWAS output with the established format
    if(missing_ref_alt==TRUE) {
      formatted_gwas <- gwas[, c("variant", "n_samples", "beta", "se", "pval",
                                 "z", "ref", "alt", "chr", "pos", "rsID", 
                                 "effect_allele", "non_effect_allele")]
    } else {
      formatted_gwas <- gwas[, c("variant", "n_samples", "beta", "se", "pval",
                                 "z", "ref", "alt", "chr", "pos", "rsID")]      
    }
  }
  
  #Remove unnecessary data
  rm(gwas)
  
  #Create output file
  write_tsv(formatted_gwas, file = paste0(gwas_dir, "/", output_name))
  
  #Pass the GWAS name (ended in .tsv) to the next process in a .txt file
  output_content <- paste(output_name, genome_version, missing_ref_alt, missing_chr_pos, sep = "\t")
  writeLines(output_content, sub("\\.tsv$", ".txt", output_name))
}