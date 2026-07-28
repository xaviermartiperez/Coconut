#!/usr/bin/env Rscript

#Suppress only automatic warnings from R
options(warn = -1)

#Load Packages
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(tidyverse))
suppressMessages(require(stringi))
suppressMessages(require(vcfR))
suppressMessages(require(data.table))

#Put the arguments passed to the script into a vector
files <- commandArgs(trailingOnly = TRUE)
gwas_info <- read.delim(files[1], header = FALSE)
gwas_dir <- files[2]
gwas_path <- paste0(gwas_dir, "/", gwas_info[1])
vcf_dir <- files[3]
n_sample <- as.numeric(gwas_info[2])
ref_effect <- gwas_info[3]
hapmap3 <- files[4]

#Define output name and path
output_name <- paste0("format_", sub("\\..*", "", gwas_info[1]), ".tsv")
output_path <- paste0(gwas_dir, "/", output_name)

if (file.exists(output_path)) {
  #If the file already exist, emit the expected output
  writeLines(output_name, sub("\\.tsv$", ".txt", output_name))
} else {
  #Function to load GWAS according extension
  load_GWAS <- function(data) {
    extension <- sub(".*\\.(.+)\\.(gz|bgz)$|.*\\.(.+)$", "\\1\\3", basename(data))
    switch(extension,
           "tsv" = fread(data, header = TRUE, sep = "\t"),
           "csv" = fread(data, header = TRUE, sep = ","),
           "txt" = fread(data),
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
  
  #Function to change header
  change_headers <- function(data) {
    #Possible names for the header
    variant <- c("^(variant|uniqID|variant_ID|markername|SNPID|ID|SNP|cptid)$")
    beta <- c("^(beta|effect|effects|b|effect_size|es)$")
    se <- c("^(se|standard_error|StandardError|StdErr|sebeta)$")
    pval <- c("^(pval|p_val|p-val||p.val|pvalue|p_value|p-value|p.value|gc_pvalue|p)$")
    LP <- c("^(lp|neg_log_10_p_value|pvalue_mlog|mlogp|log10p)$")
    chr <- c("^(chr|chromosome|chrom|snp-chr|chr_ID|#chrom)$")
    pos <- c("^(pos|base_pair_location|snppos|snp_pos|bp|position|chr_pos|pos\\(b37\\))$")
    z <- c("^(z|zscore|z-score|GC_zscore|tstat|t_stat|t-statistic)$")
    n_samples <- c("^(N|n_complete_samples|sample_size|n_total|TotalSampleSize|SS)$")
    ref <- c("^(ref|reference_allele|ref_allele)$")
    alt <- c("^(alt|alternate_allele|alt_allele)$")
    effect_allele <- c("^(effect_allele|A1|allele1|allele_1|inc_allele|EA)$")
    non_effect_allele <- c("^(other_allele|A2|allele2|allele_2|non_effect_allele|dec_allele|NEA|allele0)$")
    
    #If a column presents rsID, change the column name for rsID
    colnames(data) <- ifelse(apply(head(data, 50), 2, function(col) { 
      any(stri_detect_regex(col, "^rs", case_insensitive = TRUE))
    }), "rsID", colnames(data))
    
    #Change column names to the stablished format
    colnames(data) <- stri_replace_all_regex(colnames(data),
                                             pattern = c(variant, beta, se, z, ref, alt,
                                                         pval, LP, chr, pos, n_samples,
                                                         effect_allele, non_effect_allele),
                                             replacement = c("variant", "beta", "se", "z", "ref", "alt",
                                                             "pval", "LP", "chr", "pos", "n_samples",
                                                             "effect_allele", "non_effect_allele"),
                                             vectorize = FALSE, case_insensitive = TRUE)
    
    #Determine beta as a numeric variable
    data$beta <- as.numeric(data$beta)
    return(data)
  }
  
  #Call the function to change headers
  gwas <- change_headers(gwas)

  #If it is necessary, transform -log10(p-value) to p-value
  if (!("pval" %in% colnames(gwas)) && "LP" %in% colnames(gwas)) {
    gwas$LP <- as.numeric(gwas$LP)
    gwas$pval <- 10^(-gwas$LP)
  }
  
  #If it is necessary, calculate SE
  if (!("se" %in% colnames(gwas))) {
    if ("z" %in% colname(gwas)) {
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

  #If it is necessary, add N
  if (!("n_samples" %in% colnames(gwas))) {
    gwas$n_samples = n_sample
  }
  
  #Function for change the header names in allele columns
  alleles_f <- function(data, ref_effect) {
    if (ref_effect == TRUE) {
      
      data$beta <- -data$beta
      
      if ("effect_allele" %in% colnames(data)) {
        colnames(data)[match("effect_allele", colnames(data))] <- "ref"
        colnames(data)[match("non_effect_allele", colnames(data))] <- "alt"
        }
    
    } else {
      
      if ("effect_allele" %in% colnames(data)) {
        colnames(data)[match("effect_allele", colnames(data))] <- "alt"
        colnames(data)[match("non_effect_allele", colnames(data))] <- "ref"
        }
    }
  
  return(data)
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
    gwas <- alleles_f(gwas, ref_effect)
  }
  
  #Create function to add missing columns to the GWAS from SNPdb file (DELETE !!!!!)
  add_info_from_rsID <- function(gwas, vcf_dir, missing_columns) {
    #Initialize and empty data frame for the results
    SNPdb_info <- data.frame()
    
    #Determine the chromosomes
    if ("chr" %in% colnames(gwas)) {
      chromosome <- unique(gwas$chr)
    } else {
      chromosome <- 1:22
    }
    
    #Initialize a loop to go through all the VCF files
    for (i in chromosome) {
      #Determine file path for every chromosome
      vcf_files = paste0(vcf_dir,"/ALL.chr",i,".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
      #Load SNPdb
      vcf_data <- read.vcfR(vcf_files)
      SNPdb_chr <- getFIX(vcf_data)[, -c(6, 7)]
      #Filter SNPdb to obtain info about the variants present in the GWAS
      SNPdb_filtered_chr <- SNPdb_chr[SNPdb_chr[,"ID"] %in% unique(gwas$rsID),]
      #Create data frame with relevant info
      info_SNP <- data.frame(
        rsID = SNPdb_filtered_chr[, "ID"],
        chr =  SNPdb_filtered_chr[, "CHROM"],
        pos = SNPdb_filtered_chr[, "POS"],
        ref = SNPdb_filtered_chr[, "REF"],
        alt = SNPdb_filtered_chr[, "ALT"],
        stringsAsFactors = FALSE)
      #Filter with missing columns
      info_SNP <- info_SNP %>%
        dplyr::select(c("rsID", all_of(missing_columns)))
      #Add the result from each file to the general data frame
      SNPdb_info <- dplyr::bind_rows(SNPdb_info, info_SNP)
    }
    
    #Add missing columns to the GWAS
    gwas <- gwas %>%
      left_join(SNPdb_info, by = "rsID")
    
    return(gwas)
  }
  
  #Create function to add missing columns from variant column or from the function previously created for rsID
  #It also discard SNPs from chr X, Y, MT and variants that are not 1 substitution base, and create a variant column. 
  variant_function <- function(gwas) {
    #Determine if there is missing columns about chromosome, position or alleles.
    columns_variants <- c("chr", "pos", "ref", "alt")
    missing_columns <- setdiff(columns_variants, colnames(gwas))

    #Replace any separator (-, /, or _) in the variant variable with :
    if("variant" %in% colnames(gwas)) {
      gwas$variant <- gsub("[-/_]", ":", gwas$variant)
    }

    #Add missing columns
    if (length(missing_columns) != 0) {
      #From variant column structured as chr:pos:ref:alt
      if ("variant" %in% colnames(gwas) && all(grepl("^([0-9]+|X|Y|MT|M):[0-9]+:[A-Za-z]+:[A-Za-z]+$", head(gwas$variant, 10)))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 4)
        gwas$chr <- split_variant[, 1]
        gwas$pos <- split_variant[, 2]
        gwas$ref <- split_variant[, 3]
        gwas$alt <- split_variant[, 4]
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
      #From variant column structured as chr:pos
      } else if ("variant" %in% colnames(gwas) && any(grepl("^[0-9]+:[0-9]+$", head(gwas$variant, 10))) && identical(missing_columns, c("chr", "pos"))) {
        split_variant <- str_split_fixed(gwas$variant, ":", 2)
        gwas$chr <- split_variant[, 1]
        gwas$pos <- split_variant[, 2]
        gwas$variant <- paste0(split_variant[, 1], ":", split_variant[, 2])
      #From rsID ######## SE OBTIENEN DEL OTRO SCRIPT ASI QUE HACER LO NECESARIO !!!
      } else if ("rsID" %in% colnames(gwas)) {
        gwas <- add_info_from_rsID(gwas, vcf_dir, missing_columns)
      }
    }
    #Determine position and chromosome as numeric variables
    gwas$pos <- as.numeric(gwas$pos)
    #Discard SNPs from chr X, Y, MT and variants that are not 1 substitution base
    gwas$chr <- as.numeric(gwas$chr)
    gwas <- gwas[!is.na(gwas$chr), ]

    #Determine ref and alt columns as string and select variants with 1 substitution base.
    gwas$ref <- as.character(gwas$ref)
    gwas$alt <- as.character(gwas$alt)
    gwas <- gwas[nchar(gwas$ref) == 1 & nchar(gwas$alt) == 1, ]

    if ("variant" %in% colnames(gwas)) {
      gwas$variant <- ifelse(
        grepl("^(1[0-9]|2[0-2]|[1-9]):[0-9]+:[A-Za-z]+:[A-Za-z]+$", gwas$variant), 
        sub(":([A-Za-z]+:[A-Za-z]+)$", "", gwas$variant), gwas$variant)
    } else {
      #Add column variant if it is necessary
      gwas$variant <- paste0(gwas$chr, ":", gwas$pos)
    }

    #Determine variant as a string
    gwas$variant <- as.character(gwas$variant)

    return(gwas)
  }
  
  #Call the function
  gwas <- variant_function(gwas)
  
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
  columns_req = c("rsID", "variant", "n_samples", "beta", "se", 
                  "pval", "z", "ref", "alt", "chr", "pos")
  missing_columns = setdiff(columns_req, colnames(gwas))
  
  if (length(missing_columns) > 0) {
    message("There is one or more columns missing, revise ", gwas_name, " file to verify it. In case there is no missing column, please adjust the header to match the script requirements.")
    stop("Error: The following columns are missing: ", paste(missing_columns, collapse = ", "))
  } else {
    #Create GWAS output with the established format
    gwas_format <- gwas[, c("variant", "n_samples", "beta", "se", "pval",
                            "z", "ref", "alt", "chr", "pos", "rsID")]


  }

      #gwas_format <- gwas %>%
      #mutate(across(all_of(setdiff(columns_req, names(.))), ~NA)) %>%
      #select(all_of(columns_req))
  
  #Remove unnecessary data
  rm(gwas)
  
  #Create output file
  write_tsv(gwas_format, file = paste0(gwas_dir, "/", output_name))
  
  #Pass the GWAS name (ended in tsv) to the next process in a txt file
  writeLines(output_name, sub("\\.tsv$", ".txt", output_name))
}
