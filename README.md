This repository includes the necessary scripts to run the Coconut and helper pipelines.

Necessary software and dependancies:

- To install Nextflow, visit https://docs.seqera.io/nextflow/install and follow the instructions.

- To install LDSC, visit https://github.com/bulik/LDSC and clone the repository. 
IMPORTANTLY, clone the pull request 2-to-3: https://github.com/bulik/ldsc/pull/360 and set the Git Head on branch 2-to-3, like this:
git clone https://github.com/bulik/ldsc.git
cd ldsc
git remote add belowlab https://github.com/belowlab/ldsc.git
git fetch belowlab 2-to-3
git checkout -b 2-to-3 belowlab/2-to-3

- To install Relate, visit https://myersgroup.github.io/relate/index.html and follow the instructions.

- To install PALM, visit https://github.com/standard-aaron/palm and clone the repository

- To install R, visit https://www.r-project.org and follow the instructions.
The required packages: data.table, dplyr, stringi, stringr, tidyverse can all be installed through the R function "install.packages()".

- It is advisable that you have a dbSNP/ directory with the latest VCF split into chromosomes in case some sumstats need fixing. The latest version of the dbSNP VCF and index file (.tbi) can be downloaded here: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/

# Coconut
- coconut.nf --> Nextflow script with the Coconut pipeline. Adjust parameters as needed according to your cluster. 
- prepare_relate_data.nf --> Helper Nextflow pipeline to obtain the necessary per-population -rds file. Adjust parameters as needed according to your cluster
- nextflow.config --> Configuration file needed by Nextflow to run. Adjust parameters as needed
- nextflow.slurm --> Shell script to launch Either coconut.nf or prepare_relate_data.nf to you cluster. This is considered a the best practice when launching Nextflow jobs to a cluster
- /bin --> Directory from which Nextflow automatically draws the necessary scripts for each process. The scripts palm_custom.py and SampleBranchLengths_custom.sh should be placed inside the corresponding directories.
- /useful_files --> Directory that contains the Liftover executable (make sure it can be exectured) and the chain files from hg19 to hg38 and vice versa, and Hapmap3 SNP set for LDSC analysis from in hg19 and hg38 coordinates and a compressed file with per-chromosome LD-scores (that need to be decompressed). Make sure to put every file to the corresponding path.
