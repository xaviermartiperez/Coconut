#!/bin/bash
#SBATCH --job-name=NF-(J)PALM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -o nextflow.out
#SBATCH --mail-user=xavier.martii@upf.edu
#SBATCH --mail-type=ALL

#Load Nextflow
module load Nextflow/24.10.2

#Allocate some reasonable memory in the cluster
#export NXF_OPTS="-Xms100G -Xmx120G"

#Run the Nextflow script
nextflow run coconut.nf -resume
