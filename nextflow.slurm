#!/bin/bash
#SBATCH --job-name=NF-(J)PALM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -o nextflow.out
#SBATCH --mail-user=your.email@here.com
#SBATCH --mail-type=ALL

#Load Nextflow
module load Nextflow/24.10.2

#Run the Nextflow script
nextflow run coconut.nf -resume
