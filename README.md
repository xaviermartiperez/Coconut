# Coconut
This repository includes the necessary scripts to run the Coconut pipeline.

coconut.nf --> Nextflow script with the Coconut pipeline. Adjust parameters as needed according to your cluster
prepare_relate_data.nf --> Helper Nextflow pipeline to obtain the necessary per-population -rds file. Adjust parameters as needed according to your cluster
nextflow.config --> Configuration file needed by Nextflow to run. Adjust parameters as needed
nextflow.slurm --> Shell script to launch Either coconut.nf or prepare_relate_data.nf to you cluster. This is considered a the best practice when launching Nextflow jobs to a cluster
/bin --> Directory from which Nextflow automatically draws the necessary scripts for each process. The scripts inside should not require modifications.
