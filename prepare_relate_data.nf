#!/usr/bin/env nextflow

//This pipeline comes with a "nextflow.config" file which specifies
//that all processes are to be executed by SLURM. 

//Define parameters
params.WorkDir = "$projectDir"
params.RelateDir = "$projectDir/relate"
params.VCFDir= "/gpfs/projects/lab_dcomas/1000genomes_phase3_dcomas/vcf"
params.SampleIDs = "$projectDir/igsr-1000_genomes_phase_3_release.tsv"
params.LD_Blocks = "$projectDir/ld_blocks.bed"
params.coal = "$projectDir/relate/ALL_RESULTS/1000GP_Phase3_mask_prene.coal"
params.interdir = "$projectDir/intermediate_files_relate"
params.outdir = "$projectDir/results_relate_pipeline"

Chromosomes_ch = Channel.of(1..22)
Populations_ch = Channel.of('FIN','CEU','GBR','TSI','IBS')

POP_CHR_ch = Populations_ch.combine(Chromosomes_ch)



log.info """\
    R E L A T E   F I L E S - N F   P I P E L I N E
    ===============================================
    Relate Dir	: ${params.RelateDir}
    VCF Dir     : ${params.VCFDir}
    Sample IDs  : ${params.SampleIDs}
    Interdir    : ${params.interdir}
    Outdir      : ${params.outdir}
    """
    .stripIndent()


process Annotate_Mut {

    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path Mut_Dir
    val(Chromosome)

    output:
    path '*_mut_report.txt'


    script:
    """
    add_mut_annotation_relate.sh ${Mut_Dir} ${Chromosome}
    """

}

process Extract_Anc_Mut {

    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    //val ready
    tuple val(Population), val(Chromosome)
    path Anc_Mut_Dir

    output:
    path '*_extract_report.txt'

    script:
    """
    extract_relate_populations.sh ${Population} ${Chromosome} ${Anc_Mut_Dir}
    """

}

process Extract_VCF_and_AF {
    
    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path IDs
    path VCFDir
    path RelateDir
    tuple val(Population), val(Chromosome)

    output:
    path '*_report.txt'

    script:
    """
    extract_1000G_VCF_and_AF.sh ${IDs} ${VCFDir} ${RelateDir} ${Population} ${Chromosome}
    """

}

process Extract_coal {
    
    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path coal
    val population
    path RelateDir

    output:
    path '*_coal.report'

    script:
    """
    extract_coal_files.sh ${coal} ${population} ${RelateDir} 
    """

}

process Prepare_Relate_RDS {

    module = 'R/4.2.0-foss-2021b'
    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir params.interdir, mode:'copy' 

    input:
    val ready
    val ready
    val Populuation
    path RelateDir
    path WorkDir
    path LD

    output:
    path 'relate_snps_*_report.txt'

    script:
    """
    relate_snps_preparation.R ${Populuation} ${RelateDir} ${WorkDir} ${LD}
    """
  
}

workflow {
    Mut_ch = Annotate_Mut(params.RelateDir, Chromosomes_ch)
    Anc_Mut_ch = Extract_Anc_Mut(Mut_ch.collect(), Chromosomes_ch.flatten(), Populations_ch.flatten(), paramms.RelateDir)
    Anc_Mut_ch = Extract_Anc_Mut(POP_CHR_ch, params.RelateDir)
    VCF_AF_ch = Extract_VCF_and_AF(params.SampleIDs, params.VCFDir, params.RelateDir, POP_CHR_ch)
    Coal_ch = Extract_coal(params.coal, Populations_ch.flatten(), params.RelateDir)
    Relate_RDS_ch = Prepare_Relate_RDS(Anc_Mut_ch.collect(), VCF_AF_ch.collect(), Populations_ch.flatten(), params.RelateDir, params.WorkDir, params.LD_Blocks)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! All the necessary files should be ready now" : "Oops... something went wrong" )
}
