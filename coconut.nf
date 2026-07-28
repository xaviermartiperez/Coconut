#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define main parameters
params.user = "xmarti"
params.analysis = "PALM" // Options: PALM, JPALM
params.GWAS_list = "gwas_of_gcat_list.txt"
params.GWAS_dir = "gwas_of_gcat"
params.population = "IBS"
params.maxp = 5e-6
// The assembly version refers to the one of the Relate data, the one
// that the GWAS must be in the end
params.version = "hg19" // Options: hg19, hg38

// Define directories
params.workDir = "/homes/users/${params.user}/scratch/gcat"
params.GWAS_raw_dir = "${params.workDir}/${params.GWAS_dir}"
params.VCFDir = "/homes/users/${params.user}/scratch/dbSNP"
params.LiftOverDir = "${params.workDir}/liftover"
params.LDSCDir = "${params.workDir}/ldscore"
interDir = "${params.workDir}/intermediate_files_${params.GWAS_dir}_${params.population}_${params.version}"
likDir = "${params.workDir}/SNP_likelihoods_${params.population}_${params.version}"
palmDir = "${params.workDir}/results_marginal_palm_${params.population}_${params.version}"
jpalmDir = "${params.workDir}/results_joint_palm_${params.population}_${params.version}"

params.PackagesDir = "/homes/users/${params.user}/scratch/github_packages"
LDSC_ToolDir = "${params.PackagesDir}/ldsc"
Relate_ToolDir = "${params.PackagesDir}/relate/scripts/SampleBranchLengths"
PALM_ToolDir = "${params.PackagesDir}/palm"

// Define where the GWAS list file is
params.GWAS_list_file = "${params.workDir}/${params.GWAS_list}"

// Define LDSC-related tools, files and directories
// IT IS CRITICAL THAT LDSC DIR IS IN BRANCH 2-TO-3
params.LDSC = "${params.LDSCDir}/eur_w_ld_chr" // WHY NO hg19 - hg38 VERSION SPECIFICATION??
params.MungeTool = "${LDSC_ToolDir}/munge_sumstats.py"
params.LDSC_Tool = "${LDSC_ToolDir}/ldsc.py"
// Hapmap3 SNP file must be defined in this manner because each GWAS may require a different version
params.Hapmap3_hg38 = "${params.LDSCDir}/hm3_SNPs_variant_id_hg38.tsv"
params.Hapmap3_hg19 = "${params.LDSCDir}/hm3_SNPs_variant_id_hg19.tsv"
// The same happens with the LD blocks file
params.LD_blocks_hg38 = "${params.workDir}/ld_blocks_hg38.bed" 
params.LD_blocks_hg19 = "${params.workDir}/ld_blocks_hg19.bed"

// Define LiftOver tool and chain files
params.LiftOverTool = "${params.LiftOverDir}/liftOver"
params.LiftOverChain = (params.version == "hg38") ? "${params.LiftOverDir}/hg19ToHg38.over.chain" : 
                       (params.version == "hg19") ? "${params.LiftOverDir}/hg38ToHg19.over.chain" : null

// Define RELATE-related tools, files and directories
Relate_SNPs = (params.version == "hg19") ? "${params.workDir}/relate_snps_${params.population}.rds" :
              (params.version == "hg38") ? "${params.workDir}/relate_snps_${params.population}_hg38.rds" : null
Relate_coal = (params.version == "hg19") ? "${params.workDir}/relate/${params.population}/coal/${params.population}.coal" : 
              (params.version == "hg38") ? "${params.workDir}/relate/${params.population}_hg38/coal/${params.population}_hg38.coal" : null
Relate_anc_mut = (params.version == "hg19") ? "${params.workDir}/relate/${params.population}/anc_mut" : 
                 (params.version == "hg38") ? "${params.workDir}/relate/${params.population}_hg38/anc_mut" : null
params.SBL_Tool = "${Relate_ToolDir}/SampleBranchLengths_custom.sh"

// Define PALM-related tools and directories
params.Lik_Tool = "${PALM_ToolDir}/lik.py"
params.PALM_Tool = "${PALM_ToolDir}/palm_custom.py"

// Other numerical parameters
params.N_per_batch = 10
params.min_SI_SNPs = 25

log.info """\
    P A L M  /  J - P A L M   W O R K F L O W
    =========================================
    Mode        : ${params.analysis}
    Population  : ${params.population}
    Assembly    : ${params.version}
    GWAS list   : ${params.GWAS_list}
    Relate      : ${Relate_SNPs}

    Max. p-val  : ${params.maxp}
    """
    .stripIndent()

process Unlist_GWAS {

    clusterOptions = '--partition=haswell'

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 8

    input:
    path GWAS_list

    output:
    path 'gwas_*'

    script:
    """
    unlist_GWAS.sh ${GWAS_list}
    """

}

process Format_GWAS {

    module = 'R/4.3.2-gfbf-2023a'
    //clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 8

    input:
    path GWAS
    path GWASdir
    path Hapmap3_hg38
    path Hapmap3_hg19

    when:

    output:
    path 'formatted_*.txt'

    script:
    """
    format_gwas.R ${GWAS} ${GWASdir} ${Hapmap3_hg38} ${Hapmap3_hg19}
    """

}
// This process is not finished, beware if sumstats do not contain all
// the necessary columns!!!
process Fix_Sumstats {

    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 8

    input:
    path GWASinfo
    path GWASdir
    path VCFdir

    output:
    path 'fixed_*.txt'

    script:
    """
    fix_sumstats.sh ${GWASinfo} ${GWASdir} ${VCFdir}
    """

}

process Add_Filter_LD_blocks {

    module = 'R/4.3.2-gfbf-2023a'
    clusterOptions = '--partition=haswell'
    
    publishDir interDir, mode: 'copy'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path GWAS
    path GWASdir
    path LD_blocks_hg38
    path LD_blocks_hg19
    val MaxPval

    output:
    path '*-+*'

    script:
    """
    ld_block_filter.R ${GWAS} ${GWASdir} ${LD_blocks_hg38} ${LD_blocks_hg19} ${MaxPval}
    """

}

process LiftOver {

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 5

    input:
    path GWAS
    path GWAS_dir
    path LiftOverTool
    path LiftOverChain
    val version

    output:
    path '*_liftOver_checked.txt.txt'

    when: 
    GWAS.exists() &&
        GWAS.name.tokenize('-+').size() == 2 &&
        Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs

    script:
    """
    gwas_liftover.sh ${GWAS} ${GWAS_dir} ${LiftOverTool} ${LiftOverChain} ${version}
    """
}

process Obtain_N_Comparisons {

    clusterOptions = '--partition=haswell'

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path GWAS_list

    output:
    path 'n_pairwise_comparisons.txt'

    when:
    params.analysis == "JPALM"
 
    script:
    """
    n_comparisons.sh ${GWAS_list}
    """

}

process Obtain_GWAS_Pairs {

    clusterOptions = '--partition=haswell'

    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path GWAS_list

    output:
    path 'comparison_batch_*'

    when:
    params.analysis == "JPALM"
 
    script:
    """
    gwas_pairs_batch_obtention.sh ${GWAS_list}
    """

}

process Munge_Sumstats {
    
    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path GWAS
    path GWASdir
    val projectDir
    path ldscore

    output:
    path '*.sumstats.info'

    when:
    params.analysis == "JPALM" &&
    GWAS.exists() &&
        GWAS.name.tokenize('-+').size() == 2 &&
        Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs

    script:
    """
    munge_sumstats.sh ${GWAS} ${GWASdir} ${projectDir} ${ldscore}
    """

}

process Genetic_Correlation {

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path GWAS_Pair
    path N_Comparisons
    path projectDir
    path ldscore

    output:
    path '*-+*-+*-+*'

    when:
    params.analysis == "JPALM"

    script:
    """
    genetic_correlation_batch.sh ${GWAS_Pair} ${N_Comparisons} ${projectDir} ${ldscore}
    """

}

process SelectSNPs {

    module = 'R/4.3.2-gfbf-2023a'
    clusterOptions = '--partition=haswell'

    publishDir interDir, mode: 'copy'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path relate_SNPs
    path GWAS
    path GWAS_dir

    output:
    path '*_selected_SNPs.tsv'

    when:
    params.analysis == "PALM" //&& 
   //(
        //(params.version in ["hg19", "hg38"] && 
        //GWAS.exists() && 
        //GWAS.name.tokenize('-+').size() == 2 && 
        //Float.parseFloat(GWAS.name.tokenize('-+')[1]) >= params.min_SI_SNPs) 
        //|| 
        //params.version in ["hg19tohg38", "hg38tohg39"]
    //)

    script:
    """
    snp_selection.R ${relate_SNPs} ${GWAS} ${GWAS_dir}
    """

}

process SelectSNPs_JPALM {

    module = 'R/4.3.2-gfbf-2023a'
    clusterOptions = '--partition=haswell'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    path rg
    path relate_SNPs
    path GWAS_dir
    val MaxPval

    output:
    path '*-*_selected_SNPs.tsv'

    when:
    params.analysis == "JPALM" &&
    params.version in ["hg19", "hg38"] && 
    rg.exists()  &&
        rg.name.tokenize('-+').size() == 4 &&
        Math.abs(Float.parseFloat(rg.name.tokenize('-+')[2])) > 0.2 &&
        Float.parseFloat(rg.name.tokenize('-+')[3]) < 0.005

    script:
    """
    snp_selection_jpalm.R ${rg} ${relate_SNPs} ${GWAS_dir} ${MaxPval}
    """

}

process SelectSNPs_JPALM_lifted {

    module = 'R/4.3.2-gfbf-2023a'
    clusterOptions = '--partition=haswell'

    publishDir interDir, mode: 'copy'

    memory { 16.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path rg
    path relate_SNPs
    path GWAS_dir
    val MaxPval
    val version

    output:
    path '*-*_selected_SNPs.tsv'

    when:
    params.analysis == "JPALM" &&
    //params.version in ["hg19tohg38", "hg38tohg19"] &&
    rg.exists()  &&
        rg.name.tokenize('-+').size() == 4 &&
        Math.abs(Float.parseFloat(rg.name.tokenize('-+')[2])) > 0.2 &&
        Float.parseFloat(rg.name.tokenize('-+')[3]) < 0.005

    script:
    """
    snp_selection_jpalm_lift.R ${rg} ${relate_SNPs} ${GWAS_dir} ${MaxPval} ${version}
    """
}

process MergeSNPs {

    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path Selected_SNPs_per_GWAS

    output:
    path 'merged_SNPs.tsv'

    script:
    """
    merge_SNPs.sh ${Selected_SNPs_per_GWAS}
    """

}

process TrimSNPs {

    clusterOptions = '--partition=haswell'

    memory { 2.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..144 ? 'retry' : 'terminate' }
    maxRetries 1

    input:
    path Selected_SNPs
    val MaxPval
    // Analysis must here be specified within the script as it directly
    // influences the code
    val analysis
    val N_per_batch

    output:
    path 'batch_*_*.tsv'

    script:
    """
    trim_snps_batch.sh ${Selected_SNPs} ${MaxPval} ${analysis} ${N_per_batch}
    """

}

process SBL_and_SNP_Likelihood {

    clusterOptions = '--partition=haswell'

    time '1h'
    memory 8.GB

    input:
    path SNP_batch
    val population
    path SBL_tool
    path lik_tool
    path Relate_coal
    path Relate_anc_mut
    path LikDir

    output:
    path 'bp*.quad_fit.npy'

    script:
    """
    sbl_and_snp_likelihood.sh ${SNP_batch} ${population} ${SBL_tool} ${lik_tool} ${Relate_coal} ${Relate_anc_mut} ${LikDir}
    """

}

process Apply_PALM {

    publishDir palmDir, mode: 'copy'

    clusterOptions = '--partition=haswell'

    memory { 4.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path Filtered_GWAS
    val MaxPvalue
    path PALM_Tool
    path likDir

    output:
    path '*_marginal_PALM.txt'

    when:
    params.analysis == "PALM"

    script:
    """
    palm.sh ${Filtered_GWAS} ${MaxPvalue} ${PALM_Tool} ${likDir}
    """

}

process Apply_JPALM {
    
    publishDir jpalmDir, mode: 'copy'

    clusterOptions = '--partition=haswell'

    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus in 135..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    val ready
    path Filtered_GWAS_Pair
    val MaxPvalue
    path PALM_Tool
    path likDir

    output:
    path '*_J_PALM.txt'
    path '*_significant_independent_SNPs.txt'

    when:
    params.analysis == "JPALM"

    script:
    """
    j_palm.sh ${Filtered_GWAS_Pair} ${MaxPvalue} ${PALM_Tool} ${likDir}
    """

}

workflow {
    Unlist_GWAS_ch = Unlist_GWAS(params.GWAS_list_file)
    Format_GWAS_ch = Format_GWAS(Unlist_GWAS_ch.flatten(), params.GWAS_raw_dir, params.Hapmap3_hg38, params.Hapmap3_hg19)
    Fix_Sumstats_ch = Fix_Sumstats(Format_GWAS_ch, params.GWAS_raw_dir, params.VCFDir)
    Add_Filter_LD_blocks_ch = Add_Filter_LD_blocks(Fix_Sumstats_ch, params.GWAS_raw_dir, params.LD_blocks_hg38, params.LD_blocks_hg19, params.maxp)
    
    Lifted_GWAS_ch = LiftOver(Add_Filter_LD_blocks_ch, params.GWAS_raw_dir, params.LiftOverTool, params.LiftOverChain, params.version)
    
    //if (params.analysis == "JPALM") {
     //   GWAS_N_ch = Obtain_N_Comparisons(Add_Filter_LD_blocks_ch.collect())
      //  GWAS_Pairs_ch = Obtain_GWAS_Pairs(Add_Filter_LD_blocks_ch.collect())
       // Sumstats_ch = Munge_Sumstats(Add_Filter_LD_blocks_ch, params.GWAS_raw_dir, workflow.projectDir, params.LDSC)
        //RG_ch = Genetic_Correlation(Sumstats_ch.collect(), GWAS_Pairs_ch.flatten(), GWAS_N_ch, workflow.projectDir, params.LDSC)
        //SNPs_ch = SelectSNPs_JPALM_lifted(Lifted_GWAS_ch.collect(), RG_ch.flatten(), Relate_SNPs, params.GWAS_raw_dir, params.maxp, params.version)
    //} else if (params.analysis == "PALM") {
    //    SNPs_ch = SelectSNPs(Relate_SNPs, Lifted_GWAS_ch, params.GWAS_raw_dir)
    //} else {
    //    log.info "Analysis parameter must be either PALM or JPALM"
    //    exit 1
    //}

    //Merged_SNPs_ch = MergeSNPs(SNPs_ch.collect())
    //Trimmed_SNPs_ch = TrimSNPs(Merged_SNPs_ch, params.maxp, params.analysis, params.N_per_batch)
    //SBL_Lik_ch = SBL_and_SNP_Likelihood(Trimmed_SNPs_ch.flatten(), params.population, params.SBL_Tool, params.Lik_Tool, Relate_coal, Relate_anc_mut, likDir)
    
    //if (params.analysis == "PALM") {
    //    PALM_ch = Apply_PALM(SBL_Lik_ch.collect(), SNPs_ch, params.maxp, params.PALM_Tool, likDir)
    //} else if (params.analysis == "JPALM") {
    //    J_PALM_ch = Apply_JPALM(SBL_Lik_ch.collect(), SNPs_ch, params.maxp, params.PALM_Tool, likDir)
    //} else {
    //    log.info "Analysis parameter must be either PALM or JPALM"
    //    exit 1
    //}
    
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "\nOops... something went wrong\n" )
}