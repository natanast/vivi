// Workflow definition
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    // Use the module
    include { FASTQC } from './modules/fastqc/main.nf'

    fastqc_ch = FASTQC(read_pairs_ch)
}

workflow.onComplete {
    log.info (workflow.success ? "\nDone! Check the FastQC reports in --> $params.outdir\n" : "Oops .. something went wrong")
}
