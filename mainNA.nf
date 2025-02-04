#!/usr/bin/env nextflow

log.info """\
    F A S T Q C   P I P E L I N E
    ===============================
    reads        : ${params.samplesheet}
    outdir       : ${params.outdir}
    """
    .stripIndent()

// Include modules
include { FASTQC } from './modules/nf-core/fastqc/main.nf'

/*
 * Read the sample sheet and create a channel
 */
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample, [row.read1, row.read2]) }
    .set { read_pairs_ch }

/*
 * Workflow definition
 */
workflow {
    // Channel
    //     .fromFilePairs(params.reads, checkIfExists: true)
    //     .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Check the FastQC reports in --> $params.outdir\n" : "Oops .. something went wrong" )
}

// #!/usr/bin/env nextflow

// /*
//  * pipeline input parameters
//  */
// params.reads = "/mnt/d/Users/nanastasiadou/Bioinfo/Projects/vivi/E23221_S39_L001_R{1,2}_001.fastq.gz"
// params.outdir = "results"

// log.info """\
//     F A S T Q C   P I P E L I N E
//     ===============================
//     reads        : ${params.reads}
//     outdir       : ${params.outdir}
//     """
//     .stripIndent()

// /*
//  * define the `FASTQC` process that runs FastQC on the provided reads
//  */
// process FASTQC {
//     tag "FASTQC on $sample_id"

//     input:
//     tuple val(sample_id), path(reads)

//     output:
//     path "fastqc_${sample_id}_logs"

//     script:
//     """
//     echo "Processing $sample_id with FastQC"
//     mkdir -p fastqc_${sample_id}_logs
//     fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
//     """
// }

// /*
//  * Workflow definition
//  */
// workflow {
//     Channel
//         .fromFilePairs(params.reads, checkIfExists: true)
//         .set { read_pairs_ch }

//     fastqc_ch = FASTQC(read_pairs_ch)
// }

// workflow.onComplete {
//     log.info ( workflow.success ? "\nDone! Check the FastQC reports in --> $params.outdir\n" : "Oops .. something went wrong" )
// }
