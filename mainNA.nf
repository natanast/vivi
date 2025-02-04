nextflow.enable.dsl=2

include { TRUST4 } from './modules/nf-core/trust4/main.nf'

// Create channels from samplesheet
workflow {
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.bam))
        }
        .set { bam_ch }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.fasta))
        }
        .set { fasta_ch }

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.sample]
            tuple(meta, file(row.vdj_reference))
        }
        .set { vdj_reference_ch }

    TRUST4(bam_ch, fasta_ch, vdj_reference_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Check the TRUST4 reports in --> $params.outdir\n" : "Oops .. something went wrong" )
}


// // #!/usr/bin/env nextflow

// // log.info """\
// //     F A S T Q C   P I P E L I N E
// //     ===============================
// //     reads        : ${params.samplesheet}
// //     outdir       : ${params.outdir}
// //     """
// //     .stripIndent()

// // // Include modules
// // include { FASTQC } from './modules/nf-core/fastqc/main.nf'

// // /*
// //  * Read the sample sheet and create a channel
// //  */
// // Channel
// //     .fromPath(params.samplesheet)
// //     .splitCsv(header: true)
// //     .map { row -> tuple(row.sample, [row.read1, row.read2]) }
// //     .set { read_pairs_ch }

// // /*
// //  * Workflow definition
// //  */
// // workflow {
// //     /*
// //      * Run FastQC on the read pairs
// //      */
// //     fastqc_ch = FASTQC(read_pairs_ch)
// // }

// // workflow.onComplete {
// //     log.info ( workflow.success ? "\nDone! Check the FastQC reports in --> $params.outdir\n" : "Oops .. something went wrong" )
// // }
