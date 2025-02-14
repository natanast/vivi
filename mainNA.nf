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