process TRUST4 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::trust4=1.0.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/trust4:1.1.5--h5ca1c30_0' :
        'biocontainers/trust4:1.1.5--h5ca1c30_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(vdj_reference)

    output:
    tuple val(meta), path("*.tsv")                  , emit: tsv
    tuple val(meta), path("*_airr.tsv")             , emit: airr_files
    tuple val(meta), path("${meta.id}_airr.tsv")    , emit: airr_tsv
    tuple val(meta), path("*_report.tsv")           , emit: report_tsv
    tuple val(meta), path("*.fa")                   , emit: fasta
    tuple val(meta), path("*.out")                  , emit: out
    tuple val(meta), path("*.fq")                   , emit: fq
    tuple val(meta), path("**")                     , emit: outs
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_mode = bam ? "-b ${bam}" : ''
    def reference = vdj_reference ? "--ref ${vdj_reference}" : ""
    

    def resultsDir = "${projectDir}/results/${prefix}/" 

    """
    mkdir -p ${resultsDir}

    run-trust4 \\
        ${bam_mode} \\
        -t $task.cpus \\
        -f ${fasta} \\
        -o ${prefix} \\
        ${reference} \\
        $args

    # Copy output files to the desired directory
    cp -r ${prefix}_* ${resultsDir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^/TRUST4 using /' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_airr.tsv
    touch ${prefix}_airr_align.tsv
    touch ${prefix}_report.tsv
    touch ${prefix}_assembled_reads.fa
    touch ${prefix}_annot.fa
    touch ${prefix}_cdr3.out
    touch ${prefix}_raw.out
    touch ${prefix}_final.out
    touch ${prefix}_toassemble.fq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trust4: \$(run-trust4 2>&1 | grep -o 'v[0-9.]*-r[0-9]*' | sed 's/^/TRUST4 using /' )
    END_VERSIONS
    """
}
