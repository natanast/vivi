/*
 * define the `FASTQC` process that runs FastQC on the provided reads
 */
process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    echo "Processing $sample_id with FastQC"
    mkdir -p fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}
