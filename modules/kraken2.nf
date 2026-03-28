process KRAKEN2 {

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}.report")

    script:
    """
    kraken2 \
    --db ${params.kraken_db} \
    --paired ${reads[0]} ${reads[1]} \
    --threads ${params.threads} \
    --report ${sample_id}.report
    """
}