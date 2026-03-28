process MEGAHIT {

    publishDir "${params.outdir}/megahit", mode: 'copy'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}_assembly")

    script:
    """
    megahit \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${sample_id}_assembly \
    -t 4 \
    -m 0.5
    """
}
