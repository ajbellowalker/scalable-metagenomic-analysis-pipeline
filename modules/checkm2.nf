process CHECKM2 {

    publishDir "${params.outdir}/checkm2", mode: 'copy'

    tag "checkm2"

    input:
    path bins

    output:
    path "checkm2_out"

    script:
    """
    checkm2 predict \
    -i ${bins} \
    -x fa \
    -o checkm2_out \
    --threads ${params.threads}
    """
}