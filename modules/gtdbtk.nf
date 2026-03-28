process GTDBTK {

    publishDir "${params.outdir}/gtdbtk", mode: 'copy'

    tag "gtdbtk"

    input:
    path bins

    output:
    path "gtdbtk_out"

    script:
    """
    gtdbtk classify_wf \
    --genome_dir ${bins} \
    --out_dir gtdbtk_out \
    --cpus ${params.threads}
    """
}