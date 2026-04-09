process METABAT2 {

    publishDir "${params.outdir}/metabat2", mode: 'copy'

    tag "$sample_id"

    input:
    val sample_id
    path assembly_dir

    output:
    path "bins"

    script:
    """
    mkdir bins

    metabat2 \
    -i ${assembly_dir}/final.contigs.fa \
    -o bins/bin

    """
}