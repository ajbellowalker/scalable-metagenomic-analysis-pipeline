process MULTIQC {

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    tag "multiqc"

    input:
    path results

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc ${results} -o .
    """
}