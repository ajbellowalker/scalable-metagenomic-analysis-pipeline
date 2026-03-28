nextflow.enable.dsl=2

include { KRAKEN2 } from '../modules/kraken2'
include { MEGAHIT } from '../modules/megahit'
include { METABAT2 } from '../modules/metabat2'
include { CHECKM2 } from '../modules/checkm2'
include { GTDBTK } from '../modules/gtdbtk'
include { MULTIQC } from '../modules/multiqc'

workflow METAGENOME_PIPELINE {

    reads_ch = Channel.fromFilePairs(params.reads)

    kraken_out = KRAKEN2(reads_ch)

    if (params.run_megahit) {
    assembly_out = MEGAHIT(reads_ch)}

    if (params.run_metabat2) {
        bins = METABAT2(assembly_out)

        if (params.run_checkm2) {
            qc = CHECKM2(bins)
        }

        if (params.run_gtdbtk) {
            taxonomy = GTDBTK(bins)
        }
    }

    // Collect all outputs for MultiQC
    results_ch = Channel.fromPath("${params.outdir}")

    MULTIQC(results_ch)
}
