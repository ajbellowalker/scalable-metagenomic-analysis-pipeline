nextflow.enable.dsl=2

include { KRAKEN2 } from '../modules/kraken2'
include { MEGAHIT } from '../modules/megahit'

workflow METAGENOME_PIPELINE {

    reads_ch = Channel.fromFilePairs(params.reads)

    kraken_out = KRAKEN2(reads_ch)

    assembly = MEGAHIT(reads_ch)
}