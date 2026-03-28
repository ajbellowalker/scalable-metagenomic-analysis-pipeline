nextflow.enable.dsl=2

include { METAGENOME_PIPELINE } from './workflows/metagenome'

workflow {
    METAGENOME_PIPELINE()
}