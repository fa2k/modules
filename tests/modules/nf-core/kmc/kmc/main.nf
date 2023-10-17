#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { KMC_KMC } from '../../../../../modules/nf-core/kmc/kmc/main.nf'

workflow test_kmc_kmc {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
       
    ]

    KMC_KMC ( input )
}
