#!/usr/bin/env nextflow

include { paramsSummaryLog        } from 'plugin/nf-schema'
include { validateParameters      } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'
include { GRZQC                   } from './workflows/grzqc'

workflow {

    main:
    // Display parameters that differ from the schema default
    log.info paramsSummaryLog(workflow)

    // Validate the parameters using nf-schema
    validateParameters()

    if( params.submission && params.samplesheet ) {
        error "Please provide either --submission OR --samplesheet, not both."
    }

    if( !params.submission && !params.samplesheet ) {
        error "You must provide either --submission or --samplesheet."
    }

    if( params.submission ) {
        // create samplesheet from submission metadata.json file
        METADATA_TO_SAMPLESHEET(
            Channel.fromPath(params.submission, checkIfExists: true)
        )
        samplesheets = METADATA_TO_SAMPLESHEET.out.samplesheet
    } else if ( params.samplesheet ) {
        // Use the provided samplesheet file directly
        samplesheets = Channel.fromPath(params.samplesheet, checkIfExists: true)
    }

    // construct the samples channel using the schema
    samples = samplesheets
        .flatMap { samplesheet ->
            samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json")
        }
        .map { meta, reads1, reads2 ->
            if ( !reads2 ) {
                [ meta + [ single_end: true ], [ reads1 ] ]
            } else {
                [ meta + [ single_end: false ], [ reads1, reads2 ] ]
            }
        }
        .groupTuple()
        .map {
            meta, reads -> [ meta, reads.flatten() ]
        }
        .map{ meta, reads ->
          meta.targets = meta.targets ? meta.targets : file(params.genomes[meta.reference].targets)
          [ meta, reads ]
        }
        .dump(tag: 'samples')

    GRZQC(
      samples,
      params.genomes,
      params.reference_path,
      params.thresholds
    )

}
