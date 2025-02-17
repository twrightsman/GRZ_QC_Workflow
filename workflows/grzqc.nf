/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                 } from '../modules/nf-core/fastp/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { FASTQ_ALIGN_BWA        } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { MOSDEPTH as MOSDEPTH_TARGET   } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grzqc_pipeline'

include { COMPARE_THRESHOLD      } from '../modules/local/compare_threshold'
include { MERGE_REPORTS          } from '../modules/local/merge_reports'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GRZQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Creating merge fastq channel from samplesheet
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, fastqs)
    }.branch {
            meta, fastqs ->
                single_sample_paired_end  : fastqs.size() == 2
                    return [ meta, fastqs.flatten() ]
                multiple_sample_paired_end: fastqs.size() > 2
                    return [ meta, fastqs.flatten() ]
        }.set{
            ch_fastqs
        }

    // Creating bed_file channel from samplesheet : Only taking the first bed file in case there are multiple samples from the same patient
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, bed_file.first())
    }set{
            ch_bed_file
        }


    // Create BWA_index channels for BWA_MEM
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, reference.first())
    }.branch {
            meta, reference ->
                ref_37  : (reference == "GRCh37" || reference == "hg19")
                    return [ meta, params.bwa_index_37 ]
                ref_38  : (reference == "GRCh38" || reference == "hg38")
                    return [ meta, params.bwa_index_38 ]
        }.set{
        ch_reference_ind
        }

    ch_reference_ind.ref_37.mix(
        ch_reference_ind.ref_38
    ).set{
        ch_bwa_index
    }

    // Create bed channels for Mosdepth on ~ 400 representative genes
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, reference.first())
    }.branch {
            meta, reference ->
                rep_genes_hg19  : (reference == "GRCh37" || reference == "hg19")
                    return [ meta, params.rep_genes_hg19 ]
                rep_genes_hg38  : (reference == "GRCh38" || reference == "hg38")
                    return [ meta, params.rep_genes_hg38 ]
        }.set{
        ch_rep_genes_by_ref
        }

    ch_rep_genes_by_ref.rep_genes_hg19.mix(
        ch_rep_genes_by_ref.rep_genes_hg38
    ).set{
        ch_rep_genes
    }

    // Create fasta channels for BWA_MEM
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, reference.first())
    }.branch {
            meta, reference ->
                ref_37  : (reference == "GRCh37" || reference == "hg19")
                    return [ meta, params.fasta_37 ]
                ref_38  : (reference == "GRCh38" || reference == "hg38")
                    return [ meta, params.fasta_38 ]
        }.set{
        ch_fasta_ind
        }


    ch_fasta_ind.ref_37.mix(
        ch_fasta_ind.ref_38
    ).set{
        ch_fasta
    }


    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastqs.multiple_sample_paired_end
    )
    .reads
    .mix(ch_fastqs.single_sample_paired_end)
    .set { ch_cat_fastq }

    // 
    // MODULE: FASTP
    //
    save_trimmed_fail = false
    save_merged = false
    FASTP(
        ch_cat_fastq,
        [], // we are not using any adapter fastas at the moment
        false, // we don't use discard_trimmed_pass at the moment
        save_trimmed_fail,
        save_merged
    )
    ch_fastp_json = FASTP.out.json
        .map { meta, json -> tuple(meta, json) }

    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{ meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{ meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())


    FASTQ_ALIGN_BWA (
        ch_cat_fastq, ch_bwa_index, true, ch_fasta
    )
    ch_bam = FASTQ_ALIGN_BWA.out.bam
    ch_bai = FASTQ_ALIGN_BWA.out.bai

    // Is there a better way to do this, i.e join more than 2 channels into one in one step?
    ch_mosdepth_input_1 = ch_bam.join(ch_bai, by: 0).map { meta, file1, file2 -> tuple(meta, file1, file2) }
    ch_mosdepth_input = ch_mosdepth_input_1.join(ch_bed_file, by: 0).map { meta, file1, file2, bed_file -> tuple(meta, file1, file2, bed_file) }
    //
    // MODULE: MOSDEPTH, get average coverage
    //
    MOSDEPTH (
        ch_mosdepth_input, ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.map{meta, file -> file}.collect())

    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    //
    // MODULE: MOSDEPTH, get coverage on target genes
    //
    ch_mosdepth_target_genes_input_1 = ch_bam.join(ch_bai, by: 0).map { meta, file1, file2 -> tuple(meta, file1, file2) }
    ch_mosdepth_target_genes_input = ch_mosdepth_input_1.join(ch_rep_genes, by: 0).map { meta, file1, file2, bed_file -> tuple(meta, file1, file2, bed_file) }
    MOSDEPTH_TARGET  (
        ch_mosdepth_target_genes_input, ch_fasta
    )

    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_TARGET.out.summary_txt.map{meta, file -> file}.collect())

    //
    // MODULE: Compare coverage with thresholds: writing the results file 
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes
    //
    ch_mosdepth_summary = MOSDEPTH.out.summary_txt
        .join(MOSDEPTH_TARGET.out.regions_bed, by: 0)
        .map { meta, summary, bed -> tuple(meta, summary, bed) }
    
    ch_fastp_mosdepth_merged = ch_fastp_json
        .join(ch_mosdepth_summary, by: 0)
        .map { meta, json, summary, bed -> 
            tuple(meta, json, summary, bed)
    }

    csv_ch = COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged
    )
    
    //
    // MODULE: Merge Reports
    //
    csvs_ch = csv_ch.result_csv.collect()
    MERGE_REPORTS(csvs_ch)
    
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
