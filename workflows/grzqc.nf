/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { FASTQ_ALIGN_BWA        } from '../subworkflows/nf-core/fastq_align_bwa/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grzqc_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process ECHO_READS {
 
    debug true
 
    input:
    tuple val(meta), path(reads)
 
    script:
    """
    echo ${reads}
    """
}

process COMPARE_COVERAGE {
 
    debug true
 
    input:
    tuple val(meta), path(summary)
 
    output:
    path('*.result.txt')      , emit: result_txt

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.result.txt
    echo ${meta.id}
    if [[ "${meta.id}" =~ (TUMOR|tumor|Tumor) ]]; then
        if [[ "${meta.id}" =~ (WES|wes) ]]; then
            required_coverage=200  # Example coverage for TUMOR + WES
        elif [[ "${meta.id}" =~ (WGS|wgs) ]]; then
            required_coverage=100  # Example coverage for TUMOR + WGS
        else
            required_coverage=100  # Default for TUMOR
        fi
    elif [[ "${meta.id}" =~ (NORMAL|normal|Normal) ]]; then
        if [[ "${meta.id}" =~ (WES|wes) ]]; then
            required_coverage=60  # Example coverage for NORMAL + WES
        elif [[ "${meta.id}" =~ (WGS|wgs) ]]; then
            required_coverage=30  # Example coverage for NORMAL + WGS
        else
            required_coverage=30  # Default for NORMAL
        fi
    else
        echo "The sample ID does not contain the strings TUMOR or NORMAL." >> ${prefix}.result.txt
    exit 1
    fi

    # Extract the mean coverage from the "total" row of the mosdepth file
    mosdepth_cov=\$(awk 'NR==1 {for (i=1; i<=NF; i++) if (\$i == "mean") col=i} \$1 == "total" {print \$col}' ${summary}) 
    echo "Sample_id,Mosdepth_cov,Required_cov,Qualtiy_check" >> ${prefix}.result.txt
    echo "\$mosdepth_cov >= \$required_coverage" | bc -l
    if (( \$(echo "\$mosdepth_cov >= \$required_coverage" | bc -l) )); then
        echo "${meta.id},\$mosdepth_cov,\$required_coverage,PASS" >> ${prefix}.result.txt
    else
        echo "${meta.id},\$mosdepth_cov,\$required_coverage,FAIL" >> ${prefix}.result.txt
    fi
    """
}

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
                ref_37  : reference == "GRCh37"
                    return [ meta, params.bwa_index_37 ]
                ref_38  : reference == "GRCh38"
                    return [ meta, params.bwa_index_38 ]
        }.set{
        ch_reference_ind
        }

    ch_reference_ind.ref_37.mix(
        ch_reference_ind.ref_38
    ).set{
        ch_bwa_index
    }

    // Create fasta channels for BWA_MEM
    ch_samplesheet.map{
        meta, fastqs, bed_file, reference -> tuple(meta, reference.first())
    }.branch {
            meta, reference ->
                ref_37  : reference == "GRCh37"
                    return [ meta, params.fasta_37 ]
                ref_38  : reference == "GRCh38"
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
    // MODULE: Echo reads
    //
    ECHO_READS (
        ch_cat_fastq
    )

    FASTQ_ALIGN_BWA (
        ch_cat_fastq, ch_bwa_index, true, ch_fasta
    )
    ch_bam = FASTQ_ALIGN_BWA.out.bam
    ch_bai = FASTQ_ALIGN_BWA.out.bai

    //Is there a better way to do this, i.e combine more than 2 channels into one in one step?
    ch_mosdepth_input_1 = ch_bam.combine(ch_bai, by: 0).map { meta, file1, file2 -> tuple(meta, file1, file2) }
    ch_mosdepth_input = ch_mosdepth_input_1.combine(ch_bed_file, by: 0).map { meta, file1, file2, bed_file -> tuple(meta, file1, file2, bed_file) }

    //
    // MODULE: MOSDEPTH
    //
    MOSDEPTH (
        ch_mosdepth_input, ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.map{meta, file -> file}.collect())

    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    //
    // MODULE: Compare coverage
    //
    COMPARE_COVERAGE(
        MOSDEPTH.out.summary_txt
    )

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
