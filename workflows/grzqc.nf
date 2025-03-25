/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grzqc_pipeline'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq'
include { FASTQC                 } from '../modules/nf-core/fastqc'
include { FASTP                  } from '../modules/nf-core/fastp'
include { MULTIQC                } from '../modules/nf-core/multiqc' 
include { CONVERT_BED_CHROM      } from '../modules/local/convert_bed_chrom'
include { COMPARE_THRESHOLD      } from '../modules/local/compare_threshold'
include { MERGE_REPORTS          } from '../modules/local/merge_reports'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_37        } from '../modules/nf-core/samtools/faidx/main'    
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_38         } from '../modules/nf-core/samtools/faidx/main'   
include { FASTQ_ALIGN_BWA_MARKDUPLICATES as FASTQ_ALIGN_BWA_MARKDUPLICATES_HG38 } from '../subworkflows/local/fastq_align_bwa_markduplicates'
include { FASTQ_ALIGN_BWA_MARKDUPLICATES as FASTQ_ALIGN_BWA_MARKDUPLICATES_HG37 } from '../subworkflows/local/fastq_align_bwa_markduplicates'
include { BWAMEM2_INDEX  as  BWAMEM2_INDEX_HG38   } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX  as  BWAMEM2_INDEX_HG37   } from '../modules/nf-core/bwamem2/index/main'
include { MOSDEPTH as MOSDEPTH_HG38 } from '../modules/nf-core/mosdepth'
include { MOSDEPTH as MOSDEPTH_HG37 } from '../modules/nf-core/mosdepth'

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

    // create reference channels
    fasta_37       = Channel.fromPath(params.fasta_37, checkIfExists: true)
                    .map{ file -> tuple([id: file.getSimpleName()], file) }.collect()
    fasta_38        = Channel.fromPath(params.fasta_38, checkIfExists: true)
                    .map{ file -> tuple([id: file.getSimpleName()], file) }.collect()

    bwa_index_37    = params.bwa_index_37   ? Channel.fromPath(params.bwa_index_37 , checkIfExists: true, type: 'dir').map{ file -> tuple([id: "hg19"], file) }.collect()
                                            : Channel.empty()

    bwa_index_38    = params.bwa_index_38   ? Channel.fromPath(params.bwa_index_38 , checkIfExists: true, type: 'dir').map{ file -> tuple([id: "hg38"], file) }.collect()
                                            : Channel.empty()

    if (!params.thresholds){
        thresholds  = Channel.fromPath("${projectDir}/assets/default_files/thresholds.json", checkIfExists: true).collect()
    } else {
        thresholds  = Channel.fromPath(params.thresholds, checkIfExists: true).collect()
    }
    if (!params.rep_genes_hg38){
        rep_genes_38  = Channel.fromPath("${projectDir}/assets/default_files/hg38_440_omim_genes.bed", checkIfExists: true).collect()
    } else {
        rep_genes_38  = Channel.fromPath(params.rep_genes_hg38, checkIfExists: true).collect()
    }
    if (!params.rep_genes_hg19){
        rep_genes_37  = Channel.fromPath("${projectDir}/assets/default_files/hg19_439_omim_genes.bed", checkIfExists: true).collect()
    } else {
        rep_genes_37  = Channel.fromPath(params.rep_genes_hg19, checkIfExists: true).collect()
    }

    chrom_mapping_37 = Channel.value(file("${projectDir}/assets/default_files/hg19_NCBI2UCSC.txt"))
    chrom_mapping_38 = Channel.value(file("${projectDir}/assets/default_files/hg38_NCBI2UCSC.txt"))


    // Creating merge fastq channel from samplesheet
    ch_samplesheet.map{
        meta, fastqs, bed_file -> tuple(meta, fastqs)
    }.branch {
            meta, fastqs ->
                single_sample_paired_end  : fastqs.size() < 2 | fastqs.size() == 2
                    return [ meta, fastqs.flatten() ]
                multiple_sample_paired_end: fastqs.size() > 2
                    return [ meta, fastqs.flatten() ]
        }.set{
            ch_fastqs
        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastqs.multiple_sample_paired_end
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    CAT_FASTQ.out.reads
        .mix(ch_fastqs.single_sample_paired_end)
        .set { ch_cat_fastq }


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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

    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{ meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{ meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // branch out samples per reference
    ch_cat_fastq.branch{
        def meta = it[0]
        hg38: meta.reference == "GRCh38" || meta.reference == "hg38"
        hg37: meta.reference == "GRCh37" || meta.reference == "hg19"
        other:false
    }.set{input_samples}

    ch_bams = Channel.empty()
    // bwa index creation might be dropped

    if (!params.bwa_index_38){
        BWAMEM2_INDEX_HG38(fasta_38)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_HG38.out.versions.first())
        bwa_index_38 = BWAMEM2_INDEX_HG38.out.index
    }

    // Create sequence fai files
    SAMTOOLS_FAIDX_38(fasta_38,[[],[]],false)
    fasta_38_fai = SAMTOOLS_FAIDX_38.out.fai 
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_38.out.versions.first())

    // hg38 alignment analysis
    FASTQ_ALIGN_BWA_MARKDUPLICATES_HG38 (
        input_samples.hg38,
        bwa_index_38,
        true,
        fasta_38,
        fasta_38_fai
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG38.out.versions.first())

    // bwa index creation might be dropped
    if (!params.bwa_index_37){
        BWAMEM2_INDEX_HG37(fasta_37)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_HG38.out.versions.first())
        bwa_index_37 = BWAMEM2_INDEX_HG37.out.index
    }

    // Create sequence fai files
    SAMTOOLS_FAIDX_37(fasta_37,[[],[]],false)
    fasta_37_fai = SAMTOOLS_FAIDX_37.out.fai
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_37.out.versions.first())

    // hg37 alignment analysis
    FASTQ_ALIGN_BWA_MARKDUPLICATES_HG37 (
        input_samples.hg37,
        bwa_index_37,
        true,
        fasta_37,
        fasta_37_fai
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG37.out.versions.first())

    ch_bams = ch_bams
                .mix(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG37.out.bam
                    .join(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG37.out.bai, by:0))

    ch_bams = ch_bams
                .mix(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG38.out.bam
                    .join(FASTQ_ALIGN_BWA_MARKDUPLICATES_HG38.out.bai, by:0))

    // prepare mosdepth inputs               
    // Prepare bed files 
    // for WGS: defined bed files of ~400 genes
    // for WES and panel: extract bed from samplesheet and 
    //      do the conversion with the correct mapping file (if it is NCBI format, covert it to UCSC).
    // for both cases different file required for hg38 and hg19 
    
    ch_samplesheet_target = ch_samplesheet.filter { meta, fastqs, bed_file ->
        meta.libraryType in ["wes", "panel", "wes_lr", "panel_lr"]
    }
    ch_samplesheet_wgs = ch_samplesheet.filter { meta, fastqs, bed_file ->
        meta.libraryType in ["wgs", "wgs_lr"]
    }

    // --- For target samples: prepare bed for conversion ---
    ch_samplesheet_target.map { meta, fastqs, bed_file ->
        def bed = bed_file.first()
        def mapping = (meta.reference == "GRCh38" || meta.reference == "hg38") ? chrom_mapping_38.getVal():
                      ((meta.reference == "GRCh37" || meta.reference == "hg19") ? chrom_mapping_37.getVal(): null)
        return tuple(meta, bed, mapping)
    }.set{ch_bed_for_conversion_target}

    // for WES and panel, run the conversion process: if the bed file has NCBI-style names, they will be converted.
    CONVERT_BED_CHROM(
        ch_bed_for_conversion_target
    )
    ch_converted_bed = CONVERT_BED_CHROM.out.converted_bed
    ch_versions = ch_versions.mix(CONVERT_BED_CHROM.out.versions)

    // Prepare mosdepth inputs for target samples
    ch_mosdepth_input_target = ch_samplesheet_target
        .map { meta, fastqs, bed_file -> tuple(meta, bed_file.first()) }
        .join(ch_bams, by: 0)
        .join(ch_converted_bed, by: 0)
        .map { joined ->
            def meta = joined[0]
            // original bed is ignored; we use the converted version
            def bam = joined[2]
            def bai = joined[3]
            def converted_bed = joined[4]
            return tuple(meta, bam, bai, converted_bed)
        }

    // --- Prepare mosdepth inputs for WGS samples ---
    ch_mosdepth_input_wgs = ch_samplesheet_wgs
        .map { meta, fastqs, bed_file -> tuple(meta, bed_file.first()) }
        .join(ch_bams, by: 0)
        .map { joined ->
            def meta = joined[0]
            def bam = joined[2]
            def bai = joined[3]
            // For WGS, we use rep_genes of ~400 genes:
            def rep_genes = (meta.reference == "GRCh38" || meta.reference == "hg38") ? rep_genes_38.getVal() :
                            ((meta.reference == "GRCh37" || meta.reference == "hg19") ? rep_genes_37.getVal() : null)
            return tuple(meta, bam, bai, rep_genes)
        }

    // --- Merge the two sets of mosdepth inputs ---
    ch_mosdepth_input_target
        .mix(ch_mosdepth_input_wgs)
        .branch {
            def meta = it[0]
            hg38: meta.reference == "GRCh38" || meta.reference == "hg38"
            hg37: meta.reference == "GRCh37" || meta.reference == "hg19"
            other: false
        }.set{ch_mosdepth_input}

    // hg38 mosdepth analysis
    MOSDEPTH_HG38(
        ch_mosdepth_input.hg38,
        fasta_38,
    )
    ch_versions = ch_versions.mix(MOSDEPTH_HG38.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_HG38.out.summary_txt.map{meta, file -> file}.collect())

    // hg37 mosdepth analysis
    MOSDEPTH_HG37(
        ch_mosdepth_input.hg37,
        fasta_37,
    )
    ch_versions = ch_versions.mix(MOSDEPTH_HG37.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_HG37.out.summary_txt.map{meta, file -> file}.collect())

    // collect the results for comparison

    MOSDEPTH_HG37.out.summary_txt.mix(
        MOSDEPTH_HG38.out.summary_txt
    ).set{ch_mosdepth_summary}

    MOSDEPTH_HG37.out.regions_bed.mix(
        MOSDEPTH_HG38.out.regions_bed
    ).set{ch_mosdepth_target_regions}

    ch_mosdepth_summary.join(ch_mosdepth_target_regions, by:0)
        .join(FASTP.out.json, by:0)
        .map{meta, summary, json, bed -> tuple(meta, json, summary,bed)}
        .set{ch_fastp_mosdepth_merged}

    //
    // MODULE: Compare coverage with thresholds: writing the results file
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes
    //
    COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged,
        thresholds
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)

    //
    // MODULE: Merge Reports
    //
    MERGE_REPORTS(COMPARE_THRESHOLD.out.result_csv.collect())
    ch_versions = ch_versions.mix(MERGE_REPORTS.out.versions)

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
