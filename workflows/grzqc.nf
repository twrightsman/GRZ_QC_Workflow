include { PREPARE_REFERENCES            } from '../subworkflows/local/prepare_references'
include { CAT_FASTQ                     } from '../modules/nf-core/cat/fastq'
include { FASTQC                        } from '../modules/nf-core/fastqc'
include { FASTP                         } from '../modules/nf-core/fastp'
include { BWAMEM2_MEM                   } from '../modules/nf-core/bwamem2/mem'
include { PICARD_ADDORREPLACEREADGROUPS } from '../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES         } from '../modules/nf-core/picard/markduplicates/main'
include { CONVERT_BED_CHROM             } from '../modules/local/convert_bed_chrom'
include { MOSDEPTH                      } from '../modules/nf-core/mosdepth'
include { COMPARE_THRESHOLD             } from '../modules/local/compare_threshold'
include { MERGE_REPORTS                 } from '../modules/local/merge_reports'
include { paramsSummaryMap              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MULTIQC                       } from '../modules/nf-core/multiqc'

workflow GRZQC {

    take:
    samples
    genomes
    reference_path
    thresholds

    main:
    thresholds = channel.fromPath(thresholds, checkIfExists: true)
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    PREPARE_REFERENCES(
      samples,
      genomes,
      reference_path
    )

    // Concatenate multiple read files from same sample if required
    CAT_FASTQ(samples)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    FASTQC(CAT_FASTQ.out.reads)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ meta, zip -> zip })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    FASTP(
      CAT_FASTQ.out.reads,
      [],    // adapter sequences
      false, // discard passed reads
      false, // save failed reads
      false  // save reads to merged file
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map{ meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.map{ meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions)

    BWAMEM2_MEM(
      FASTP.out.reads,
      PREPARE_REFERENCES.out.bwamem2
        .map{ meta, index -> [meta.id, index] }
        .combine(FASTP.out.reads.map{ meta, reads -> [meta.reference] }, by: 0)
        .map{ id, index -> [[id: id], index] },
      channel.value([[:], []]),  // no reference needed for BAM output
      true  // sort alignments
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    PICARD_ADDORREPLACEREADGROUPS(
      BWAMEM2_MEM.out.bam,
      channel.value([[:], []]),  // no reference needed for BAM input
      channel.value([[:], []])   // ^
    )
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions)

    PICARD_MARKDUPLICATES(
      PICARD_ADDORREPLACEREADGROUPS.out.bam,
      channel.value([[:], []]),  // no reference needed for BAM input
      channel.value([[:], []])   // ^
    )
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.map{ meta, metrics -> metrics })
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    CONVERT_BED_CHROM(
      PICARD_MARKDUPLICATES.out.bam
        .map{ meta, bam -> [ meta, meta.targets ] },
      PICARD_MARKDUPLICATES.out.bam
        .map{ meta, bam -> [ meta.reference ] }
        .combine(
          channel.of(genomes)
            .flatMap{ it.entrySet() }
            .map{ [it.getKey(), it.getValue() ] },
          by: 0
        )
        .map{ id, paths -> file(paths.name_map) }
    )
    ch_versions = ch_versions.mix(CONVERT_BED_CHROM.out.versions)

    MOSDEPTH(
      PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai)
        .join(CONVERT_BED_CHROM.out.converted_bed),
      channel.value([[:], []])   // no reference needed for BAM input
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.map{ meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.map{ meta, file -> file })

    COMPARE_THRESHOLD(
      MOSDEPTH.out.regions_bed
        .join(MOSDEPTH.out.summary_txt)
        .join(FASTP.out.json),
      thresholds
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)

    MERGE_REPORTS(
      COMPARE_THRESHOLD.out.result_csv.collect()
    )
    ch_versions = ch_versions.mix(MERGE_REPORTS.out.versions)

    //
    // Collate and save software versions
    //
    ch_multiqc_files = ch_multiqc_files.mix(
      softwareVersionsToYAML(ch_versions)
        .collectFile(
          storeDir: "${params.outdir}/pipeline_info",
          name: 'pipeline_software_mqc_versions.yml',
          sort: true,
          newLine: true
        )
    )

    MULTIQC(
      ch_multiqc_files.collect(),
      Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true),
      [], // extra config
      [], // custom logo
      [], // replace names
      []  // sample names
    )
}
