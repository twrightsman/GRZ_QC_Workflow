include { SAMTOOLS_BGZIP } from '../../../modules/nf-core/samtools/bgzip'
include { BWAMEM2_INDEX  } from '../../../modules/nf-core/bwamem2/index'

workflow PREPARE_REFERENCES {
  take:
  samples
  genomes
  reference_path

  main:
  def necessary_references = samples
    .map{ meta, reads -> meta.reference }
    .unique()

  def reference_fasta_state = necessary_references
    .branch{ id ->
      existent: reference_path ? (file(reference_path) / id / "${id}.fasta.gz").exists() : false
      nonexistent: true
    }

  fastas_existent = reference_fasta_state.existent
    .map{ id -> [[id: id], file(reference_path) / id / "${id}.fasta.gz" ] }

  reference_fasta_state.nonexistent
    .join(channel.of(genomes)
      .flatMap{ it.entrySet() }
      .map{ [it.getKey(), it.getValue()] }
    )
    .map { id, paths -> [[id: id], file(paths.fasta)] }
    | SAMTOOLS_BGZIP

  fastas_nonexistent = SAMTOOLS_BGZIP.out.fasta

  fastas = fastas_existent.mix(fastas_nonexistent)

  def bwamem2_index_state = necessary_references
    .branch{ id ->
      existent: reference_path ? (file(reference_path) / id / 'bwamem2/').exists() : false
      nonexistent: true
    }

  bwamem2_indices_existent = bwamem2_index_state.existent
    .map{ id -> [[id: id], file(reference_path) / id / 'bwamem2/' ] }

  BWAMEM2_INDEX(
    bwamem2_index_state
      .nonexistent
      .map{ id -> [id: id] }
      .join( fastas )
  )

  emit:
  bwamem2 = bwamem2_indices_existent.mix(BWAMEM2_INDEX.out.index)

}
