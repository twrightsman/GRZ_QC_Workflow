
process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data' :
        'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:cd97ba68cc5b8463' }"

    input:
    path submission_basepath

    output:
    path("*samplesheet.csv"), emit: samplesheet
    env(genome),  emit: genome

    // when:
    // task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metadata_to_samplesheet.py "${submission_basepath}" 
    genome="GrCh37"
    """

}