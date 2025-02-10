
process COMPARE_THRESHOLD {
    publishDir "${params.outdir}/results", mode: 'copy'
    debug true
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"
 
    input:
    tuple val(meta), path(fastp_json), path(summary), path(bed)
 
    output:
    path('*.result.csv')      , emit: result_csv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.result.csv
    echo ${meta.id}

    # Extract Q30 rate from the fastp JSON
    q30_rate=\$(jq '.summary.before_filtering.q30_rate' ${fastp_json})
    
    # Extract the required coverage from the thresholds JSON
    THRESHOLDS=\$(jq --arg lt '${meta.libraryType}' \
                     --arg st '${meta.sequenceSubtype}' \
                     --arg gst '${meta.genomicStudySubtype}' \
                     '.[] | select(.libraryType == \$lt and .sequenceSubtype == \$st and .genomicStudySubtype == \$gst) | .thresholds' '${params.thresholds_json}')

    MEAN_DEPTH_THRESHOLD=\$(echo \$THRESHOLDS | jq -r '.meanDepthOfCoverage')
    QUALITY_THRESHOLD=\$(echo \$THRESHOLDS | jq -r '.fractionBasesAboveQualityThreshold.qualityThreshold')
    if [[ \$QUALITY_THRESHOLD == "30" ]]; then
        Q30_THRESHOLD=\$(echo \$THRESHOLDS | jq -r '.fractionBasesAboveQualityThreshold.fractionBasesAbove')
    fi
    TARGET_MIN_COVERAGE=\$(echo \$THRESHOLDS | jq -r '.targetedRegionsAboveMinCoverage.minCoverage')
    TARGET_FRACTION_ABOVE_THRESHOLD=\$(echo \$THRESHOLDS | jq -r '.targetedRegionsAboveMinCoverage.fractionAbove')

    # For WGS, extract the mean coverage from the "total" row of the mosdepth file
    # For WES, extract the mean coverage from the "total_region" row of the mosdepth file
    if [[ "${meta.libraryType}" =~ (panel|wes) ]]; then
        mosdepth_cov=\$(awk 'NR==1 {for (i=1; i<=NF; i++) if (\$i == "mean") col=i} \$1 == "total_region" {print \$col}' ${summary}) 
    else
        mosdepth_cov=\$(awk 'NR==1 {for (i=1; i<=NF; i++) if (\$i == "mean") col=i} \$1 == "total" {print \$col}' ${summary}) 
    fi

    mosdepth_cov_rate_target=\$(gunzip -c ${bed} | awk -v req="\${TARGET_MIN_COVERAGE}" '
        BEGIN { count=0; total=0 }
        { total++; if (\$4 >= req) count++ }
        END { if (total > 0) print count/total; else print 0 }
    ' )
      
    echo "Sample_id, libraryType, sequenceSubtype, genomicStudySubtype, q30_rate,Q30_THRESHOLD,Mosdepth_cov,MEAN_DEPTH_THRESHOLD,Mosdepth_cov_ratio_target_genes,TARGET_FRACTION_ABOVE_THRESHOLD,Qualtiy_check" >> ${prefix}.result.csv
       
    if (( \$(echo "(\$mosdepth_cov >= \$MEAN_DEPTH_THRESHOLD)*(\$q30_rate >= \$Q30_THRESHOLD)*(\$mosdepth_cov_rate_target >= \$TARGET_FRACTION_ABOVE_THRESHOLD)" | bc -l) )); then
        echo "${meta.id},${meta.libraryType},${meta.sequenceSubtype},${meta.genomicStudySubtype},\$q30_rate,\$Q30_THRESHOLD,\$mosdepth_cov,\$MEAN_DEPTH_THRESHOLD,\$mosdepth_cov_rate_target,\$TARGET_FRACTION_ABOVE_THRESHOLD,PASS" >> ${prefix}.result.csv
    else
        echo "${meta.id},${meta.libraryType},${meta.sequenceSubtype},${meta.genomicStudySubtype},\$q30_rate,\$Q30_THRESHOLD,\$mosdepth_cov,\$MEAN_DEPTH_THRESHOLD,\$mosdepth_cov_rate_target,\$TARGET_FRACTION_ABOVE_THRESHOLD,FAIL" >> ${prefix}.result.csv
    fi
    """
}
