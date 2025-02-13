
process COMPARE_THRESHOLD {
    publishDir "${params.outdir}/results", mode: 'copy'
    debug true
    
    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/pip_gzip-utils_json-advanced_openpyxl_pandas:96fa2b102489abbc"
  
    input:
    tuple val(meta), path(fastp_json), path(summary), path(bed)
 
    output:
    path('*.result.csv')      , emit: result_csv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3  
    import pandas as pd
    import gzip
    import json

    sample_id = "${meta.id}"
    prefix = "${prefix}"
    libraryType = "${meta.libraryType}"
    sequenceSubtype = "${meta.sequenceSubtype}"
    genomicStudySubtype = "${meta.genomicStudySubtype}"

    ### Extract results and thresholds
    # Extract Q30 rate from the fastp JSON.
    with open("${fastp_json}", "r") as f:
        fastp_data = json.load(f)
    q30_rate = fastp_data["summary"]["before_filtering"]["q30_rate"]

    # Extract thresholds JSON.
    with open("${params.thresholds_json}", "r") as f:
        thresholds_data = json.load(f)

    THRESHOLDS = None
    for item in thresholds_data:
        if (item["libraryType"] == libraryType and
            item["sequenceSubtype"] == sequenceSubtype and
            item["genomicStudySubtype"] == genomicStudySubtype):
            THRESHOLDS = item["thresholds"]
            break

    if THRESHOLDS is None:
        raise ValueError("No matching thresholds found for meta values: " + libraryType + ", " + sequenceSubtype + ", " + genomicStudySubtype)

    MEAN_DEPTH_THRESHOLD = float(THRESHOLDS["meanDepthOfCoverage"])
    QUALITY_THRESHOLD = THRESHOLDS["fractionBasesAboveQualityThreshold"]["qualityThreshold"]

    if str(QUALITY_THRESHOLD) == "30":
        Q30_THRESHOLD = float(THRESHOLDS["fractionBasesAboveQualityThreshold"]["fractionBasesAbove"])

    TARGET_MIN_COVERAGE = float(THRESHOLDS["targetedRegionsAboveMinCoverage"]["minCoverage"])
    TARGET_FRACTION_ABOVE_THRESHOLD = float(THRESHOLDS["targetedRegionsAboveMinCoverage"]["fractionAbove"])

    ### compare results with thresholds

    mosdepth_summary = "${summary}"
    bed_file = "${bed}"

    # 1. Read mosdepth summary file
    df = pd.read_csv(mosdepth_summary, sep="\\t")
    row_name = "total_region" if libraryType in ["panel", "wes"] else "total"
    mosdepth_cov = float(df.loc[df.iloc[:,0] == row_name, "mean"].values[0])

    # 2. parse the mosdepth target gene result 
    count = 0
    total = 0
    with gzip.open(bed_file, "rt") as f:
        for line in f:
            total += 1
            parts = line.strip().split()
            if float(parts[3]) >= TARGET_MIN_COVERAGE:
                count += 1

    mosdepth_cov_rate_target = count / total if total > 0 else 0
    ##    # per base solution
    ##    # ! note if we change this , we also need to change the output channel of MOSDEPTH_TARGET
    ##    total_length = 0
    ##    filtered_length = 0
    ##    with gzip.open(bed_file, "rt") as f:
    ##        for line in f:
    ##            parts = line.strip().split()
    ##            # Calculate interval length: col3 - col2
    ##            start = float(parts[1])
    ##            end   = float(parts[2])
    ##            interval_length = end - start
    ##            total_length += interval_length
    ##            # If column 4 > TARGET_MIN_COVERAGE, add to filtered length.
    ##            if float(parts[3]) > TARGET_MIN_COVERAGE:
    ##                filtered_length += interval_length
    ##
    ##    mosdepth_cov_rate_target = filtered_length / total_length if total_length > 0 else 0

    # 3. Perform the quality check.
    quality_check = "PASS" if (mosdepth_cov >= MEAN_DEPTH_THRESHOLD and q30_rate >= Q30_THRESHOLD and mosdepth_cov_rate_target >= TARGET_FRACTION_ABOVE_THRESHOLD) else "FAIL"

    # 4. Write the results to a CSV file.
    with open(f"{prefix}.result.csv", "w") as f:
        f.write("Sample_id, libraryType, sequenceSubtype, genomicStudySubtype, q30_rate, Q30_THRESHOLD, Mosdepth_cov, MEAN_DEPTH_THRESHOLD, Mosdepth_cov_ratio_target_genes, TARGET_FRACTION_ABOVE_THRESHOLD, Quality_check\\n")
        f.write(f"{sample_id},{libraryType},{sequenceSubtype},{genomicStudySubtype},{q30_rate},{Q30_THRESHOLD},{mosdepth_cov},{MEAN_DEPTH_THRESHOLD},{mosdepth_cov_rate_target},{TARGET_FRACTION_ABOVE_THRESHOLD},{quality_check}\\n")

    """
}
