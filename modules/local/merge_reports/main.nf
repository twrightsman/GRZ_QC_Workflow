//
// Collect all result files produced by COMPARE_THRESHOLD into an array
//
process MERGE_REPORTS {
    publishDir "${params.outdir}/results", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"
 
    input:
      path csv_files 

    output:
      path "merged_result.csv"
      path "merged_result.xlsx"

    script:
    """
    # Use awk to skip the header of every file except the first
    awk 'FNR==1 && NR!=1 {next} 1' ${csv_files.join(' ')} > merged_result.csv
    
    # Convert the CSV output to Excel 
    python3 - <<EOF
import pandas as pd
df = pd.read_csv("merged_result.csv")
df.to_excel("merged_result.xlsx", index=False)
EOF
    """
}
