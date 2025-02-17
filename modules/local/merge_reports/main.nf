//
// Collect all result files produced by COMPARE_THRESHOLD into an array
//
process MERGE_REPORTS {
    publishDir "${params.outdir}/results", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/pip_openpyxl_pandas:abdd97f3f86df268"

    input:
      path csv_files 

    output:
      path "merged_result.csv"
      path "merged_result.xlsx"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    # concat all csv files
    files = "${csv_files.join(' ')}".split()
    dfs = [pd.read_csv(f) for f in files]
    df_merged = pd.concat(dfs, ignore_index=True)
    df_merged.to_csv("merged_result.csv", index=False)
    df_merged.to_excel("merged_result.xlsx", index=False)
    """
}
