#!/bin/bash

set -e

submission_basepath="$1"
output_basepath="$2"

# create output directories
mkdir "${output_basepath}/"
mkdir "${output_basepath}/work/"
mkdir "${output_basepath}/grzqc_output/"

# create samplesheet
tmp_samplesheet_path="$(mktemp)"

python3 bin/metadata_to_samplesheet.py \
    "${submission_basepath}" \
    "$tmp_samplesheet_path"

# check if the contents are equal, otherwise copy the new samplesheet to the output directory
cmp --silent "$tmp_samplesheet_path" "${output_basepath}/grzqc_output/grzqc_samplesheet.csv" || \
    cp "$tmp_samplesheet_path" "${output_basepath}/grzqc_output/grzqc_samplesheet.csv"

# run nextflow
nextflow run main.nf \
    -profile grzqc,conda \
    --outdir "${output_basepath}/grzqc_output/" \
    -work-dir "${output_basepath}/work/" \
    --input "${output_basepath}/grzqc_output/grzqc_samplesheet.csv" \
    -resume
