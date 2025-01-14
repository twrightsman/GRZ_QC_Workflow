#!/bin/bash


submission_base_path=$1
output_path=$2
mkdir -p $output_path"/grzqc_output/"


python metadata_to_samplesheet.py $submission_base_path $output_path"/grzqc_output/"

nextflow run main.nf -profile docker,grzqc --outdir $output_path"/grzqc_output/" --input $output_path"/grzqc_output/grzqc_samplesheet.csv" -resume
