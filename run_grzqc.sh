#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks 1 --cpus-per-task 50
#SBATCH --partition=urgent

submission_base_path=/s/public_webshare/public/grz-example-submissions/WES_tumor+germline/
output_path=/data/ceph/ssd/scratch/tmp/wayu0/grz_qc_workflow/
mkdir -p $output_path"/grzqc_output/"

source /data/nasif12/modules_if12/SL7/i12g/miniforge/24.9.0-0/etc/profile.d/conda.sh

# conda activate nextflow_flo
python3 metadata_to_samplesheet.py $submission_base_path $output_path"/grzqc_output/"

export NXF_HOME=/data/ceph/ssd/scratch/tmp/wayu0/grz_qc_workflow/nextflow_cache_WES

nextflow run main.nf -profile grzqc,conda --outdir $output_path"/grzqc_output/" -work-dir "${output_path}/work/" --input $output_path"/grzqc_output/grzqc_samplesheet.csv" -resume
conda deactivate