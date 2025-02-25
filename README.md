[![GitHub Actions CI Status](https://github.com/nf-core/grzqc/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/grzqc/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BfArM-MVH/GRZ_QC_Workflow** performs extended quality control of GRZ submissions according to the defined thresholds.

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`FASTP`](https://github.com/OpenGene/fastp))
2. Alignment using ([`BWAMEM2`](https://github.com/bwa-mem2/bwa-mem2))
3. Coverage calculation by ([`Mosdepth`](https://github.com/brentp/mosdepth))
4. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))


## Setup

- Install nextflow (and dependencies)
- Make sure to have either conda, docker or singularity.
- Clone the github repository

```bash
git clone https://github.com/BfArM-MVH/GRZ_QC_Workflow.git
$output_path = "path/to/analysis/dir"
```

### Setting up reference files

> [!WARNING]
BWAMEM2 index folder is required for this folder can be either used directly or can be produced through the first run (using grzqc profile) of this pipeline.

- You can also run Download necessary reference fasta files and place into /references directory.

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz $output_path/references
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
mv hg38.fa.gz $output_path/references
```


## Usage

This pipeline needs a samplesheet which is generated automatically from the metadata.json file included in the submission base directory. Please make sure that the submission base directory has the required folder structure. The script run_grzqc.sh parses the metadata.json file to create a nextflow samplesheet:

```bash
python3 metadata_to_samplesheet.py $submission_base_path $output_path
```

Now, you can run the pipeline using:

```bash
nextflow run GRZ_QC_Workflow/main.nf -profile grzqc,docker --outdir $output_path --input $output_path"/grzqc_samplesheet.csv"
```

For your next run, you can use prebuild references. Please prepare your own config file to do so.

## Pipeline output

Output :

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Sample_id`  | Sample id                                                          |
| `libraryType` | ......                                                            |
| `sequenceSubtype` | .......                                                       |
| `genomicStudySubtype` | .......                                                       |
| `q30_rate` | .......                                                       |
| `Q30_THRESHOLD` | .......                                                       |
| `Mosdepth_cov` | .......                                                       |
| `MEAN_DEPTH_THRESHOLD` | .......                                                       |
| `Mosdepth_cov_ratio_target_genes` | .......                                                       |
| `TARGET_FRACTION_ABOVE_THRESHOLD` | .......                                                       |
| `Quality_check` | .......                                                       |


## Credits

This project is under project of ..

## Contributions and Support

nf-core/grzqc was originally written by Shounak Chakraborty, Yun Wang, Kuebra Narci and Florian R. Hölzlwimmer
 .



## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/grzqc for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
>
