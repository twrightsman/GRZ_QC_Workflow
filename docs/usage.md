# Usage

The pipeline input can be specified in one of two ways:

1. `--submission_basepath path/to/submission` (to use a submission folder directly)
2. `--input samplesheet.csv` (to use a manually-prepared samplesheet)

## Samplesheets

Instead of parsing the submission metadata, a CSV samplesheet containing the required information for the pipeline can be provided.
The samplesheet provides a bit more flexibility at the cost of manual preparation.
For example, starting directly from aligned reads is only possible using a samplesheet.

When running the pipeline from unaligned reads each row in the sample sheet is a "run".
Each run consists of either a pair of FASTQ files for paired-end data or a single FASTQ file for single-end data.
A "sample" is equivalent to a "lab datum" in the [official GRZ submission metadata schema](https://github.com/BfArM-MVH/MVGenomseq_GRZ/blob/main/GRZ/grz-schema.json) and can contain multiple runs.

### Multiple runs of the same sample

Runs from the same sample must contain the same sample ID in the `sample` column to tie them together.
The `runId` must be unique to a run within a sample to differentiate them.
Runs in different `sample`s may therefore have the same `runId`.
In short, each row in the samplesheet must have a unique combination of `sample` and `runId`.

Below is an example for the same sample sequenced with 3 runs:

```
sample,runId,libraryType,sequenceSubtype,genomicStudySubtype,reads1,reads2
CONTROL_REP1,run1,wgs,somatic,tumor+germline,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,run2,wgs,somatic,tumor+germline,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,run3,wgs,somatic,tumor+germline,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

A more complete example samplesheet is available [here](../tests/data/test-datasets_sarek/samplesheet.example.csv).

### Example pipeline execution

```bash
nextflow run main.nf \
    -profile conda \
    --outdir "path/to/outdir" \
    --input "samplesheet.csv" \
    --genome "GRCh37"  # or "GRCh38"
```

### Column descriptions

| Column                  | Required                       | Description                                                                                                                                              |
| ----------------------- | ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                | yes                            | Unique sample identifier. Must be identical for multiple sequencing runs from the same sample and not contain spaces.                                    |
| `runId`                 | for unaligned input            | Uniquely identifies a run within a sample. Must not contain spaces.                                                                                      |
| `donorPseudonym`        | no                             | Unique identifier of the donor.                                                                                                                          |
| `laneId`                | no                             | Lane ID of run.                                                                                                                                          |
| `flowcellId`            | no                             | Flowcell ID of run.                                                                                                                                      |
| `labDataName`           | no                             | Sample name or description.                                                                                                                              |
| `libraryType`           | yes                            | `panel`, `wgs`, `wes`, `panel_lr`, `wgs_lr`, or `wes_lr`                                                                                                 |
| `sequenceSubtype`       | yes                            | `somatic` or `germline`                                                                                                                                  |
| `sequencerManufacturer` | no                             | Sequencing platform manufacturer (e.g. Illumina).                                                                                                        |
| `genomicStudySubtype`   | yes                            | `tumor+germline`, `tumor-only`, or `germline-only`                                                                                                       |
| `reads1`                | for unaligned short read input | Full path to FASTQ file for Illumina short reads pair 1 (R1). Must be gzipped and have the extension ".fastq.gz" or ".fq.gz".                            |
| `reads2`                | for unaligned short read input | Full path to FASTQ file for Illumina short reads pair 2 (R2). Must be gzipped and have the extension ".fastq.gz" or ".fq.gz".                            |
| `reads_long`            | for unaligned long read input  | Full path to FASTQ or BAM (PacBio) file for long reads. Must have the extension ".fastq.gz", ".fq.gz", or ".bam".                                        |
| `aligned_reads`         | for aligned input              | Full path to aligned reads (BAM file). Can be used as an alternative to FASTQ reads. See starting from aligned reads section below for more information. |
| `bed_file`              | for WES and panel              | Target region BED for WES and panels with the extension ".bed.gz" or ".bed". Empty for WGS.                                                              |
| `fastp_json`            | no                             | Corresponding FASTP JSON report for `aligned_reads`.                                                                                                     |

## Reference files

You can use reference files with `nextflow run main.nf --reference_path "your/reference/path/references"`. Your reference folder needs a specific directory stucture, with subdirectory of both GRCh37 and GRCh38. In each subdirectory, it contains a genome file, a genome index file and a folder with bwa-mem index. The advantage to use `--reference_path` is that the pipeline can automatically use the right genome reference, you do not have to check the genome version in your GRZ submission beforehand.

```bash
$ tree .
.
├── GRCh37
│   ├── bwamem2
│   │   ├── genome.0123
│   │   ├── genome.amb
│   │   ├── genome.ann
│   │   ├── genome.bwt.2bit.64
│   │   └── genome.pac
│   ├── genome.fa
│   ├── genome.fa.fai
│   └── genome.mmi
└── GRCh38
    ├── bwamem2
    │   ├── genome.0123
    │   ├── genome.amb
    │   ├── genome.ann
    │   ├── genome.bwt.2bit.64
    │   └── genome.pac
    ├── genome.fa
    ├── genome.fa.fai
    └── genome.mmi
```

The other option is to set `--fasta`, `--fai`, `--bwa` individually, or prepare config a file like this:

```bash
    fasta = "your/path/to/reference/GRCh37/genome.fa"
    fai   = "your/path/to/reference/GRCh37/genome.fa.fai"
    bwa   = "your/path/to/reference/GRCh37/bwamem2"
    mmi   = "your/path/to/reference/GRCh37/genome.mmi"
```

You can also set only the genome file with `--fasta <genome file>`. The pipeline will prepare the genome index and bwa index automatically.

Of note, `--fasta`, `--fai`, `--bwa` will only be considered when `--reference_path` is not given.

| Parameters            | Description                                                                                       |
| --------------------- | ------------------------------------------------------------------------------------------------- |
| `save_reference`      | save reference when `--save_reference true` , default false                                       |
| `save_reference_path` | save reference path, default `${outdir}`                                                          |
| `reference_path`      | reference path , default null                                                                     |
| `fasta`               | genome fasta path , only use when reference path is null , default null                           |
| `fai`                 | genome fai path , only use when reference path is null and fasta is also given, default null      |
| `bwa`                 | bwamem index path , only use when reference path is null and fasta is also given , default null   |
| `mmi`                 | minimap2 index path , only use when reference path is null and fasta is also given , default null |

## Starting from aligned reads (BAM)

When running the pipeline directly from aligned reads each row in the samplesheet is a sample, instead of a run.
This assumes that you have merged the alignments from each run already into a single alignment file.

If `fastp_json` is not provided, base quality calculation will be carried out using a [python script](../bin/calculate_basequality.py).

```
sample,aligned_reads,fastp_json
CONTROL_REP1,AEG588A1_S1001.bam,AEG588A1_S1_L002_R2_001.fastp.json
CONTROL_REP2,AEG588A1_S2001.bam,
```

Find a complete example samplesheet [here](../tests/data/test-dataset-alignments/samplesheet_alignment.csv).

Note: running with `--submission_basepath` is not possible when using alignments as input.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run BfArM-MVH/GRZ_QC_Workflow -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull BfArM-MVH/GRZ_QC_Workflow
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [releases page](https://github.com/BfArM-MVH/GRZ_QC_Workflow/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

If you want to use multiple config profiles for various institutional clusters, they are available at nf-core/configs. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation) and download manually to your conf directory.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [documentation](https://www.nextflow.io/docs/stable/index.html) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default this pipeline uses containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
