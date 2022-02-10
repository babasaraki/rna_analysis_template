# RNA analysis template

---

Analysis template for analysing, visualising and communicating the findings of [Illumina RNA sequencing](https://www.illumina.com/techniques/sequencing/rna-sequencing.html) data. This is set up to analyse human data, but can be adapted to other species with some tweaks to the code, namely the pipeline inputs. This template uses the count data and other QC from both the [rnaseq](https://github.com/nf-core/rnaseq/) and [parabricks rna](https://docs.nvidia.com/clara/parabricks/v3.6.1/text/rna_pipeline.html) pipelines, so I'm assuming you've run these pipelines on your data before using this template. See a live example github page [here](https://leahkemp.github.io/rna_analysis_template/).

---

## Table of contents

- [RNA analysis template](#rna-analysis-template)
  - [Table of contents](#table-of-contents)
  - [What this template can do](#what-this-template-can-do)
  - [What this template can't do](#what-this-template-cant-do)
  - [What's this template gonna do?](#whats-this-template-gonna-do)
  - [Prerequisite software](#prerequisite-software)
  - [Testing](#testing)
  - [How to use this template](#how-to-use-this-template)
    - [1. Fork the template repo to a personal or lab account](#1-fork-the-template-repo-to-a-personal-or-lab-account)
    - [2. Take this template to the data on your local machine](#2-take-this-template-to-the-data-on-your-local-machine)
    - [3. Format your input files](#3-format-your-input-files)
      - [3.1 Fastq naming convention](#31-fastq-naming-convention)
      - [3.2 Metadata file](#32-metadata-file)
      - [3.3 Configuration file](#33-configuration-file)
    - [4. Run the template](#4-run-the-template)
    - [5. Commit and push to your forked version of the github repo](#5-commit-and-push-to-your-forked-version-of-the-github-repo)
    - [6. Repeat step 5 each time you re-run the analysis with different parameters](#6-repeat-step-5-each-time-you-re-run-the-analysis-with-different-parameters)
    - [7. Create a github page (optional)](#7-create-a-github-page-optional)
    - [8. Contribute back!](#8-contribute-back)

## What this template can do

This template uses open source tools and includes several scripts for researchers to analyse, explore and communicate findings to other researchers through interactive tables and plots. The results of this template can be served as a [github page](https://pages.github.com/) that renders html files and provides links to RShiny apps hosted on [shinyapps.io](https://www.shinyapps.io/) - this means a single weblink can be given to your collaborators to provide them with all your analysis code and results. Most importantly, this template puts you in good steed to ensure your analysis is reproducible!

## What this template can't do

- Tell you what analysis tools and parameters are appropriate for your data or research question, the assumption is that the tools this template uses are tools you've intentionally chosen to use and that you will actively adapt this template for your use-case
- Account for different operating systems and compute infrastructures - this means some UNIX experience may be required to run the pipelines/scripts on your operating system or job scheduler. I won't tell you how to do this here, but the pipelines and tools used here are generally portable (ie. able to be run on different operating systems) and I've used [renv](https://rstudio.github.io/renv/articles/renv.html) environments to make the R code more portable

## What's this template gonna do?

Beyond the QC the pipelines undertake, additional QC is undertaken to summarise the read counts and mapping rates of the data.

These counts datasets output by the pipelines are analysed in [R](https://www.r-project.org/) to undertake a differential expression analysis of all these RNA species to find differently expressed RNA's. Two methods were employed to undertake a differential expression analysis, namely [limma/voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29) and [deseq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).

Beyond a traditional differential expression analysis, the data is prepared and presented in an interactive [RShiny](https://shiny.rstudio.com/) app that allows the user to explore RNA expression (both raw counts and counts per million).

Interactive MDS and PCA plots are also created to explore clusters of RNA's/samples in the data.

Lastly, the composition of the RNA species that are identified in each treatment group are compared.

## Prerequisite software

[python 3](https://www.python.org/), [cut](https://www.man7.org/linux/man-pages/man1/cut.1.html), [zgrep](https://linux.die.net/man/1/zgrep), [GNU software](https://www.gnu.org/software/), [pandas](https://pypi.org/project/pandas/), [PyYAML](https://pypi.org/project/PyYAML/), [R](https://www.r-project.org/), [git](https://git-scm.com/)

## Testing

This template has been validated to work on:

- Pipeline outputs created using [nextflow 21.10.6](https://github.com/nextflow-io/nextflow/tree/v21.10.6), [singularity 3.7.2](https://github.com/hpcng/singularity/tree/v3.7.2), [rnaseq version 3.5](https://github.com/nf-core/rnaseq/tree/3.5), [parabricks rna pipeline version 3.6.1](https://docs.nvidia.com/clara/parabricks/v3.6.1/text/rna_pipeline.html)
- R version 4.0.5
- CentOS Linux 7

Test data available in the [test directory](./test/), including [fastq data](./test/fastq/), associated [metadata](./test/metadata.csv), [smrnaseq pipeline outputs](./test/smrnaseq_pipeline_run/) and [excerpt pipeline outputs](./test/excerpt_pipeline_run/) run on this test fastq data

## How to use this template

### 1. Fork the template repo to a personal or lab account

### 2. Take this template to the data on your local machine

Clone the forked [rna_analysis_template](https://github.com/leahkemp/rna_analysis_template) repo to the machine you'd like to analyse the data on

```bash
git clone https://github.com/leahkemp/rna_analysis_template.git
```

### 3. Format your input files

#### 3.1 Fastq naming convention

```bash
sample.fastq.gz
```

- one fastq file per sample
- sample name matching the sample names in the metadata file
- ".fastq.gz" extension

For example see the test fastq files [here](./test/fastq/)

#### 3.2 Metadata file

Required columns:

- "sample"
  - must be titled "sample"
  - must contain a row with a unique sample name/id for each fastq file present in the directory of fastq files to be analysed

- "treatment"
  - must be titled "treatment"

Other notes:

- you can have additional columns
- you can't have any duplicate column names eg. two columns named "sample" and "Sample"
- make sure every sample in the fastq directory to be analysed is included in the metadata file and is associated with a treatment group

For example see the test metadata file [here](./config/metadata.csv)

#### 3.3 Configuration file

Set up [./config/config.yaml](config/config.yaml)

### 4. Run the template

Run/work through the [master RMarkdown file](./master.Rmd), this will do the bulk of the analyses and generate several html file for data visualisation and csv files with processed data

### 5. Commit and push to your forked version of the github repo

To maintain reproducibility of your analysis, commit and push:

- All configuration files
- All run scripts
- All your documentation/notes

### 6. Repeat step 5 each time you re-run the analysis with different parameters

### 7. Create a github page (optional)

Push all the all the html/Rmd/txt files that you'd like to be rendered into a github page to github. Then create a [github page](https://guides.github.com/features/pages/). See this example [live page](https://leahkemp.github.io/rna_analysis_template/) and the underlying [index.md file](index.md) that links all the html/Rmd/txt files included in the github page based on the files produced by this template.

### 8. Contribute back!

- Raise issues in [the issues page](https://github.com/leahkemp/rna_analysis_template/issues)
- Create feature requests in [the issues page](https://github.com/leahkemp/rna_analysis_template/issues)
- Ask questions or start general discussion about this repo in [the discussions page](https://github.com/leahkemp/rna_analysis_template/discussions)
- Contribute your code! Create your own branch from the [development branch](https://github.com/leahkemp/rna_analysis_template/tree/dev) and create a pull request to the [development branch](https://github.com/leahkemp/rna_analysis_template/tree/dev) once the code is on point!

Contributions and feedback are always welcome! :blush:
