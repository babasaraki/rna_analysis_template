# RNA analysis

## Overview of the analysis

The RNA-seq data was analysed with [rnaseq version 3.5](https://github.com/nf-core/rnaseq/tree/3.5). The outputs of these pipelines were then analysed with the [rna_analysis_template 1.0.0](https://github.com/leahkemp/rna_analysis_template/tree/1.0.0). The following treatment groups were compared for analysis:

- treatment1 - treatment2

The full analysis has been documented so others can take a "deep dive" into the analysis to check/reproduce the work.

## Quality control reports

### rnaseq

- [MultiQC report](./test/rnaseq_pipeline_run/results/multiqc/star_salmon/multiqc_report.html)

## RNA counts/expression

- [RNA expression plotting](https://esr-cri.shinyapps.io/rnaseq_expression_plotting_example/)

## MDS/PCA plotting

- [MDS plots](./example_webpage/mds.html)
- [PCA plots](https://esr-cri.shinyapps.io/rnaseq_pca_example/)
  
## Differential expression

- [Differential expression (limma/voom)](./example_webpage/diff_expression_limma_voom.html)
- [Differential expression (deseq2)](./example_webpage/diff_expression_deseq.html)
- [Differential expression - all results](./example_webpage/diff_expression_both.html)

## Heatmaps

- [Heatmaps](./example_webpage/heatmaps.html)

## Presence/absence

- [Presence/absence](./example_webpage/presence_absence.html)

## Raw count data

- [rnaseq transcripts](./test/rnaseq_pipeline_run/results/star_salmon/salmon.merged.transcript_counts.tsv)
- [rnaseq genes](./test/rnaseq_pipeline_run/results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv)

## Analysis workflow for reproducing this analysis

- [master.Rmd](./master.Rmd)
