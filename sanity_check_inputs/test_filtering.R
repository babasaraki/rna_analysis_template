# load libraries
library(dplyr)
library(limma)
library(DESeq2)
library(janitor)
library(edgeR)
library(DT)
library(apeglm)
library(plotly)
library(heatmaply)
library(gtools)
library(textshape)
library(tidyr)

# read in yaml config file
config <- yaml::yaml.load_file("./config/config.yaml")

# read in metadata
metadata <- utils::read.csv(base::file.path(config$metadata))

# load all the count datasets
raw_transcript_rnaseq_data <- utils::read.table(base::file.path(config$rnaseq_results_dir,
                                                                "star_salmon/salmon.merged.transcript_counts.tsv"),
                                                header = TRUE,
                                                stringsAsFactors = FALSE,
                                                check.names = FALSE) %>%
  # rename columns
  dplyr::rename(gene_transcript = tx) %>%
  # remove columns before normalization step
  dplyr::select(-gene_id) %>%
  # convert column to rowname so the data can be normalised later
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var="gene_transcript")

raw_gene_rnaseq_data <- utils::read.table(base::file.path(config$rnaseq_results_dir,
                                                          "star_salmon/salmon.merged.gene_counts_length_scaled.tsv"),
                                          header = TRUE,
                                          stringsAsFactors = FALSE,
                                          check.names = FALSE) %>%
  # rename columns
  dplyr::rename(gene_transcript = gene_id) %>%
  # remove columns before normalization step
  dplyr::select(-gene_name) %>%
  # convert column to rowname so the data can be normalised later
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var="gene_transcript")

# remove rows that are all 0 (no counts for that given gene/transcript found in any sample)
raw_transcript_rnaseq_data <- raw_transcript_rnaseq_data[rowSums(raw_transcript_rnaseq_data[])>0,]
raw_gene_rnaseq_data <- raw_gene_rnaseq_data[rowSums(raw_gene_rnaseq_data[])>0,]

# create a vector defining the count datasets to analyse (that are set to TRUE) based on the yaml user configuration file
to_analyse <- config[c("transcript_rnaseq",
                       "gene_rnaseq")] %>%
  base::as.data.frame() %>%
  base::t() %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("dataset") %>%
  dplyr::rename("analyse" = "V1") %>%
  tidyr::separate(dataset, c("dataset", "pipeline")) %>%
  dplyr::filter(analyse == TRUE) %>%
  dplyr::mutate(dataset_name = base::paste0(dataset, "_", pipeline)) %>%
  dplyr::pull(dataset_name)

# set up a vector of all the possible datasets to analyse
count_datasets <- base::list(transcript_rnaseq = raw_transcript_rnaseq_data,
                             gene_rnaseq = raw_gene_rnaseq_data)

# filter all possible datasets by the datasets the user has specified to analyse
count_datasets <- count_datasets[to_analyse]

# evalutate/setup minimum logFC threshold
min_logfc <- base::eval(base::parse(text = config$min_lfc))

# specify treatments by creating a string of conditions that match the order of the columns/samples in the count data
# get the treatments and samples names from the metadata file
treatments <- metadata %>%
  dplyr::select(sample, treatment)

# also sort by the sample column (important so it matches the order of the samples count datasets)
# this is critical for DESeq2 - it assumes they are in the same order
treatments <- treatments[gtools::mixedorder(base::as.character(treatments$sample)),]

# extract only the conditions/groups and create a list out of it
ordered_treatments <- treatments %>%
  dplyr::pull(treatment)

# make sure the samples (columns) in all the datasets are in the same order for specifying the treatments downstream
# (that depends on the columns being in the correct order)  (loop over all count datasets)
count_datasets <- base::lapply(count_datasets, function(x) {
  
  x[ , gtools::mixedsort(names(x))]
  
})

# convert datasets to matrices (loop over all count datasets)
# also round read counts to the nearest integer to avoid a downstream error with DESeqDataSetFromMatrix() (see this discussion https://www.biostars.org/p/368158/)
count_datasets <- base::lapply(count_datasets, function(x) {
  
  x %>%
    base::round() %>%
    base::as.matrix(sep = "\t", row.names = "gene_transcript_id")
  
})

# create DGEList objects for the datasets (loop over all count datasets)
dge <- base::lapply(count_datasets, function(x) {
  
  edgeR::DGEList(counts = x, group = ordered_treatments)
  
})

# filter lowly expressed RNA
keep <- base::lapply(dge, function(x) {
  
  edgeR::filterByExpr(x, group=ordered_treatments, min.count = config$min_count, min.total.count = config$min_total_count)
  
})

dge_transcript_rnaseq_filtered <- dge$transcript_rnaseq[keep$transcript_rnaseq,, keep.lib.sizes=FALSE]
dge_gene_rnaseq_filtered <- dge$gene_rnaseq[keep$gene_rnaseq,, keep.lib.sizes=FALSE]

# establish if there are 2 or more rows of data in the count datasets
gene_rnaseq_enough_rows <- nrow(dge_gene_rnaseq_filtered$counts) >=2
transcript_rnaseq_enough_rows <- nrow(dge_transcript_rnaseq_filtered$counts) >=2
