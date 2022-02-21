# load libraries
library(dplyr)
library(textshape)
library(edgeR)
library(tibble)
library(gtools)
library(janitor)

# read in yaml config file
config <- yaml::yaml.load_file("./config/config.yaml")

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

# calculate counts per million (loop over all the count datasets the user has specified)
counts_cpm <- base::lapply(base::seq_along(count_datasets), function(y, pipeline, dataset, i){
  
  # calculate counts per million for all datasets in the list
  edgeR::cpm(y[[i]]) %>%
    base::as.data.frame() %>%
    # sort all the columns in the different datasets so they are consistent among each other
    # (important for specifying the treatment groups downstream)
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    tibble::rownames_to_column("gene_transcript") %>%
    # make data long
    tidyr::pivot_longer(-gene_transcript, names_to = "sample", values_to = "counts_per_million") %>%
    # create a column that has defines the dataset the data has come from
    # based on the names of the "count_dataset" list, grab everything BEFORE 
    # the underscore
    dplyr::mutate(dataset = dataset[[i]]) %>%
    # also create a column that defines the pipeline the data has come from
    # based on the names of the "count_dataset" list, grab everything AFTER
    # the underscore
    dplyr::mutate(pipeline = pipeline[[i]])
}, y=count_datasets,
dataset=base::sub("\\_.*", "", base::paste(base::names(count_datasets))),
pipeline=base::sub(".*\\_", "", base::paste(base::names(count_datasets)))
)

# calculate log counts per million (loop over all the count datasets the user has specified)
counts_lcpm <- base::lapply(seq_along(count_datasets), function(y, pipeline, dataset, i){
  
  # calculate log counts per million for all datasets in the list
  edgeR::cpm(y[[i]], log = TRUE) %>%
    base::as.data.frame() %>%
    # sort all the columns in the different datasets so they are consistent among each other
    # (important for specifying the treatment groups downstream)
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    tibble::rownames_to_column("gene_transcript") %>%
    # make data long
    tidyr::pivot_longer(-gene_transcript, names_to = "sample", values_to = "log_counts_per_million") %>%
    # create a column that has defines the dataset the data has come from
    # based on the names of the "count_dataset" list, grab everything BEFORE 
    # the underscore
    dplyr::mutate(dataset = dataset[[i]]) %>%
    # also create a column that defines the pipeline the data has come from
    # based on the names of the "count_dataset" list, grab everything AFTER
    # the underscore
    dplyr::mutate(pipeline = pipeline[[i]])
}, y=count_datasets,
dataset=base::sub("\\_.*", "", base::paste(base::names(count_datasets))),
pipeline=base::sub(".*\\_", "", base::paste(base::names(count_datasets)))
)

# prepare raw counts (loop over all the count datasets the user has specified)
counts_raw <- base::lapply(seq_along(count_datasets), function(y, pipeline, dataset, i){
  
  # prepare raw counts for all datasets in the list
  y[[i]] %>%
    base::as.data.frame() %>%
    # sort all the columns in the different datasets so they are consistent among each other
    # (important for specifying the treatment groups downstream)
    dplyr::select(gtools::mixedsort(tidyselect::peek_vars())) %>%
    tibble::rownames_to_column("gene_transcript") %>%
    # make data long
    tidyr::pivot_longer(-gene_transcript, names_to = "sample", values_to = "raw_counts") %>%
    # create a column that has defines the dataset the data has come from
    # based on the names of the "count_dataset" list, grab everything BEFORE 
    # the underscore
    dplyr::mutate(dataset = dataset[[i]]) %>%
    # also create a column that defines the pipeline the data has come from
    # based on the names of the "count_dataset" list, grab everything AFTER
    # the underscore
    dplyr::mutate(pipeline = pipeline[[i]])
}, y=count_datasets,
dataset=base::sub("\\_.*", "", base::paste(base::names(count_datasets))),
pipeline=base::sub(".*\\_", "", base::paste(base::names(count_datasets)))
)

# collapse my list of dataframes into single dataframes
counts_cpm <- base::Reduce (rbind, counts_cpm)
counts_lcpm <- base::Reduce(rbind, counts_lcpm)
counts_raw <- base::Reduce(rbind, counts_raw)

# join all three dataframes into one large dataframe of all count data!
counts <- dplyr::full_join(counts_cpm, counts_lcpm, by = c("gene_transcript", "sample", "pipeline", "dataset")) 
counts <- dplyr::full_join(counts, counts_raw, by = c("gene_transcript", "sample", "pipeline", "dataset"))

# fix for css highlighting in expression plotting shiny app not working for genes/transcripts/rows with ":", "|" or "."
# doing this here so the genes/transcripts are named consistently throughout all documents
# | needed to be escaped with \\ in order to be interpreted correctly
# also need to use gsub instead of sub to replace all occurrences instead of just the first one
counts <- counts%>%
  dplyr::mutate(gene_transcript = base::gsub(":|\\||\\.", "_", gene_transcript))

# write the data to a csv file so I can use it in other documents
utils::write.csv(counts, "./prepare_counts/counts.csv", row.names = FALSE)

# also save some rds objects
dir.create("./prepare_counts/rds_objects/", showWarnings = FALSE)

# natural sorting of columns/sample names
raw_transcript_rnaseq_data <- raw_transcript_rnaseq_data %>%
  dplyr::select(gtools::mixedsort(tidyselect::peek_vars()))

raw_gene_rnaseq_data <- raw_gene_rnaseq_data %>%
  dplyr::select(gtools::mixedsort(tidyselect::peek_vars()))

# also apply changes to gene/transcript names here
rownames(raw_transcript_rnaseq_data) <- base::gsub(":|\\||\\.", "_", rownames(raw_transcript_rnaseq_data))
rownames(raw_gene_rnaseq_data) <- base::gsub(":|\\||\\.", "_", rownames(raw_gene_rnaseq_data))

# calculate counts per million
cpm_transcript_rnaseq_data <- raw_transcript_rnaseq_data %>%
  edgeR::cpm()
cpm_gene_rnaseq_data <- raw_gene_rnaseq_data %>%
  edgeR::cpm()

# calculate log counts per million
lcpm_transcript_rnaseq_data <- raw_transcript_rnaseq_data %>%
  edgeR::cpm(log = TRUE)
lcpm_gene_rnaseq_data <- raw_gene_rnaseq_data %>%
  edgeR::cpm(log = TRUE)

# save all rds objects to file
base::saveRDS(raw_transcript_rnaseq_data, file = "./prepare_counts/rds_objects/raw_transcript_rnaseq_counts.rds")
base::saveRDS(raw_gene_rnaseq_data, file = "./prepare_counts/rds_objects/raw_gene_rnaseq_counts.rds")
base::saveRDS(cpm_transcript_rnaseq_data, file = "./prepare_counts/rds_objects/cpm_transcript_rnaseq_counts.rds")
base::saveRDS(cpm_gene_rnaseq_data, file = "./prepare_counts/rds_objects/cpm_gene_rnaseq_counts.rds")
base::saveRDS(lcpm_transcript_rnaseq_data, file = "./prepare_counts/rds_objects/lcpm_transcript_rnaseq_counts.rds")
base::saveRDS(lcpm_gene_rnaseq_data, file = "./prepare_counts/rds_objects/lcpm_gene_rnaseq_counts.rds")

# clean up
rm(list = ls())
