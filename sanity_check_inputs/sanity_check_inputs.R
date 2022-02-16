##################################################################
##                          Setup                               ##
##################################################################

# load libraries
library(dplyr)

# read in yaml config file
config <- yaml::yaml.load_file("./config/config.yaml")

# read in metadata
metadata <- utils::read.csv(file.path(config$metadata))

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

# create an empty list to put all the results of the issues to later return to user
issues <- base::list()

##################################################################
##                   Check config file                          ##
##################################################################

# check for missing values
missing_values_in_config <- base::lapply(config, is.null)

# get the values in the list (corresponding to values in the config file) that are missing
missing_values_in_config <- base::which(base::unlist(missing_values_in_config))

# if there is a missing value in the config file, create a message to append to my list of issues
if (base::length(missing_values_in_config >= 1)) {
  
  issues <- c(issues, base::paste0("The following values are missing in the configuaration file: '", base::names(missing_values_in_config)))
  
}

# check the values that should be read into R as TRUE/FALSE (logical) have been entered correctly in the config file
# if they haven't been entered correctly, create a message to append to my list of issues
if (base::is.logical(config$gene_rnaseq) == FALSE) {
  
  issues <- c(issues, base::paste0("'gene_rnaseq' is incorrectly formatted in the configuaration file. Value needs to be logical, either TRUE or FALSE (and not surrounded in quotation marks)."))
  
}

if (base::is.logical(config$transcript_rnaseq) == FALSE) {
  
  issues <- c(issues, base::paste0("'transcript_rnaseq' is incorrectly formatted in the configuaration file. Value needs to be logical, either TRUE or FALSE (and not surrounded in quotation marks)."))
  
}

# check the values that should be read into R as integers have been entered correctly in the config file
# if they haven't been entered correctly, create a message to append to my list of issues
if (base::is.integer(config$min_count) == FALSE) {
  
  issues <- c(issues, base::paste0("'min_count' is incorrectly formatted in the configuaration file. Value needs to be an integer (and not surrounded in quotation marks)."))
  
}

if (base::is.integer(config$min_total_count) == FALSE) {
  
  issues <- c(issues, base::paste0("'min_total_count' is incorrectly formatted in the configuaration file. Value needs to be an integer (and not surrounded in quotation marks)."))
  
}

# check the values that should be read into R as characters have been entered correctly in the config file
# if they haven't been entered correctly, create a message to append to my list of issues
if (base::is.character(config$template_dir) == FALSE) {
  
  issues <- c(issues, base::paste0("'template_dir' is incorrectly formatted in the configuaration file. Value needs to be a character string (either surrounded in or not surrounded in quotation marks)."))
  
}

if (base::is.character(config$metadata_path) == FALSE) {
  
  issues <- c(issues, base::paste0("'metadata_path' is incorrectly formatted in the configuaration file. Value needs to be a character string (either surrounded in or not surrounded in quotation marks)."))
  
}

if (base::is.character(config$contrasts) == FALSE) {
  
  issues <- c(issues, base::paste0("'contrasts' is incorrectly formatted in the configuaration file. Value needs to be a character string (either surrounded in or not surrounded in quotation marks)."))
  
}

# check the file paths the user have passed exists
if (file.exists(config$fastq_dir) == FALSE) {
  
  issues <- c(issues, base::paste0("Directory passed to 'fastq_dir' in the configuration file doesn't exist."))
  
}

if (file.exists(config$template_dir) == FALSE) {
  
  issues <- c(issues, base::paste0("Directory passed to 'template_dir' in the configuration file doesn't exist."))

}

if (file.exists(config$metadata_path) == FALSE) {
  
  issues <- c(issues, base::paste0("File path passed to 'metadata_path' in the configuration file doesn't exist."))
  
}

if (file.exists(config$rnaseq_results_dir) == FALSE) {
  
  issues <- c(issues, base::paste0("Directory passed to 'rnaseq_results_dir' in the configuration file doesn't exist."))
  
}

##################################################################
##                  Check metadata file                         ##
##################################################################

# check the required column (sample) is present in the metadata
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)
if ((file.exists(config$metadata_path) == TRUE) & ("sample" %in% base::colnames(metadata)) == FALSE) {
  
  issues <- c(issues, base::paste0("No column named 'sample' in the metadata file at ",
                                   config$metadata_path,
                                   ". This column is required and must be all lower case."))
  
}

# check the required column (treatment) is present in the metadata
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)
if ((file.exists(config$metadata_path) == TRUE) & ("treatment" %in% base::colnames(metadata)) == FALSE) {
  
  issues <- c(issues, base::paste0("No column named 'treatment' in the metadata file at ",
                                   config$metadata_path,
                                   ". This column is required and must be all lower case."))
  
}

# check that all values of the 'sample' column in the metadata are unique
# ie. each row in the metadata represents a unique sample
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)
if ((file.exists(config$metadata_path) == TRUE) & (length(unique(metadata$sample)) != length(metadata$sample))) {
  
  issues <- c(issues, base::paste0("There appears to be duplicate samples in the 'sample' column in the metadata file at  ",
                                   config$metadata_path,
                                   ". Please ensure all sample ID's (values in the 'sample' column) are unique."))
  
}

# get the names of the fastq files in the fastq directory
fastq_file_names <- base::list.files(config$fastq_dir)

# strip .fastq.gz suffix to get just sample names
fastq_file_sample_names <- base::lapply(fastq_file_names, function(x)
  
  base::gsub(".fastq.gz", "", x)
  
)

# check the sample names of the fastq files correspond to the sample names in the metadata file
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)
if ((file.exists(config$metadata_path) == TRUE) & (any(fastq_file_sample_names %in% metadata$sample == FALSE) == TRUE)) {
  
  issues <- c(issues, base::paste0("Not all fastq file names in the fastq directory at ",
                                   config$fastq_dir,
                                   " correspond to the samples defined in the 'sample' column of the metadata file at ",
                                   config$metadata_path,
                                   ". Please ensure every sample in the fastq directory is defined in the metadata file."))
  
}

if ((file.exists(config$metadata_path) == TRUE) & (any(metadata$sample %in% fastq_file_sample_names == FALSE) == TRUE)) {
  
  issues <- c(issues, base::paste0("Not all the samples defined in the 'sample' column of the metadata file at ",
                                   config$metadata_path,
                                   " correspond to the fastq file names in the fastq directory at ",
                                   config$fastq_dir,
                                   ". Please ensure every sample defined in the metadata file correspond to the filenames of the files in the fastq directory."))
  
}

# check there are two or more samples in each treatment group (required for differential expression analysis)
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)

# get the number of samples per treatment group (based on the metadata file)
samples_per_treatment <- metadata %>%
  dplyr::group_by(treatment) %>%
  dplyr::summarise(n_samples = n()) %>%
  dplyr::pull(n_samples)

# append message if requirement not satisfied
if ((file.exists(config$metadata_path) == TRUE) & (all(samples_per_treatment >= 2) == FALSE)) {
  
  issues <- c(issues, base::paste0("There appears to be treatment groups with less than 2 samples per treatment (no replicates) based on the metadata file at ",
                                   config$metadata_path,
                                   " . Please ensure each treatment group contains 2 or more samples."))
}

##################################################################
##                 Check raw count data                         ##
##################################################################

# check the sample names in the raw count data correspond to the sample names in the metadata file
# (don't print message if the metadata file path isn't correct since this error is already accounted for above)
if ((file.exists(config$metadata_path) == TRUE) & (any(base::colnames(raw_gene_rnaseq_data) %in% metadata$sample == FALSE) == TRUE)) {
  
  issues <- c(issues, base::paste0("Not all samples in the raw gene count data at ",
                                   base::file.path(config$rnaseq_results_dir, "star_salmon/salmon.merged.gene_counts_length_scaled.tsv"),
                                   " correspond to the samples defined in the 'sample' column of the metadata file at ",
                                   config$metadata_path,
                                   ". Please ensure every sample in the raw gene count data is defined in the metadata file."))
  
}

if ((file.exists(config$metadata_path) == TRUE) & (any(base::colnames(raw_transcript_rnaseq_data) %in% metadata$sample == FALSE) == TRUE)) {
  
  issues <- c(issues, base::paste0("Not all samples in the raw transcript count data at ",
                                   base::file.path(config$rnaseq_results_dir, "star_salmon/salmon.merged.transcript_counts.tsv"),
                                   " correspond to the samples defined in the 'sample' column of the metadata file at ",
                                   config$metadata_path,
                                   ". Please ensure every sample in the raw transcript count data is defined in the metadata file."))
  
}

# check that there is more than one row of data in the count matrices after filtering
# (this is required by the differential expression analysis)

# run a little script that filters the data
base::source("./sanity_check_inputs/test_filtering.R")

if ((file.exists(config$metadata_path) == TRUE) & (gene_rnaseq_enough_rows == FALSE)) {
  
  issues <- c(issues, base::paste0("There are not enough rows of data (genes) available after filtering lowly expressed genes to carry out a differential expression analysis of the gene_rnaseq dataset. Please set 'gene_rnaseq' in the metadata file at ",
                                   config$metadata_path,
                                   " to FALSE to not analyse this dataset."))
  
}

if ((file.exists(config$metadata_path) == TRUE) & (transcript_rnaseq_enough_rows == FALSE)) {
  
  issues <- c(issues, base::paste0("There are not enough rows of data (transcripts) available after filtering lowly expressed transcripts to carry out a differential expression analysis of the transcript_rnaseq dataset. Please set 'transcript_rnaseq' in the metadata file at ",
                                   config$metadata_path,
                                   " to FALSE to not analyse this dataset."))
  
}

##################################################################
##                  Send messages to user                       ##
##################################################################

# return the list of issues identified
# send message when no issues found or list of issues identified
if (base::length(issues) == 0) {
  
  message_for_user <- "No issues found in the input configuration/metadata/raw count data files! Woohoo!"
  
}

if (base::length(issues) >= 1) {
  
  message_for_user <- issues
  
}

# return message to user
print(message_for_user)


# clean up
rm(list = ls())
