##################################################################
##                          Setup                               ##
##################################################################

# load libraries
library(dplyr)
library(DBI)
library(yaml)

# read in processed count data (processed in the differential expression document)
counts <- utils::read.csv("./prepare_counts/counts.csv")

# read in yaml config file
config <- yaml::yaml.load_file("./config/config.yaml")

# read in metadata
metadata <- utils::read.csv(file.path(config$metadata))

# create output directory
dir.create("./expression_plotting/expr_plotting_results/", showWarnings = FALSE)

# subset all count data to the count data specified to be analysed by the user in the config file
# this will remove rows of data matching the said conditions

if(config$transcript_rnaseq == "FALSE") {
  
  counts <- counts %>%
    dplyr::filter(!((dataset == "transcript") & (pipeline == "rnaseq")))
  
}

if(config$gene_rnaseq == "FALSE") {
  
  counts <- counts %>%
    dplyr::filter(!((dataset == "gene") & (pipeline == "rnaseq")))
  
}

# read in metadata
metadata <- utils::read.csv(config$metadata_path)

# read in differential expression data
diff_expr_data <- utils::read.table("./diff_expression/diff_expr_results/diff_expr_results.tsv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

##################################################################
##                   Prep of count data                         ##
##################################################################

# append metadata to the count data
counts <- dplyr::left_join(counts, metadata, by = "sample")

##################################################################
##           Prep of differential expression data               ##
##################################################################

# extract a list of the significant genes/transcripts (at three significance levels)
sig_diff_expr_data_1 <- diff_expr_data %>%
  dplyr::filter(significance == "significant_1%") %>%
  dplyr::select(gene_transcript) %>%
  dplyr::distinct(gene_transcript)

sig_diff_expr_data_5 <- diff_expr_data %>%
  dplyr::filter(significance == "significant_5%") %>%
  dplyr::select(gene_transcript) %>%
  dplyr::distinct(gene_transcript)

sig_diff_expr_data_10 <- diff_expr_data %>%
  dplyr::filter(significance == "significant_10%") %>%
  dplyr::select(gene_transcript) %>%
  dplyr::distinct(gene_transcript)

#################################################################
##                        Write to file                        ##
#################################################################

# differential expression data
utils::write.csv(diff_expr_data, file = "./expression_plotting/expr_plotting_results/master_diff_expr_data.csv", row.names = FALSE)
utils::write.csv(sig_diff_expr_data_1, file = "./expression_plotting/expr_plotting_results/master_sig_diff_expr_data_1.csv", row.names = FALSE)
utils::write.csv(sig_diff_expr_data_5, file = "./expression_plotting/expr_plotting_results/master_sig_diff_expr_data_5.csv", row.names = FALSE)
utils::write.csv(sig_diff_expr_data_10, file = "./expression_plotting/expr_plotting_results/master_sig_diff_expr_data_10.csv", row.names = FALSE)

# write the config file needed for app to file within this project directory
# so it can be sent to shinyappsio and used
yaml::write_yaml(config, "./expression_plotting/expr_plotting_results/config.yaml")

# write the metadata file needed for app to file within this project directory
# so it can be sent to shinyappsio and used
utils::write.csv(metadata, file = "./expression_plotting/expr_plotting_results/metadata.csv", row.names = FALSE)

# count data
# this data was written to an sql lite database
# this significantly speeds up the app compared to reading a csv file of the data and filtering/subsetting the data

# initialise a database
db <- DBI::dbConnect(RSQLite::SQLite(), "./expression_plotting/expr_plotting_results/master-count.sqlite")

# get all unique datasets
dataset_choices <- base::data.frame(base::unique(base::as.character(counts$dataset)))

# create dataframe of unique combinations of datasets and gene_transcript
gene_transcript_choices <- base::data.frame(base::unique(counts[c("dataset", "gene_transcript")]))

# write dataframes and count data to database tables
DBI::dbWriteTable(db, "my_dataset", dataset_choices, overwrite=TRUE)
DBI::dbWriteTable(db, "gene_transcript_choices", gene_transcript_choices, overwrite=TRUE)
DBI::dbWriteTable(db, "counts", counts, overwrite=TRUE)

# create indexes for common data lookups the app with do (will further speed up app, particularly index_gene_transcript)
# note. I tried to make more indexes, but it made the database too large to deploy to shiny apps io
DBI::dbSendQuery(db,"CREATE INDEX index_gene_transcript ON counts (gene_transcript, dataset)")

# close the connection to the database
DBI::dbDisconnect(db)

# clean up
rm(list = ls())
