##################################################################
##                  Setup & calculate PCA                       ##
##################################################################

# load libraries
library(dplyr)
library(FactoMineR)
library(textshape)

# read in processed count data (processed in the differential expression document)
counts <- utils::read.csv("./prepare_counts/counts.csv")

# read in yaml config file
config <- yaml::yaml.load_file("./config/config.yaml")

# read in metadata
metadata <- utils::read.csv(config$metadata_path)

# create output directory
dir.create("./pca/pca_results/", showWarnings = FALSE)

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

# calculate PCA
pca <- counts %>%
  dplyr::filter(pipeline == "rnaseq") %>%
  dplyr::select(sample, gene_transcript, counts_per_million) %>%
  tidyr::pivot_wider(names_from = sample, values_from = counts_per_million) %>%
  textshape::column_to_rownames("gene_transcript") %>%
  base::as.matrix() %>%
  FactoMineR::PCA(graph = FALSE)

# extract the plotting coordinates from the PCA
scree_data <- pca$eig %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("component")

scree_data$component <- base::gsub("comp ", "", scree_data$component)

# extract the plotting coordinates from the PCA
variables_coords <- pca$var$coord %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::rename_with(function(x) base::paste0(x, "_coord"), where(is.numeric))

# extract the cos values from the PCA
variables_cos <- pca$var$cos2 %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::rename_with(function(x) base::paste0(x, "_cos"), where(is.numeric))

# join the two dataframes (plotting coordinates and cos values)
variables_data <- dplyr::left_join(variables_coords, variables_cos, by = "sample")

# append the metadata
variables_data <- dplyr::left_join(variables_data, metadata, by = "sample")

# extract the plotting coordinates from the PCA
individuals_coords <- pca$ind$coord %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("gene_transcript") %>%
  dplyr::rename_with(function(x) base::paste0(x, "_coord"), where(is.numeric))

# extract the cos values from the PCA
individuals_cos <- pca$ind$cos2 %>%
  base::as.data.frame() %>%
  tibble::rownames_to_column("gene_transcript") %>%
  dplyr::rename_with(function(x) base::paste0(x, "_cos"), where(is.numeric))

# join the two dataframes (plotting coordinates and cos values)
individuals_data <- dplyr::left_join(individuals_coords, individuals_cos, by = "gene_transcript")

##################################################################
##                Prepare PCA data for app                      ##
##################################################################

# make matrices for all combinations of dimensions so they are pre-calculated and make the shiny app quicker
# grab both the dimension and cos data for a given dimension pair and rename columns to generic names for the app to grab

# 1st dimension compared to everything
individuals_1_2 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.1"), contains("Dim.2"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "y")

individuals_1_3 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.1"), contains("Dim.3"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "y")

individuals_1_4 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.1"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

individuals_1_5 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.1"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 2nd dimension compared to everything
individuals_2_3 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.2"), contains("Dim.3"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "y")

individuals_2_4 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.2"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

individuals_2_5 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.2"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 3rd dimension compared to everything
individuals_3_4 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.3"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

individuals_3_5 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.3"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 4th dimension compared to everything
individuals_4_5 <- individuals_data %>%
  dplyr::select(c(gene_transcript, contains("Dim.4"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 1st dimension compared to everything
variables_1_2 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.1"), contains("Dim.2"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "y")

variables_1_3 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.1"), contains("Dim.3"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "y")

variables_1_4 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.1"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

variables_1_5 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.1"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.1", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 2nd dimension compared to everything
variables_2_3 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.2"), contains("Dim.3"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "y")

variables_2_4 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.2"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

variables_2_5 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.2"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.2", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 3rd dimension compared to everything
variables_3_4 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.3"), contains("Dim.4"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "y")

variables_3_5 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.3"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.3", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

# 4th dimension compared to everything
variables_4_5 <- variables_data %>%
  dplyr::select(c(sample, treatment, contains("Dim.4"), contains("Dim.5"))) %>%
  dplyr::rename_with(gsub, pattern = "Dim.4", replacement = "x") %>%
  dplyr::rename_with(gsub, pattern = "Dim.5", replacement = "y")

#################################################################
##                        Write to file                        ##
#################################################################

# scree data
utils::write.csv(scree_data, "./pca/pca_results/scree.csv", row.names = FALSE)

# individuals (genes/transcripts) data
utils::write.csv(individuals_1_2, "./pca/pca_results/individuals_1_2.csv", row.names = FALSE)
utils::write.csv(individuals_1_3, "./pca/pca_results/individuals_1_3.csv", row.names = FALSE)
utils::write.csv(individuals_1_4, "./pca/pca_results/individuals_1_4.csv", row.names = FALSE)
utils::write.csv(individuals_1_5, "./pca/pca_results/individuals_1_5.csv", row.names = FALSE)
utils::write.csv(individuals_2_3, "./pca/pca_results/individuals_2_3.csv", row.names = FALSE)
utils::write.csv(individuals_2_4, "./pca/pca_results/individuals_2_4.csv", row.names = FALSE)
utils::write.csv(individuals_2_5, "./pca/pca_results/individuals_2_5.csv", row.names = FALSE)
utils::write.csv(individuals_3_4, "./pca/pca_results/individuals_3_4.csv", row.names = FALSE)
utils::write.csv(individuals_3_5, "./pca/pca_results/individuals_3_5.csv", row.names = FALSE)
utils::write.csv(individuals_4_5, "./pca/pca_results/individuals_4_5.csv", row.names = FALSE)

# variables (samples) data
utils::write.csv(variables_1_2, "./pca/pca_results/variables_1_2.csv", row.names = FALSE)
utils::write.csv(variables_1_3, "./pca/pca_results/variables_1_3.csv", row.names = FALSE)
utils::write.csv(variables_1_4, "./pca/pca_results/variables_1_4.csv", row.names = FALSE)
utils::write.csv(variables_1_5, "./pca/pca_results/variables_1_5.csv", row.names = FALSE)
utils::write.csv(variables_2_3, "./pca/pca_results/variables_2_3.csv", row.names = FALSE)
utils::write.csv(variables_2_4, "./pca/pca_results/variables_2_4.csv", row.names = FALSE)
utils::write.csv(variables_2_5, "./pca/pca_results/variables_2_5.csv", row.names = FALSE)
utils::write.csv(variables_3_4, "./pca/pca_results/variables_3_4.csv", row.names = FALSE)
utils::write.csv(variables_3_5, "./pca/pca_results/variables_3_5.csv", row.names = FALSE)
utils::write.csv(variables_4_5, "./pca/pca_results/variables_4_5.csv", row.names = FALSE)

# clean up
rm(list = ls())
