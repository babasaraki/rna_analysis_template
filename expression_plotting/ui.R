# load libraries
library(shiny)
library(shinythemes)
library(shinydashboard)
library(dplyr)
library(plotly)
library(shinycssloaders)
library(DBI)
library(dbplyr)
library(DT)

# load significant differential expression data
sig_diff_expr_data_1 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_1.csv", header = TRUE, stringsAsFactors = FALSE)
sig_diff_expr_data_5 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_5.csv", header = TRUE, stringsAsFactors = FALSE)
sig_diff_expr_data_10 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_10.csv", header = TRUE, stringsAsFactors = FALSE)

# define UI for app
ui <- shiny::fluidPage(
  
  # app theme
  theme = shinythemes::shinytheme("yeti"),
  
  # css to highlight the genes/transcripts in the drop down box (defined by gene_transcript_choice) that were found to be significantly differentially expressed
  # a different highlight is used for the three significance levels (1%, 5%, 10%)
  tags$head(
    base::lapply(sig_diff_expr_data_1, function(x){
      tags$style(HTML(
        base::paste0(".option[data-value=", x, "], .item[data-value=", x, "] {background: #85C659 !important;}")
      ))
    })
  ),
  
  tags$head(
    base::lapply(sig_diff_expr_data_5, function(x){
      tags$style(HTML(
        base::paste0(".option[data-value=", x, "], .item[data-value=", x, "] {background: #febf2a !important;}")
      ))
    })
  ),
  
  tags$head(
    base::lapply(sig_diff_expr_data_10, function(x){
      tags$style(HTML(
        base::paste0(".option[data-value=", x, "], .item[data-value=", x, "] {background: #ec1515 !important;}")
      ))
    })
  ),
  
  # define page layout
  shiny::verticalLayout(
    
    # render user options panel
    shiny::wellPanel(
      
      # user input to choose the dataset of interest to plot
      shiny::selectizeInput("dataset_choice",
                            h4("Select/search for a dataset:"),
                            choices = NULL,
                            options = list(maxItems = 1)),
      
      # user input to choose gene/transcript of interest to plot (options dependent on the first input - dataset_choice)
      shiny::selectizeInput("gene_transcript_choice",
                            h4("Select/search for an gene/transcript:"),
                            choices = NULL,
                            options = list(maxItems = 1)),
      
      # user input to choose the main variable of interest to compare
      shiny::selectizeInput("main_variable",
                            h4("Select the main variable of interest to compare:"),
                            choices = NULL,
                            options = list(maxItems = 1)),
      
      # add help text that reminds the user that all this differential expression data is based on
      # on only one of the variables they can choose to explore
      shiny::helpText(tags$i(
        h5("Explore the differential expression results, raw counts and counts per million by gene/transcript. Genes/transcripts highlighted in green, yellow and red were found to be significantly differentially expressed by at least one differential expression analysis.")))
    ),
    
    # plot all the tables and boxplots        
    shiny::wellPanel(
      
      # plot differential expression table
      DT::dataTableOutput("table_diff_expr") %>%
        shinycssloaders::withSpinner(color="#0097db")
    ),
    
    shiny::wellPanel(
      
      plotly::plotlyOutput(outputId = "scatterplot_raw_by_sample") %>%
        shinycssloaders::withSpinner(color="#0097db")
    ), 
    
    shiny::wellPanel(
      
      plotly::plotlyOutput(outputId = "box_plot_raw") %>%
        shinycssloaders::withSpinner(color="#0097db")
    ),
    
    shiny::wellPanel(
      
      plotly::plotlyOutput(outputId = "scatterplot_cpm_by_sample") %>%
        shinycssloaders::withSpinner(color="#0097db")
    ),
    
    shiny::wellPanel(
      
      plotly::plotlyOutput(outputId = "box_plot_cpm") %>%
        shinycssloaders::withSpinner(color="#0097db")
    ),
    
    shiny::wellPanel(
      
      DT::dataTableOutput("table") %>%
        shinycssloaders::withSpinner(color="#0097db")
      
    )
  )
)

