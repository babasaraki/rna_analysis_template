# load libraries
library(shiny)
library(shinythemes)
library(shinydashboard)
library(dplyr)
library(plotly)
library(shinycssloaders)
library(DT)
library(crosstalk)
library(gtools)

# define UI for app ----
ui <- shiny::fluidPage(
  
  # app theme ----
  theme = shinythemes::shinytheme("yeti"),
  
  # define page layout ----
  shiny::verticalLayout(
    
    # render user options panel ----
    shiny::wellPanel(
      
      # user input to choose the dimensions to plot ----
      shiny::selectizeInput("dimension_choice",
                            h4("Select the dimensions to plot:"),
                            choices = c("1 and 2"),
                            selected = "1 and 2",
                            options = list(maxItems = 1)),
      
      # add help text that reminds the user that this pca data is based on all datasets
      shiny::helpText(tags$i(
        h5("Note. these plots explore all genes/transcripts identified in the samples. Box/lasso select data in the plots to view those data points in the table. Alternatively, click to select rows in the table to show that value in the plot. Hold the select key when selecting rows in tables or points in graphs for multiselect mode. Double click on a plot (in 'Box Select' or 'Lasso Select' mode) to reset selection. Zoom in on points in the PCA using the 'zoom' mode"))),
      
    ),
    
    shiny::wellPanel(style = "padding: 80px;",
                     
                     crosstalk::bscols(plotly::plotlyOutput(outputId = "scree_plot"), DT::dataTableOutput("table_scree"), widths = c(5,6)) %>%
                       shinycssloaders::withSpinner(color="#0097db")
    ),
    
    shiny::wellPanel(style = "padding: 80px;",
                     
                     crosstalk::bscols(plotly::plotlyOutput(outputId = "variables_plot"), DT::dataTableOutput("table_variables"), widths = c(5,7)) %>%
                       shinycssloaders::withSpinner(color="#0097db")
                     
    ),
    
    shiny::wellPanel(style = "padding: 80px;",
                     
                     # add help text that reminds the user that this pca data is based on all datasets
                     helpText(tags$i(
                       h5("Please be patient, I'm a little slow because there are many data points!"))),
                     
                     crosstalk::bscols(plotly::plotlyOutput(outputId = "individuals_plot"), DT::dataTableOutput("table_individuals"), widths = c(5,7)) %>%
                       shinycssloaders::withSpinner(color="#0097db")
    )
    
  )
)
