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

# load database data
db <- DBI::dbConnect(RSQLite::SQLite(), "./expr_plotting_results/master-count.sqlite")
counts <- dplyr::tbl(db, "counts")
my_dataset <- dplyr::tbl(db, "my_dataset")
gene_transcript_dataset_pairs <- dplyr::tbl(db, "gene_transcript_choices")

# load differential expression data
diff_expr_data <- utils::read.csv("./expr_plotting_results/master_diff_expr_data.csv", header = TRUE, stringsAsFactors = FALSE)
sig_diff_expr_data_1 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_1.csv", header = TRUE, stringsAsFactors = FALSE)
sig_diff_expr_data_5 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_5.csv", header = TRUE, stringsAsFactors = FALSE)
sig_diff_expr_data_10 <- utils::read.csv("./expr_plotting_results/master_sig_diff_expr_data_10.csv", header = TRUE, stringsAsFactors = FALSE)

# read in yaml config file
config <- yaml::yaml.load_file("./expr_plotting_results/config.yaml")

# read in metadata
metadata <- utils::read.csv("./expr_plotting_results/metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# get list of datasets
dataset_choices <- my_dataset %>%
  dplyr::distinct() %>%
  dplyr::pull()

# get list of variables of possible interest to plot by
variables_of_interest <- metadata %>%
  dplyr::select(-sample) %>%
  base::colnames()

# split the metadata columns into a list
# do for all the columns/variables in the metadata that will be used (everything except the sample column)
metadata_split <- metadata %>%
  select(-sample) %>%
  as.list()

# for each value in list (originally metadata columns)...
widths <- lapply(metadata_split, function(x)
  
  # ...get number of observations of unique values...
  base::table(x) %>%
    base::as.data.frame() %>%
    # ...and divide this by the length of the list
    dplyr::pull(Freq)/base::length(x)
  
)

# define server logic for app
server <- function(input, output, session) {
  
  observe({
    
    # get the current datasets the user has chosen
    shiny::updateSelectizeInput(session, "dataset_choice", choices = dataset_choices, server = TRUE)
    
  })
  
  observe({
    
    # create a user input switch that will let the user view the data based on their main variable of interest
    shiny::updateSelectizeInput(session, "main_variable", choices = variables_of_interest, selected = "treatment", server = TRUE)
    
  })
  
  # return genes/transcripts available for the user to choose from (depends on the users choice of datasets)
  gene_transcript_choices <- shiny::reactive({
    
    if(base::is.null(input$dataset_choice)){
      return(gene_transcript_dataset_pairs) %>%
        dplyr::distinct() %>%
        dplyr::pull()
    }
    
    gene_transcript_dataset_pairs %>%
      dplyr::filter(dataset==local(input$dataset_choice)) %>%
      dplyr::distinct() %>%
      dplyr::pull()
  })
  
  
  # get the current gene_transcript the user has chosen
  observe({
    
    shiny::updateSelectizeInput(session, "gene_transcript_choice", choices = gene_transcript_choices(), server=TRUE)
    
  })
  
  # subset count data based on the user choice of dataset and gene/transcript
  # the dplyr::collect() function is required to force evaluation of "input$gene_transcript_choice" and "input$dataset_choice" - otherwise is does lazy evaluation
  # the !! syntax is required to force dplyr to evaluate the user inputs (eg. input$dataset_choice) before operating on it (eg. dplyr::filter())
  subset_count_data <- shiny::reactive({
    
    shiny::req(input$dataset_choice, input$gene_transcript_choice)
    
    counts %>%
      dplyr::filter(dataset==!!input$dataset_choice & gene_transcript %in% !!input$gene_transcript_choice) %>%
      dplyr::collect() %>%
      # make variables factors so they plot correctly (all columns except several aforementioned columns)
      dplyr::mutate(across(!gene_transcript &
                             !dataset &
                             !sample &
                             !raw_counts &
                             !counts_per_million, factor))
    
  })
  
  # subset differential expression data based on the user choice of dataset and gene/transcript
  # the dplyr::collect() function is required to force evaluation of "input$gene_transcript_choice" and "input$dataset_choice" - otherwise is does lazy evaluation
  # the !! syntax is required to force dplyr to evaluate the user inputs (eg. input$dataset_choice) before operating on it (eg. dplyr::filter())
  subset_diff_expr_data <- shiny::reactive({
    
    shiny::req(input$dataset_choice, input$gene_transcript_choice)
    
    diff_expr_data %>% 
      dplyr::filter(dataset==!!input$dataset_choice & gene_transcript %in% !!input$gene_transcript_choice) %>%
      dplyr::collect()
  })
  
  # setup a reactive function that grabs the plot widths based on the current main variable the user has selected to view
  subplot_width <- shiny::reactive({
    
    widths[[input$main_variable]]
    
  })
  
  # create interactive table of the users current selection of gene/transcript and dataset that show differential expressions results
  output$table_diff_expr <- DT::renderDataTable({
    
    DT::datatable(subset_diff_expr_data() %>%
                    dplyr::select(gene_transcript,
                                  dataset,
                                  comparison,
                                  pipeline,
                                  diff_expr_method,
                                  log_fc,
                                  p_value,
                                  adj_p_value,
                                  significance) %>%
                    dplyr::mutate(across(c(dataset,
                                           pipeline,
                                           diff_expr_method,
                                           comparison,
                                           significance), base::as.factor)) %>%
                    dplyr::mutate(across(c(log_fc,
                                           p_value,
                                           adj_p_value), base::as.double)) %>%
                    dplyr::mutate(across(c(p_value,
                                           adj_p_value), ~base::round(.x, digits = 6))) %>%
                    dplyr::mutate(across(c(log_fc), ~base::round(.x, digits = 4))),
                  filter = "top",
                  rownames = FALSE,
                  colnames = c("Gene/transcript",
                               "Dataset",
                               "Treatment comparison",
                               "Pipeline",
                               "Differential expression method",
                               "Log fold change",
                               "p-value",
                               "Adjusted p-value",
                               "Significance"),
                  extensions = base::list("ColReorder" = NULL,
                                          "Buttons" = NULL,
                                          "FixedColumns" = base::list(leftColumns=1)),
                  options = base::list(
                    dom = "BRrltpi",
                    autoWidth = TRUE,
                    columnDefs = base::list(base::list(width = "200px", targets = 0)),
                    lengthMenu = base::list(c(10, 50, -1), c("10", "50", "All")),
                    ColReorder = TRUE,
                    buttons =
                      base::list("copy",
                                 base::list(extend = "collection",
                                            buttons = c("csv", "excel", "pdf"),
                                            text = "Download"),
                                 I("colvis")))) %>%
      # Highlight the genes/transcripts/rows that are significantly differentially expressed
      # A different highlight is used for the three significance levels (1%, 5%, 10%)
      # This matches the css highlighting in the genes/transcripts drop down box
    DT::formatStyle("significance",
                    target = "row",
                    backgroundColor = styleEqual(c("significant_1%",
                                                   "significant_5%",
                                                   "significant_10%"),
                                                 c("#85C659",
                                                   "#febf2a",
                                                   "#ec1515")))
    
  })
  
  # generate boxplot of normalised counts per million averaged over samples in the levels/groups of the main variable of interest the user chooses to view
  output$box_plot_cpm <- plotly::renderPlotly({
    
    # print a message if the user hasn't yet selected a gene/transcript
    shiny::validate(need(input$gene_transcript_choice, "Select an gene/transcript from the drop down box above to create gene/transcript expression plots!"))
    
    # print a message if no count per million data could be calculated for this data
    shiny::validate(need(subset_count_data()$counts_per_million, "Count per million data could not be calculated for this gene/transcript - the counts were too small"))
    
    subset_count_data() %>%
      base::split(base::list(subset_count_data() %>% dplyr::pull(get(input$main_variable)), subset_count_data()$gene_transcript)) %>%
      base::lapply(function(x) {
        plotly::plot_ly(data = x,
                        x = ~interaction(get(input$main_variable), gene_transcript),
                        y = ~counts_per_million,
                        split = ~pipeline,
                        color = ~get(input$main_variable),
                        colors = c("#0097db", "#85C659", "#ec1515", "#febf2a", "#784f96"),
                        type = "box",
                        hoverinfo = "y") %>%
          plotly::layout(yaxis = base::list(title = "Counts per million"),
                         xaxis = base::list(title = "", tickangle = 270, type = "category")) }) %>%
      plotly::subplot(shareY = TRUE, widths = c(subplot_width()))
  })
  
  # generate scatterplot of normalised counts per million for all samples in the levels/groups of the main variable of interest the user chooses to view
  output$scatterplot_cpm_by_sample <- plotly::renderPlotly({
    
    # print a message if no count per million data could be calculated for this data
    shiny::validate(need(subset_count_data()$counts_per_million, "Count per million data could not be calculated for this gene/transcript"))
    
    subset_count_data() %>%
      base::split(base::list(subset_count_data() %>% dplyr::pull(get(input$main_variable)), subset_count_data()$gene_transcript)) %>%
      base::lapply(function(x) {
        plotly::plot_ly(data = x,
                        x = x$sample,
                        y = ~counts_per_million,
                        split = ~pipeline,
                        color = ~get(input$main_variable),
                        colors = c("#0097db", "#85C659", "#ec1515", "#febf2a", "#784f96"),
                        type = "scatter",
                        mode  = "markers",
                        marker = base::list(opacity = 0.5),
                        hoverinfo = "text",
                        text = ~paste("<b>", input$gene_transcript_choice,"</b>",
                                      "</br><i>", stringr::str_to_title(input$main_variable), ": ", get(input$main_variable), "</i>",
                                      "</br><br> Counts per million:", base::format(counts_per_million, big.mark = ",", scientific = FALSE, digits = 2),
                                      "</br> Sample:", x$sample,
                                      "</br> Pipeline:", pipeline)) %>%
          plotly::layout(yaxis = base::list(title = "Counts per million"),
                         xaxis = base::list(title = "", tickangle = 270, type = "category")) }) %>%
      plotly::subplot(shareY = TRUE, widths = c(subplot_width()))
  })
  
  # generate boxplot of raw counts for all samples in the levels/groups of the main variable of interest the user chooses to view
  output$box_plot_raw <- plotly::renderPlotly({
    
    # generate boxplot of raw counts per averaged over samples in the levels/groups of the main variable of interest the user chooses to view
    subset_count_data() %>%
      base::split(base::list(subset_count_data() %>% dplyr::pull(get(input$main_variable)), subset_count_data()$gene_transcript)) %>%
      base::lapply(function(x) {
        plotly::plot_ly(data = x,
                        x = ~interaction(get(input$main_variable), gene_transcript),
                        y = ~raw_counts,
                        split = ~pipeline,
                        color = ~get(input$main_variable),
                        colors = c("#0097db", "#85C659", "#ec1515", "#febf2a", "#784f96"),
                        type = "box",
                        hoverinfo = "y") %>%
          plotly::layout(yaxis = base::list(title = "Raw counts"),
                         xaxis = base::list(title = "", tickangle = 270, type = "category")) }) %>%
      plotly::subplot(shareY = TRUE, widths = c(subplot_width()))
  })
  
  # generate scatterplot of raw counts for all sample in the levels/groups of the main variable of interest the user chooses to view
  output$scatterplot_raw_by_sample <- plotly::renderPlotly({
    
    subset_count_data() %>%
      base::split(base::list(subset_count_data() %>% dplyr::pull(get(input$main_variable)), subset_count_data()$gene_transcript)) %>%
      base::lapply(function(x) {
        plotly::plot_ly(data = x,
                        x = x$sample,
                        y = ~raw_counts,
                        split = ~pipeline,
                        color = ~get(input$main_variable),
                        colors = c("#0097db", "#85C659", "#ec1515", "#febf2a", "#784f96"),
                        type = "scatter",
                        mode  = "markers",
                        marker = base::list(opacity = 0.5),
                        hoverinfo = "text",
                        text = ~paste("<b>", input$gene_transcript_choice,"</b>",
                                      "</br><i>", stringr::str_to_title(input$main_variable), ": ", get(input$main_variable), "</i>",
                                      "</br><br> Raw count:", base::format(raw_counts, big.mark = ",", scientific = FALSE, digits = 2),
                                      "</br> Sample:", x$sample,
                                      "</br> Pipeline:", pipeline)) %>%
          plotly::layout(yaxis = base::list(title = "Raw counts"),
                         xaxis = base::list(title = "", tickangle = 270, type = "category")) }) %>%
      plotly::subplot(shareY = TRUE, widths = c(subplot_width()))
  })
  
  # create interactive table of the users current selection of gene/transcript and dataset that includes gene_transcript counts etc
  output$table <- DT::renderDataTable({
    
    DT::datatable(subset_count_data() %>% dplyr::select(gene_transcript,
                                                        dataset,
                                                        sample,
                                                        pipeline,
                                                        raw_counts,
                                                        counts_per_million,
                                                        treatment) %>%
                    dplyr::mutate(across(c(dataset,
                                           sample,
                                           pipeline,
                                           treatment), base::as.factor)) %>%
                    dplyr::mutate(across(c(raw_counts,
                                           counts_per_million), base::as.integer)),
                  selection = "single",
                  filter = "top",
                  rownames = FALSE,
                  colnames = c("Gene/transcript",
                               "Dataset",
                               "Sample",
                               "Pipeline",
                               "Raw counts",
                               "Counts per million",
                               "Treatment"),
                  extensions = base::list("ColReorder" = NULL,
                                          "Buttons" = NULL,
                                          "FixedColumns" = base::list(leftColumns=1)),
                  options = base::list(
                    dom = "BRrltpi",
                    autoWidth = TRUE,
                    columnDefs = base::list(base::list(width = "200px", targets = 0)),
                    lengthMenu = base::list(c(10, 50, -1), c("10", "50", "All")),
                    ColReorder = TRUE,
                    buttons =
                      base::list("copy",
                                 base::list(extend = "collection",
                                            buttons = c("csv", "excel", "pdf"),
                                            text = "Download"),
                                 I("colvis"))))
  })
}
