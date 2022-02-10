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

# linked plotly and datatable based on this stackoverfow question and answer:
# https://stackoverflow.com/questions/54344375/how-do-i-use-crosstalk-in-shiny-to-filter-a-reactive-data-table

# load data ----

# scree plot data 
scree <- utils::read.csv("./pca_results/scree.csv", header = TRUE, stringsAsFactors = FALSE)

# individuals (genes/transcripts)
individuals_1_2 <- utils::read.csv("./pca_results/individuals_1_2.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_1_3 <- utils::read.csv("./pca_results/individuals_1_3.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_1_4 <- utils::read.csv("./pca_results/individuals_1_4.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_1_5 <- utils::read.csv("./pca_results/individuals_1_5.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_2_3 <- utils::read.csv("./pca_results/individuals_2_3.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_2_4 <- utils::read.csv("./pca_results/individuals_2_4.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_2_5 <- utils::read.csv("./pca_results/individuals_2_5.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_3_4 <- utils::read.csv("./pca_results/individuals_3_4.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_3_5 <- utils::read.csv("./pca_results/individuals_3_5.csv", header = TRUE, stringsAsFactors = FALSE)
individuals_4_5 <- utils::read.csv("./pca_results/individuals_4_5.csv", header = TRUE, stringsAsFactors = FALSE)

# variables (samples)
variables_1_2 <- utils::read.csv("./pca_results/variables_1_2.csv", header = TRUE, stringsAsFactors = FALSE)
variables_1_3 <- utils::read.csv("./pca_results/variables_1_3.csv", header = TRUE, stringsAsFactors = FALSE)
variables_1_4 <- utils::read.csv("./pca_results/variables_1_4.csv", header = TRUE, stringsAsFactors = FALSE)
variables_1_5 <- utils::read.csv("./pca_results/variables_1_5.csv", header = TRUE, stringsAsFactors = FALSE)
variables_2_3 <- utils::read.csv("./pca_results/variables_2_3.csv", header = TRUE, stringsAsFactors = FALSE)
variables_2_4 <- utils::read.csv("./pca_results/variables_2_4.csv", header = TRUE, stringsAsFactors = FALSE)
variables_2_5 <- utils::read.csv("./pca_results/variables_2_5.csv", header = TRUE, stringsAsFactors = FALSE)
variables_3_4 <- utils::read.csv("./pca_results/variables_3_4.csv", header = TRUE, stringsAsFactors = FALSE)
variables_3_5 <- utils::read.csv("./pca_results/variables_3_5.csv", header = TRUE, stringsAsFactors = FALSE)
variables_4_5 <- utils::read.csv("./pca_results/variables_4_5.csv", header = TRUE, stringsAsFactors = FALSE)

# make sure the dimensions are read as a character
scree$component <- base::as.character(scree$component)

# natural sorting of xaxis
scree$component <- base::factor(scree$component, levels = gtools::mixedsort(scree$component))

# define server logic for app ----
server <- function(input, output, session) {
  
  # creating a shared dataframe that links the plotly scatterplot and datatable table (variables/samples)
  linked_data_scree <- shiny::reactive({
    
    sd_scree <- SharedData$new(scree)
    
  })
  
  # create a user input switch that will select the pca individuals (genes/transcripts) data based on the users selection of dimensions to view ----
  # also create a shared dataframe that links the plotly scatterplot and datatable table (individuals/genes/transcripts)
  linked_data_individuals <- shiny::reactive({
    
    data_ind <- base::switch(input$dimension_choice,
                             "1 and 2" = individuals_1_2,
                             "1 and 3" = individuals_1_3,
                             "1 and 4" = individuals_1_4,
                             "1 and 5" = individuals_1_5,
                             "2 and 3" = individuals_2_3,
                             "2 and 4" = individuals_2_4,
                             "2 and 5" = individuals_2_5,
                             "3 and 4" = individuals_3_4,
                             "3 and 5" = individuals_3_5,
                             "4 and 5" = individuals_4_5)
    
    sd_individuals <- SharedData$new(data_ind)
    
  })
  
  # create a user input switch that will select the pca variables (samples) data based on the users selection of dimensions to view ----
  # also create a shared dataframe that links the plotly scatterplot and datatable table (variables/samples)
  linked_data_variables <- shiny::reactive({
    
    data_var <- base::switch(input$dimension_choice,
                             "1 and 2" = variables_1_2,
                             "1 and 3" = variables_1_3,
                             "1 and 4" = variables_1_4,
                             "1 and 5" = variables_1_5,
                             "2 and 3" = variables_2_3,
                             "2 and 4" = variables_2_4,
                             "2 and 5" = variables_2_5,
                             "3 and 4" = variables_3_4,
                             "3 and 5" = variables_3_5,
                             "4 and 5" = variables_4_5)
    
    sd_variables <- SharedData$new(data_var)
    
  })
  
  # create a switch that will output the dimensions we're viewing data based on the users selection of dimensions to view ----
  x_dim <- shiny::reactive({
    
    base::switch(input$dimension_choice,
                 "1 and 2" = "1",
                 "1 and 3" = "1",
                 "1 and 4" = "1",
                 "1 and 5" = "1",
                 "2 and 3" = "2",
                 "2 and 4" = "2",
                 "2 and 5" = "2",
                 "3 and 4" = "3",
                 "3 and 5" = "3",
                 "4 and 5" = "4")
  })
  
  y_dim <- shiny::reactive({
    
    base::switch(input$dimension_choice,
                 "1 and 2" = "2",
                 "1 and 3" = "3",
                 "1 and 4" = "4",
                 "1 and 5" = "5",
                 "2 and 3" = "3",
                 "2 and 4" = "4",
                 "2 and 5" = "5",
                 "3 and 4" = "4",
                 "3 and 5" = "5",
                 "4 and 5" = "5")
  })
  
  # grab variance explained by different dimensions based on the dimensions the user chooses to view ----
  dim_x_variance_explained <- shiny::reactive({
    
    scree %>%
      dplyr::filter(component == x_dim()) %>%
      dplyr::pull(percentage.of.variance) %>%
      base::format(big.mark = ",", scientific = FALSE, digits = 2)
    
  })
  
  dim_y_variance_explained <- shiny::reactive({
    
    scree %>%
      dplyr::filter(component == y_dim()) %>%
      dplyr::pull(percentage.of.variance) %>%
      base::format(big.mark = ",", scientific = FALSE, digits = 2)
    
  })
  
  # generate scree plot of the pca ----
  output$scree_plot <- plotly::renderPlotly({
    
    plotly::plot_ly(data = linked_data_scree(),
                    type = "bar",
                    x = ~component,
                    y = ~percentage.of.variance,
                    color = ~component,
                    colors = c("#0097db"),
                    textposition = "none",
                    height = "100%",
                    hoverinfo = "text",
                    text = ~paste("</br> Percentage of variance:", base::format(percentage.of.variance, big.mark = ",", scientific = FALSE, digits = 2), "%",
                                  "</br> Dimension:", component)) %>%
      plotly::layout(showlegend = FALSE,
                     barmode = "overlay",
                     xaxis = base::list(title = "Dimension"),
                     yaxis = base::list(title = "Percentage of variance", range = c(0, 100))) %>%
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("zoomIn2d",
                                                "zoomOut2d",
                                                "autoScale2d",
                                                "hoverClosestCartesian",
                                                "hoverCompareCartesian",
                                                "lasso2d")) %>%
      plotly::highlight(on = "plotly_selected",
                        opacityDim = 0.2,
                        selected = attrs_selected(opacity = 1))
    
  })
  
  # create interactive table of the scree data ----
  output$table_scree <- DT::renderDataTable({
    
    DT::datatable(linked_data_scree(),
                  filter = "top",
                  rownames = FALSE,
                  colnames = c("Dimension",
                               "Eigenvalue",
                               "Percentage of variance (%)",
                               "Cumulative percentage of variance (%)"),
                  extensions = base::list("ColReorder" = NULL,
                                          "Buttons" = NULL,
                                          "FixedColumns" = base::list(leftColumns=1)),
                  options = list(
                    dom = "BRrltpi",
                    autoWidth = TRUE,
                    buttons =
                      list("copy",
                           list(extend = "collection",
                                buttons = c("csv", "excel", "pdf"),
                                text = "Download"),
                           I("colvis")))) %>%
      plotly::highlight(on = "plotly_click",
                        opacityDim = 0.2,
                        selected = attrs_selected(opacity = 1))
    
  }, server = FALSE)
  
  # generate scatterplot of the pca (individuals - genes/transcripts) based on the dimensions the user chooses to view ----
  output$individuals_plot <- renderPlotly({
    
    p <- plotly::plot_ly(data = linked_data_individuals(),
                         x = ~x_coord,
                         y = ~y_coord,
                         color = ~x_cos,
                         colors = c("#febf2a", "#ec1515"),
                         height = "100%",
                         type = "scatter",
                         hoverinfo = "text",
                         text = ~paste("</br> Gene/transcript:", gene_transcript,
                                       "</br> Dim ", x_dim(), ":", base::format(x_coord, big.mark = ",", scientific = FALSE, digits = 0),
                                       "</br> Dim ", y_dim(), ":", base::format(y_coord, big.mark = ",", scientific = FALSE, digits = 0),
                                       "</br> Dim ", x_dim(), "cos2:", base::format(x_cos, big.mark = ",", scientific = FALSE, digits = 2))) %>%
      plotly::add_markers() %>%
      plotly::layout(xaxis = base::list(title = paste0("Dim ", x_dim(), " (", dim_x_variance_explained(), "% of variance)")),
                     yaxis = base::list(title = paste0("Dim ", y_dim(), " (", dim_y_variance_explained(), "% of variance)"))) %>%
      plotly::colorbar(title = paste0("Dim ", x_dim(), " cos2")) %>%
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("zoomIn2d",
                                                "zoomOut2d",
                                                "select2d",
                                                "autoScale2d",
                                                "hoverClosestCartesian",
                                                "hoverCompareCartesian")) %>%
      plotly::layout(dragmode = "lasso", showlegend = FALSE) %>%
      plotly::highlight(on = "plotly_selected",
                        opacityDim = 0.2,
                        selected = attrs_selected(opacity = 1))
    
    toWebGL(p)
    
  })
  
  # create interactive table of pca (individuals - genes/transcripts) based on the dimensions the user chooses to view ----
  output$table_individuals <- DT::renderDataTable({
    
    DT::datatable(linked_data_individuals(),
                  filter = "top",
                  rownames = FALSE,
                  colnames = c("Gene/transcript",
                               paste0("Dim ", x_dim(), " coord"),
                               paste0("Dim ", y_dim(), " coord"),
                               paste0("Dim ", x_dim(), " cos"),
                               paste0("Dim ", y_dim(), " cos")),
                  extensions = base::list("ColReorder" = NULL,
                                          "Buttons" = NULL,
                                          "FixedColumns" = base::list(leftColumns=1)),
                  options = list(
                    dom = "BRrltpi",
                    autoWidth = TRUE,
                    buttons =
                      list("copy",
                           list(extend = "collection",
                                buttons = c("csv", "excel", "pdf"),
                                text = "Download"),
                           I("colvis")))) %>%
      plotly::highlight(on = "plotly_selected",
                        opacityDim = 0.2,
                        selected = attrs_selected(opacity = 1))
    
  }, server = FALSE)
  
  # generate scatterplot of the pca (variables - samples) based on the dimensions the user chooses to view ----
  output$variables_plot <- plotly::renderPlotly({
    
    plotly::plot_ly(data = (linked_data_variables()),
                    x = ~x_coord,
                    y = ~y_coord,
                    color = ~x_cos,
                    colors = c("#febf2a", "#ec1515"),
                    type = "scatter",
                    height = "100%",
                    hoverinfo = "text",
                    text = ~paste("</br> Sample:", sample,
                                  "</br> Treatment:", treatment,
                                  "</br> Dim ", x_dim(), ":", base::format(x_coord, big.mark = ",", scientific = FALSE, digits = 2),
                                  "</br> Dim ", y_dim(), ":", base::format(y_coord, big.mark = ",", scientific = FALSE, digits = 2),
                                  "</br> Dim ", x_dim(), "cos2:", base::format(x_cos, big.mark = ",", scientific = FALSE, digits = 2))) %>%
      plotly::layout(xaxis = base::list(title = paste0("Dim ", x_dim(), " (", dim_x_variance_explained(), "% of variance)"), range = c(-1.2, 1.2)),
                     yaxis = base::list(title = paste0("Dim ", y_dim(), " (", dim_y_variance_explained(), "% of variance)"), range = c(-1.2, 1.2))) %>%
      plotly::colorbar(title = paste0("Dim ", x_dim(), " cos2")) %>%
      plotly::config(displaylogo = FALSE,
                     modeBarButtonsToRemove = c("zoomIn2d",
                                                "zoomOut2d",
                                                "select2d",
                                                "autoScale2d",
                                                "hoverClosestCartesian",
                                                "hoverCompareCartesian")) %>%
      plotly::layout(dragmode = "lasso", showlegend = FALSE) %>%
      plotly::highlight(on = "plotly_selected",
                        opacityDim = 0.1,
                        selected = attrs_selected(opacity = 1))
  })
  
  # create interactive table of pca (variables - genes/transcripts) based on the dimensions the user chooses to view ----
  output$table_variables <- DT::renderDataTable({
    
    DT::datatable(linked_data_variables(),
                  filter = "top",
                  rownames = FALSE,
                  colnames = c("Sample",
                               "Treatment",
                               paste0("Dim ", x_dim(), " coord"),
                               paste0("Dim ", y_dim(), " coord"),
                               paste0("Dim ", x_dim(), " cos"),
                               paste0("Dim ", y_dim(), " cos")),
                  extensions = base::list("ColReorder" = NULL,
                                          "Buttons" = NULL,
                                          "FixedColumns" = base::list(leftColumns=1)),
                  options = list(
                    dom = "BRrltpi",
                    autoWidth = TRUE,
                    buttons =
                      list("copy",
                           list(extend = "collection",
                                buttons = c("csv", "excel", "pdf"),
                                text = "Download"),
                           I("colvis")))) %>%
      plotly::highlight(on = "plotly_selected",
                        opacityDim = 0.01,
                        selected = attrs_selected(opacity = 1))
    
  }, server = FALSE)
}
