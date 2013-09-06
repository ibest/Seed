# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("microbePlot"),
  
  tabsetPanel(
    
    tabPanel("Data", 
#      uiOutput("fileSelection"),
      sidebarPanel(
        fileInput("metaFilename", "Select metadata file", 
                  accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        fileInput("microbeFilename", "Select microbe file", 
                  accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("relativize", "Convert microbeData to relative abundance", FALSE),
        HTML('<br>'),
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The microbe data file should include the abundances of the microbes in each sample."),
        helpText("The files must be in CSV format with equal numbers of samples. Samples must be in rows, with the same order in both files.")

      ),
             
      mainPanel(
        tableOutput("viewMetaData"),
        tableOutput("viewMicrobeData")
      )
    ),
    
    # Histogram tab
    tabPanel("Histogram", 

      uiOutput("histVariableSelection"),
  
      mainPanel(    
        plotOutput("histPlot", height="600px")
      )
    ),
    
    # Scatterplot tab
    tabPanel("Scatter",
      uiOutput("scatterVariableSelection"),
      
      mainPanel(
        plotOutput("scatterPlot", height="600px")
      )
    ),
    
    # PCA tab
    tabPanel("PCA",
             uiOutput("pcaVariableSelection"),
             
             mainPanel(
               plotOutput("pcaPlot", height="600px")
             )
    ),
    
    # bar charts
    tabPanel("Bar plot",
             uiOutput("barVariableSelection"),
             
             mainPanel(
               plotOutput("barPlot", height="600px")
             )
    ),
    tabPanel("Cluster", 
             uiOutput("clusterVariableSelection"),
             
             mainPanel(
               tabsetPanel(
                 id = "clusterTab",
                 tabPanel("Complete", 
                          plotOutput("clusterPlot", height="600px")
                 ),
                 tabPanel("Subtrees",
                          plotOutput("clusterGroupPlot", height="600px")
                 )
              )
             )
             
    ),
    # wgcna
    tabPanel("WGCNA",
             uiOutput("wgcnaVariableSelection"),
             mainPanel(
               tabsetPanel(
                 id="wgcnaTab",
                 tabPanel("Dendrogram", value="ndendrogram",
                          plotOutput("dendroPlot", height="600px")
                 ),
                 tabPanel("Heatmap", value="nheatmap",
                          plotOutput("htmpPlot", height="600px")
                 ),
                 tabPanel("Correlations", value="ncorrelations",
                          plotOutput("corPlot", height="600px")
                 )
                 
               )
             )
    )
  ),  
  mainPanel()
))
        