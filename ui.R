# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("microbePlot"),
  
  tabsetPanel(
    
    tabPanel("Data", 
      sidebarPanel(
        fileInput("metaFilename", "Select metadata file", accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("metaOptions", "Show file options"),
        conditionalPanel(
          condition = "input.metaOptions == true",
          checkboxInput('metaHeader', 'Header', TRUE),
          selectInput('metaSep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), 'Comma'),
          selectInput('metaQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), 'Double Quote')
        ),
        HTML('<hr>'),
        
        fileInput("microbeFilename", "Select microbe file", accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("microbeOptions", "Show file options"),
        conditionalPanel(
          condition = "input.microbeOptions == true",
          checkboxInput('microbeHeader', 'Header', TRUE),
          selectInput('microbeSep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), 'Comma'),
          selectInput('microbeQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), 'Double Quote'),
          checkboxInput("relativize", "Convert microbeData to relative abundance", FALSE)
        ),
        HTML('<hr>'),
        
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The microbe data file should include the abundances of the microbes in each sample."),
        helpText("Samples must be in rows, with the same order in both files.")

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
    ),
    # help
    tabPanel("Help",
      sidebarPanel(
        selectInput("helpTopic", "Help topic", choices=c("About", "Data", "Histogram", "Scatter", "PCA", "Bar plot", "Cluster", "WGCNA"))
      ),
      mainPanel(
        
      )
    )
  ),  

  mainPanel()
))
        