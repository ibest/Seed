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
          selectInput('microbeQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), 'Double Quote')
        ),
        HTML('<hr>'),
        
        radioButtons("dataConvert", "Convert data to:", 
                     list("None"="none",
                          "Relative abundance" = "relative",
                          "Presence/absence" = "presence")
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
                 tabPanel("Complete", value="complete",
                          plotOutput("clusterPlot", height="600px")
                 ),
                 tabPanel("Subtrees", value="subtree",
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
    # heatmap
    tabPanel("Heatmap",
      uiOutput("heatmapVariableSelection"),
      mainPanel(
        plotOutput("heatmapPlot", height="600px")
      )
    ),
    
    # help
    tabPanel("Help",
      sidebarPanel(
        selectInput("helpTopic", "Help topic", choices=c("About", "Data", "Histogram", "Scatter", "PCA", "Bar plot", "Cluster", "WGCNA", "Heatmap"))
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.helpTopic == 'About'",
          helpText("microbePlot is a tool for visualizing microbial community data. It is currently being developed
                   and may include errors in both plots and analyses. Any results provided by microbePlot should be 
                   used with caution. The source code and install instructions are available at https://github.com/danlbek/microbePlot. 
                   Comments, suggestions, or bug reports are welcome.")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Data'",
          helpText("About microbePlot")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Histogram'",
          helpText("About microbePlot")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Scatter'",
          helpText("About microbePlot")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'PCA'",
          helpText("About microbePlot")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Bar plot'",
          helpText("About microbePlot")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Cluster'",
          helpText("About microbePlot")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'WGCNA'",
          helpText("About microbePlot")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Heatmap'",
          helpText("The heatmap and associated sample clustering are based on the microbial data. 
                   Only a subset of microbes are used to generate the heatmap. The microbes are ranked by
                   the sum of their abundance data. The highest ranked X microbes are used in the heatmap, where
                   X is selected by the user.")
        )
      )
    )
  ),  

  mainPanel()
))
        
