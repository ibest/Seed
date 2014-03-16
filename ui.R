# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)
library(shinyIncubator)

plotHeight="600px"
heatmapPlotHeight="1000px"

shinyUI(
  pageWithSidebar(
  

  # Application title
  headerPanel("Seed",   tags$link(rel = 'stylesheet', type = 'text/css', href = 'style.css')),
  
  
  tabsetPanel( 
    
    tabPanel("Data", 
      progressInit(), 
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
        
        fileInput("microbeFilename", "Select taxa file", accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("microbeOptions", "Show file options"),
        conditionalPanel(
          condition = "input.microbeOptions == true",
          checkboxInput('microbeHeader', 'Header', TRUE),
          selectInput('microbeSep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), 'Comma'),
          selectInput('microbeQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), 'Double Quote')
        ),
        HTML('<br>'),
       
        selectInput("dataTransform", "Transform data:", 
                     list("None"="none",
                          "Relative abundance" = "total",
                          "Presence/absence" = "pa",
                          "Hellinger" = "hellinger"
                     )
        ),
        checkboxInput("loadDemo", "Load demonstration dataset", value=FALSE),
        checkboxInput("advancedOptions", "Show advanced options", value=FALSE),
        
        conditionalPanel(
          condition = "input.advancedOptions == true",
          checkboxInput("dataNamesLimit", "Limit plot feature options to most abundant taxa", value=TRUE),
          numericInput("cutoffPercent", "Combine taxa representing less than a given percentage of total counts:",
                       value = 0, min = 0, max = 100, step = 0.00001),
          helpText("Combined taxa will be labelled 'other_combined'"),            
          radioButtons("saveType", "Save plots as:", 
                       list("PDF" = "pdf",
                            "PNG" = "png")
                      )
        ),
        HTML('<hr>'),
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The taxa file should include the abundances of the taxa in each sample. 
                  Samples must be in rows."
                  ),
	      helpText("Note: Diversity indices are calculated using a relative abundance transformation of the original data.")
        ),
             
      mainPanel(
        tableOutput("viewMetaData"),
        tableOutput("viewMicrobeData"),
        HTML("<br>"),
        helpText("This Venn diagram shows the number of samples in each file. 
                  Only the overlapping samples are retained for use by Seed."
                  ),
        plotOutput("vennPlot", height="300px"),
        textOutput("dimRawMeta"),
        textOutput("dimRawMicrobe"),
        textOutput("dimPostMeta"),
        textOutput("dimPostMicrobe")
      )
    ),
    
    # Histogram tab
    tabPanel("Histogram", 
      uiOutput("histVariableSelection"),
      mainPanel(    
        plotOutput("histPlot", height=plotHeight)
      )
    ),
    
    # Scatterplot tab
    tabPanel("Scatter",
      uiOutput("scatterVariableSelection"),
      mainPanel(
        plotOutput("scatterPlot", height=plotHeight)
      )
    ),
    


      
    # PCoA tab
    tabPanel("PCoA", 
             uiOutput("pcoaVariableSelection"),
             mainPanel(
                 plotOutput("pcoaPlot", height=plotHeight)
             )
             
    ),
    
    
    # bar charts
    tabPanel("Bar plot",
             uiOutput("barVariableSelection"),
             mainPanel(
               plotOutput("barPlot", height=plotHeight)
             )
    ),
    tabPanel("Cluster", 
             uiOutput("clusterVariableSelection"),
             mainPanel(
               tabsetPanel(
                 id = "clusterTab",
                 tabPanel("Complete", value="complete",
                          plotOutput("clusterPlot", height=plotHeight)
                 ),
                 tabPanel("Subtrees", value="subtree",
                          plotOutput("clusterGroupPlot", height=plotHeight)
                 ),
                 tabPanel("Silhouette", value="silhouette",
                          plotOutput("silhouettePlot", height=plotHeight)
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
                          plotOutput("dendroPlot", height=plotHeight)
                 ),
                 tabPanel("Heatmap", value="nheatmap",
                          plotOutput("htmpPlot", height=plotHeight)
                 ),
                 tabPanel("Correlations", value="ncorrelations",
                          plotOutput("corPlot", height=plotHeight)
                 )
                 
               )
             )
    ),
    
        # stacked bar plot
    tabPanel("Stacked bar plot",
             uiOutput("stackedbarVariableSelection"),
             mainPanel(
               plotOutput("stackedBarPlot", height=plotHeight)
             )
    ),    
    
    # heatmap
    tabPanel("Heatmap",
      uiOutput("heatmapVariableSelection"),
      mainPanel(
        plotOutput("heatmapPlot", height=heatmapPlotHeight)
      )
    ),

 
    # help
    tabPanel("Help",
      sidebarPanel(
        selectInput(
              "helpTopic", 
              "Help topic", 
              choices=c("About", "Data", "Histogram", "Scatter", 
                    "PCoA", "Bar plot", "Cluster", "WGCNA", "Heatmap", "Color")
              )
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.helpTopic == 'About'",
          includeHTML("./html/help_about.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Data'",
          includeHTML("./html/help_data.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Histogram'",
          helpText("./html/help_hist.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Scatter'",
          includeHTML("./html/help_scatter.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'PCoA'",
          includeHTML("./html/help_PCoA.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Bar plot'",
          includeHTML("./html/help_barplot.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Cluster'",
          includeHTML("./html/help_cluster.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'WGCNA'",
          includeHTML("./html/help_WGCNA.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Heatmap'",
          includeHTML("./html/help_heatmap.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Color'",
          includeHTML("./html/help_color.html")
        )
      )
    )
  ),  

  mainPanel()

))
        
