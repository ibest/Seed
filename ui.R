# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)
library(shinyIncubator)

plotHeight="600px"
heatmapPlotHeight="1000px"

shinyUI(
  pageWithSidebar(
  

  # Application title
  headerPanel("Seed"),
  
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

       #checkboxInput("advancedOptions", "Show advanced options", value=FALSE),
       # conditionalPanel(
       #   condition = "input.advancedOptions == true",
       #   checkboxInput("dataNamesLimit", "Use only the 20 most abundant taxa for plot annotation
       #                 (metadata variables will always be used)", value=TRUE),
          numericInput("cutoffPercent", "Combine taxa representing less than a of 
                        total counts (default 0 does nothing):",
                        value = 0, min = 0, max = 100, step = 0.00001),
          helpText("Combined taxa will be labelled 'other_combined'"),            
          radioButtons("saveType", "Save plots as:", 
                       list("PDF" = "pdf",
                            "PNG" = "png")
        #              )
        ),
        HTML('<hr>'),
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The taxa file should include the abundances of the taxa in each sample. 
                  Samples must be in rows."
                  ),
	      helpText("Note: Diversity indices are calculated using a relative abundance transformation of the original data."),
          checkboxInput("loadDemo", "Load a demonstration dataset", value=FALSE)
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
                 tabPanel("Dendrogram", 
                          value="ndendrogram",
                       #  textOutput("wgcnaErrorOut"),
                          plotOutput("dendroPlot", height=plotHeight)
                 ),
                 tabPanel("Heatmap", 
                          value="nheatmap",
                          plotOutput("htmpPlot", height=plotHeight)
                        
                 ),
                 tabPanel("Correlations", 
                          value="ncorrelations",
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
                    "PCoA", "Bar plot", "Cluster", "WGCNA", "Heatmap", "Color",
                    "AE35")
              ),
         # the following whimsy actually serves a purpose; for an unknown reason
         # the help page is truncated without something to fill the space downward
         conditionalPanel(condition = "input.helpTopic == 'AE35'",
            img(src = "HAL9000.svg"),
            helpText(""),
            helpText("It can only be attributable to human error.")
         ),
         conditionalPanel(condition = "input.helpTopic != 'AE35'",
            img(src = "clippy.jpg")
         )

      ),
      mainPanel(
        conditionalPanel(
          condition = "input.helpTopic == 'About'",
          includeHTML("./www/help_about.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Data'",
          includeHTML("./www/help_data.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Histogram'",
          includeHTML("./www/help_hist.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Scatter'",
          includeHTML("./www/help_scatter.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'PCoA'",
          includeHTML("./www/help_PCoA.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Bar plot'",
          includeHTML("./www/help_barplot.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Cluster'",
          includeHTML("./www/help_cluster.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'WGCNA'",
          includeHTML("./www/help_WGCNA.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Heatmap'",
          includeHTML("./www/help_heatmap.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Color'",
          includeHTML("./www/help_color.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'AE35'",
          helpText("I've just picked up a fault in the AE-35 unit.")
        )
      )
    )
  ),  

  mainPanel()

))
        
