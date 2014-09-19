# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)


plotHeight="600px"
heatmapPlotHeight="1000px"

shinyUI(
  pageWithSidebar(
  

  # Application title
  headerPanel("Seed"),
  
  tabsetPanel( 
    
    tabPanel("Data", 

      sidebarPanel(

        fileInput("metaFilename", "Select metadata file", accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("metaOptions", "Show file options"),
        conditionalPanel(
          condition = "input.metaOptions == true",
          checkboxInput('metaHeader', 'Header', TRUE),
          radioButtons('metaSep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
          radioButtons('metaQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"')
        ),
        HTML('<hr>'),
        
        fileInput("microbeFilename", "Select taxa file", accept=c('text/csv', 'text/comma-separated-values,text/plain')),
        checkboxInput("microbeOptions", "Show file options"),
        conditionalPanel(
          condition = "input.microbeOptions == true",
          checkboxInput('microbeHeader', 'Header', TRUE),
          radioButtons('microbeSep', 'Separator', c(Comma=',', Semicolon=';', Tab='\t'), ','),
          radioButtons('microbeQuote', 'Quote', c(None='', 'Double Quote'='"', 'Single Quote'="'"), '"')
        ),
        HTML('<br>'),
       
        selectInput("dataTransform", "Transform data:", 
                     list("None"="none",
                          "Relative abundance" = "total",
                          "Presence/absence" = "pa",
                          "Hellinger" = "hellinger"
                     )
        ),

        checkboxInput("advancedOptions", "Show advanced options", value=FALSE),
        
        conditionalPanel(
          condition = "input.advancedOptions == true",
          radioButtons("saveType", "Save plots as:", 
                       list("PDF" = "pdf",
                            "PNG" = "png")
                      ),
          downloadLink('downloadMetaData', 'Download metadata'),
	  HTML('<br>'),
          downloadLink('downloadAbundanceData', 'Download abundance data')

        ),
        HTML('<hr>'),
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The taxa file should include the abundances of the taxa in each sample. 
                  Samples must be in rows."
                  ),
	      helpText("Note: Diversity indices are calculated using a relative abundance transformation of the original data."),
          checkboxInput("loadDemo", "Load a demonstration dataset (Ravel, et al. 2011)", value=FALSE)
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
    
    # histogram
    tabPanel("Histogram", 
      uiOutput("histVariableSelection"),
      mainPanel(    
        plotOutput("histPlot", height=plotHeight)
      )
    ),
    
    # scatterplot
    tabPanel("Scatter",
      uiOutput("scatterVariableSelection"),
      mainPanel(
        plotOutput("scatterPlot", height=plotHeight)
      )
    ),
    

    # bar charts
    tabPanel("Bar plot",
             uiOutput("barVariableSelection"),
             mainPanel(
               plotOutput("barPlot", height=plotHeight)
             )
    ),

    # stacked bar plot
    tabPanel("Stacked bar plot",
             uiOutput("stackedbarVariableSelection"),
             mainPanel(
               plotOutput("stackedBarPlot", height=plotHeight)
             )
    ),    
   
    # PCoA tab
    tabPanel("PCoA", 
             uiOutput("pcoaVariableSelection"),
             mainPanel(
                 plotOutput("pcoaPlot", height=plotHeight)
             )
             
    ),

    # cluster dendrogram
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
    
    # heatmap
    tabPanel("Heatmap",
      uiOutput("heatmapVariableSelection"),
      mainPanel(
        plotOutput("heatmapPlot", height=heatmapPlotHeight)
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
    
    # help
    tabPanel("Help",
      sidebarPanel(
        height = 12,
        selectInput(
              "helpTopic", 
              "Help topic", 
              choices=c("About", "Data", "Histogram", "Scatter", 
                    "PCoA", "Bar plot", "Cluster", "WGCNA", "Heatmap", "Color",
                    "AE35")
              ),
              
         conditionalPanel(condition = "input.helpTopic == 'AE35'",
           img(src = "HAL9000.svg"),
           helpText(""),
           helpText("It can only be attributable to human error.")
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
          includeHTML("./www/help_hist.html"),
          HTML('<br><br><br><br><br><br><br>
                <br><br><br><br><br><br><br><br>') 
            # allows menus to display properly given short main body
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Scatter'",
          includeHTML("./www/help_scatter.html"),
          HTML('<br><br><br><br><br><br><br>
                <br><br><br><br><br><br><br><br>') 
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'PCoA'",
          includeHTML("./www/help_PCoA.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Bar plot'",
          includeHTML("./www/help_barplot.html"),
          HTML('<br><br><br><br><br><br><br>
                <br><br><br><br><br><br><br><br>') 
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Cluster'",
          includeHTML("./www/help_cluster.html")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'WGCNA'",
          includeHTML("./www/help_wgcna.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Heatmap'",
          includeHTML("./www/help_heatmap.html")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Color'",
          includeHTML("./www/help_color.html")
        )
      )
    )
  ),  

  mainPanel()

))
        
