# Created 3/6/2013 by Daniel Beck
# User interface file for shiny

library(shiny)

plotHeight="600px"

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
        HTML('<br>'),
        
        selectInput("dataTransform", "Transform data:", 
                     list("None"="none",
                          "Relative abundance" = "total",
                          "Presence/absence" = "pa",
                          "Hellinger" = "hellinger"
                     )
        ),
        
        radioButtons("saveType", "Save plots as:", 
                     list("PDF" = "pdf",
                          "PNG" = "png")
        ),
        HTML('<hr>'),
        helpText("Input data must be in two files. The metadata file should include sample information. 
                  The microbe file should include the abundances of the microbes in each sample. 
                  Samples must be in rows.")
      ),
             
      mainPanel(
        tableOutput("viewMetaData"),
        tableOutput("viewMicrobeData"),
        HTML("<br>"),
        helpText("This Venn diagram shows the number of samples in each file. Only the overlapping 
                 samples are retained for use by microbePlot."),
        plotOutput("vennPlot", height="300px")
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
    
    # PCA tab
    tabPanel("PCA",
             uiOutput("pcaVariableSelection"),
             
             mainPanel(
               plotOutput("pcaPlot", height=plotHeight)
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
    # heatmap
    tabPanel("Heatmap",
      uiOutput("heatmapVariableSelection"),
      mainPanel(
        plotOutput("heatmapPlot", height=plotHeight)
      )
    ),
    # stacked bar plot
    tabPanel("Stacked bar plot",
             uiOutput("stackedbarVariableSelection"),
             mainPanel(
               plotOutput("stackedbarPlot", height=plotHeight)
             )
    ),    
 
    # help
    tabPanel("Help",
      sidebarPanel(
        selectInput("helpTopic", "Help topic", choices=c("About", "Data", "Histogram", "Scatter", 
                    "PCA", "Bar plot", "Cluster", "WGCNA", "Heatmap", "Color"))
      ),
      mainPanel(
        conditionalPanel(
          condition = "input.helpTopic == 'About'",
          helpText("microbePlot is a tool for visualizing microbial community data. It is currently being developed
                   and may include errors in both plots and analyses. Any results provided by microbePlot should be 
                   used with caution."),
          helpText("Comments, suggestions, or bug reports are welcome. I can be reached by email at danlbek@gmail.com
                   The source code and install instructions for microbePlot are available at https://github.com/danlbek/microbePlot."),
          helpText("microbePlot is licensed under GPLv3")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Data'",
          helpText("microbePlot requires input data in specific formats. The data should be in two separate files.
                   The metaData file should include sample information. Samples should be in rows, with sample names
                   in the first column. Additional columns may include any information about the samples. The microbeData
                   file should include microbial abundance data. Samples should be in rows, with sample names in the
                   first column. Additional columns should list taxa abundance. The abundance may be presence/absence,
                   relative abundance, or count data. microbePlot assumes files are in CSV format by default. Other
                   file types may be selected by checking the 'Show file options' box."),
          helpText("Files that are correctly read by microbePlot will be shown in the main pannel. Only the first five
                   lines will be displayed. When the microbeData file is loaded, microbePlot automatically calculates 
                   several diversity indices and adds them to the metaData information. These diversity calculations include
                   Shannon, Simpson, and Inverse Simpson indices. They are calculated using the vegan R package and may
                   not be appropriate for every type of data."),
          helpText("microbePlot includes options for converting abundance data to relative abundance or to presence/absence.
                   Relative abundance is calculated by dividing each taxa count by the sum of all taxa counts for each sample.
                   Presence/absence converts every non-zero taxa to 1. The converted data is used for all visualizations and
                   summaries."),
          helpText("microbePlot supports saving plots in two different file formats. The user may select either PDF or PNG formats.
                   All saved plots will be in the format selected in the 'Data' tab.")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Histogram'",
          helpText("The histogram panel displays a histogram of the selected variable. Additional summary information about
                   the variable is shown in the sidebar. The 'breaks' slider controls the aproximate number of histogram bars.")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Scatter'",
          helpText("The scatter panel displays a scatter plot of two selected variables. The X variable is displayed on the X
                   axis and the Y variable on the Y axis. The color variable and color options control the color of the plotted 
                   points. See the 'Color' help topic for details.")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'PCA'",
          helpText("The PCA panel displayes a scatter plot of the samples using PCA to define point positions. microbePlot 
                   calculates the principal components using R's princomp command. Additionally, the percent variation explained
                   by each principal component axis is shown in parentheses. Any principal component may be plotted on the X or Y
                   axes by selecting the appropriate number in the sidebar. The color variable and color options control the 
                   color of the plotted points. See the 'Color' help topic for details.")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Bar plot'",
          helpText("The bar plot panel displayes a bar plot using variables selected in the sidebar. The value variable defines 
                   the height of the bars. The bar variable defines the bar categories. The bar variable may be reduced to categories
                   by selecting 'Categorize bar variable'. This option breaks the bar variable into categories of equal width.
                   The number of samples represented by each bar is shown in blue. Error bars show 95% confidence intervals, assuming
                   that the samples are independent and normally distributed. NA values are removed from the analysis.")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'Cluster'",
          helpText("The cluster panel displays a clustering dendrogram of the samples. The distance between samples is calculated 
                   using the selected 'Distance method'. The samples are then clustered using the selected 'Cluster method'.
                   The color variable and color options control the color bar below the dendrogram. See the 'Color' help topic for details.
                   In order to focus on specific branches of the dendrogram, the subtree cutoff may be used. The 'Subtree cut height'
                   can be used to set the level at which the tree should be separated. The subtrees may then be viewed in the 'Subtrees' tab.
                   By default microbePlot uses the taxa abundances to calculate the differences between samples. This may be changed by selecting
                   different features to define the samples. Options include using the metadata, using a subset of the taxa abundances, or
                   selecting some combination of the two.")
        ),        
        conditionalPanel(
          condition = "input.helpTopic == 'WGCNA'",
          helpText("WGCNA is an analysis technique developed to detect genetic interactions that may be related to phenotypic characteristics.
                   microbePlot adopts the same analysis for detecting groups of microbes that may be related to sample information. The first
                   step in the WGCNA analysis is to determine the pairwise correlation between each taxa. The user may select from three 
                   correlation methods. Correlations that result in NA values are set to 0. The correlations are then converted to distances by 
                   subtracting the magnitude of the correlation from 1. The microbes are then clustered using average distances. The resulting 
                   dendrogram is cut into groups at the 'Cut-off level'. The correlated groups of taxa may be viewed as a heatmap in the 'Heatmap'
                   panel. The groups of taxa are next reduced to a single feature. Uncorrelated microbes are also lumped into a single group. 
                   The Pearson correlation between these groups and the sample information is then calculated. A heatmap of these correlations may be
                   viewed in the 'Correlations' panel. It should be noted that the validity of this analysis on ecological data has not been verified.
                   ")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Heatmap'",
          helpText("The heatmap and associated sample clustering are based on the microbial data. 
                   Only a subset of microbes are used to generate the heatmap. The microbes are ranked by
                   the sum of their abundance data. The highest ranked X microbes are used in the heatmap, where
                   X is selected by the user.")
        ),
        conditionalPanel(
          condition = "input.helpTopic == 'Color'",
          helpText("Many plot color points or bars based on a user selected variable. This allows the user to view the distribution of the color
                   variable across the samples, effectively adding another dimension to the plot. There are three main options, each of which 
                   may be appropriate in different situations. The 'Unique' option generates a different color for each different value of the variable.
                   This option works well up to about ten values. When there are many more values the selected colors are still unique, however, they
                   become difficult to tell apart. The 'Gradient' option may be used when the color variable is continuous. This option selects colors
                   from a gradient. This makes it easy to distinguish between low and high valued colors. The third option is to 
                   use color 'Categories'. This option breaks the color variable up into a user selected number of groups. These groups are then 
                   given unique colors.")
        )
      )
    )
  ),  

  mainPanel()
))
        
