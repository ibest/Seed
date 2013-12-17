
# Created 3/6/2013 by Daniel Beck
# Server file for shiny

library(shiny)
library(vegan)
library(WGCNA)
library(gplots)

shinyServer(function(input, output) {

##################################################################################################*
###################################### FUNCTIONS #################################################*
##################################################################################################*

  # functions used by multiple plots

  # color function wrappers
  gradientColors<-function(n){
    colorRampPalette(c("blue", "red"))(n)
  }
  discreteColors<-function(n){
    trim<-as.integer(n*0.1)+3   # prevents color overlap at ends of rainbow
    rainbow(n+trim)[1:n]
  }
  
  # generate color vector for plots. Returns both color vector and value vector 
  #(for use in legends).
  # numeric values are truncated to 6 digits
  getColor<-function(featureVector, type="unique", numCat=1){
    if (is.numeric(featureVector)) featureVector<-round(featureVector, 3)
    if (type == "unique"){
      valueVector<-as.factor(featureVector)
      if (length(levels(valueVector))<=8){
        colorVector<-as.numeric(featureVector)
      }else{
        colorVector<-discreteColors(length(levels(valueVector)))[as.numeric(valueVector)]
      }
    }
    if (type == "gradient"){
      # assign colors to variables in a linear manner
      valueVector<-featureVector
      colorVector<-gradientColors(100)[as.numeric(cut(as.numeric(featureVector), 100))]
    }
    if (type == "category"){
      valueVector<-cut(as.numeric(featureVector), numCat)
      if (numCat==1){colorVector<-1}
      if (numCat<=8 && numCat>1){ colorVector<- as.numeric(valueVector) }
      if (numCat>8){ colorVector<-discreteColors(numCat)[as.numeric(valueVector)] }
    }
    return(list(colorVector, valueVector))
  }

  # plot legend
  plotLegend<-function(colorVector, valueVector, gradient=F, 
                       title="", min=0, max=1, cex=1, keyCol=3) {
    if (gradient){
      par(mar=c(6,4,2,2)+0.1)
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = title, cex.main=cex)
      legend_image <- as.raster(matrix(gradientColors(100), ncol=100))
      rasterImage(legend_image, 0, 0, 1, 1)
      axis(side=1, at=seq(0,1,l=5), labels=seq(min,max,l=5),col.axis="black", cex.axis=cex)
      mtext("Key", side=2, las=2, cex=cex)
      
    }else{
      par(mar=c(0,2,0,2)+0.1)
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
      uniquePairs<-match(levels(as.factor(valueVector)), valueVector)
      values<-levels(as.factor(valueVector))
      colors<-colorVector[uniquePairs]
      legend("center", legend=values, fill=colors, 
             ncol=keyCol, 
             box.col="white", title=title,
             cex=cex)
    }
  }

  # tests to see if file is square (if it has column names for every column)
  isSquare<-function(filename, sep){
    test<-sapply(sapply(readLines(filename,2), function(i) strsplit(i,sep)), length)
    return(test[1]==test[2])
  }

# This isn't working for some reason
#   # save plot function
#   savePlot<-function(fileName="error", plotFunction, fileExt){
#       fileName<-paste(fileName, fileExt, sep=".")
#       if (fileExt=="png"){
#         png(fileName, width=2000, height=2000, units="px")
#         plotFunction()
#         print(fileName)
#         dev.off()
#       }
#       if (fileExt=="pdf"){
#         pdf(fileName, width=10, height=10)
#         plotFunction()
#         print(fileName)
#         dev.off()
#       }
#   }


###################################################################################################
############################################# DATA ################################################
###################################################################################################
  output$vennPlot <- renderPlot({          
    par(mar=c(0,0,0,0))
    vl<-list(microbeData=row.names(inputMicrobeData()), metaData=row.names(inputMetaData()))
    names(vl)<-c("Samples in microbe file", "Samples in metadata file")
    venn(vl)
  })

  # microbeData will contain relative abundances
  inputMicrobeData<-reactive({
    microbeFile <- input$microbeFilename$datapath
    if (is.null(microbeFile)) return(NULL)
    microbeData <- read.csv(microbeFile, header=input$microbeHeader, 
                            sep=input$microbeSep, quote=input$microbeQuote) 
    if (isSquare(microbeFile, input$microbeSep)){
      # This is a strange way to do this, but it fixes single column file quirks
      rn<-microbeData[,1]
      cn<-colnames(microbeData)
      microbeData<-as.data.frame(microbeData[,-1])
      row.names(microbeData)<-rn
      colnames(microbeData)<-cn[-1]
    }
    microbeData
  })

  microbeData <- reactive({ 
    microbeData<-inputMicrobeData()
    # if metaData is available, use only samples that overlap
    if (!is.null(input$metaFilename$datapath)){
      microbeData<-microbeData[row.names(microbeData)%in%row.names(inputMetaData()),]
      microbeData<-microbeData[na.exclude(match(row.names(inputMetaData()),
                               row.names(microbeData))),  ]
    }
    if (input$dataTransform!="none"){
      microbeData <- decostand(microbeData, method=input$dataTransform)
    }

    microbeData
  })
  
  # metaData will contain all sample information other than microbial abundances
  # When reading in metadata, also calculate diversity metrics from microbeData and add to metaData
  inputMetaData<-reactive({
    metaFile <- input$metaFilename$datapath
    if (is.null(metaFile)) return(NULL)
    metaData <- read.csv(metaFile, header=input$metaHeader, 
                         sep=input$metaSep, quote=input$metaQuote)
    if (isSquare(metaFile, input$metaSep)){
      # elaborate to prevent single column file errors
      rn<-metaData[,1]
      cn<-colnames(metaData)
      metaData<-as.data.frame(metaData[,-1])
      row.names(metaData)<-rn
      colnames(metaData)<-cn[-1]
    }
    metaData
  })

  metaData <- reactive({ 
    metaData<-inputMetaData()
    if (is.null(input$microbeFilename$datapath)) return(metaData)

    metaData<-metaData[row.names(metaData)%in%row.names(inputMicrobeData()),]

    metaData<-cbind(metaData,
      Shannon.diversity=diversity(microbeData(), index="shannon"),
      Simpson.diversity=diversity(microbeData(), index="simpson"),
      inverse.Simpson.diversity=diversity(microbeData(), index="invsimpson")
    )
    metaData
  })
  
  # include all data combined in order to produce comprehensive lists of features
  allData <- reactive({ 
    if (is.null(input$microbeFilename$datapath)||is.null(input$metaFilename$datapath)) return(NULL)
    cbind(metaData(), microbeData()) 
  })

  # extract feature names from allData
  features <- reactive({ 
    colnames(allData) 
  })
  
  # display top five lines of metaData file
  output$viewMetaData <- renderTable({
    head(metaData(), n=5)
  })
  
  # display top five lines of microbeData file
  output$viewMicrobeData <- renderTable({
    head(microbeData(), n=5)
  })

###################################################################################################
########################################### HISTOGRAM #############################################
###################################################################################################

  # dynamically generate histogram UI
  output$histVariableSelection <- renderUI({  
    # generate sidebar
    sidebarPanel(
      selectInput("variable", "Variable:", choices = colnames(allData())),
      sliderInput("breaks", "Breaks:", min=2, max=20, value=10),
      HTML('<br><br>'),
      HTML('<div align="center">'),
      tableOutput("varSum"),
      HTML('<div align="right">'),
      downloadButton("saveHist", "Save Plot"),
      HTML('</div>'),
      HTML('<br><hr><div align="left">'),
      uiOutput("histPlotOptions")
    )
  })
  
  output$histPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("histPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.histPlotOptions == true",
        sliderInput("histFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("histMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("histMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("histMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("histMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("histXlab", "X label", value=input$variable),
        textInput("histYlab", "Y label", value="Number of occurances"),
        textInput("histTitle", "Title", value=paste("Histogram of ", input$variable, sep=""))
      )
    )
  })

  # generate summary of histogram variable for display in sidebar
  output$varSum <- renderTable({
    s<-as.matrix(summary(allData()[,which(colnames(allData())==input$variable)]))
    colnames(s)<-input$variable
    s
  })
  
  # histogram plot
  plotHistogram <- function(){
    if (is.null(input$breaks)) return(NULL)
    par(mar=c(input$histMarBottom,input$histMarLeft,input$histMarTop,input$histMarRight))
    hist(as.numeric(allData()[,which(colnames(allData())==input$variable)]), 
         breaks=input$breaks, 
         xlab=input$histXlab, 
         ylab=input$histYlab,
         main=input$histTitle,
         cex.axis=input$histFontSize, cex.main=input$histFontSize, cex.lab=input$histFontSize
    )
  }
  
  # save histogram plot
  output$saveHist <- downloadHandler(
    filename = function() { paste("histogramPlot", fileExtension(), sep=".") },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$histFontSize)
        plotHistogram()
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        plotHistogram()
        dev.off()
      }
    }
  )
  
  # display histogram of the requested variable
  output$histPlot <- renderPlot({
    plotHistogram()
  })
  
###################################################################################################
########################################### SCATTER PLOT ##########################################
###################################################################################################

  # generate scatter plot UI
  output$scatterVariableSelection <- renderUI({
    # generate sidebar
    sidebarPanel(
      selectInput("variable1", "X:", choices = colnames(allData())),
      selectInput("variable2", "Y:", choices = colnames(allData())),
      selectInput("scatterColorVariable", "Color variable:", choices = colnames(allData())),
      radioButtons("scatterColorType", "Color options:", 
                  list("Unique" = "unique",
                       "Gradient" = "gradient",
                       "Categories" = "category")
      ),
      conditionalPanel(
        condition = 'input.scatterColorType == "category"',
        numericInput("nscatterColorCat", "Number of categories:", 4)
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveScatter", "Save Plot"),
      HTML('</div>'),
      uiOutput("scatterPlotOptions")
    )
  })
  
  output$scatterPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("scatterPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.scatterPlotOptions == true",
        sliderInput("scatterFontSize", "Font size", min=0.01, max=3.01, value=1.5),
	      sliderInput("scatterPointSize", "Point size", min=0.01, max=3.01, value=1.0),
        sliderInput("scatterMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("scatterMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("scatterMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("scatterMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("scatterXlab", "X label", value=input$variable1),
        textInput("scatterYlab", "Y label", value=input$variable2),
        textInput("scatterTitle", "Title", 
                  value=paste("Scatterplot of", input$variable2, 
                              "vs.", input$variable1, sep=" ")),
        textInput("scatterKeyTitle", "Legend title", value=input$scatterColorVariable),
        sliderInput("scatterKeyFontSize", "Legend font size", min=0.01, max=3.01, value=1.5),
        sliderInput("scatterKeyColumns", "Number of legend columns", min=1, max=15, value=3)

      )
    )
  })

  scatterCVlist<-reactive({
      colorVariable<-which(colnames(allData())==input$scatterColorVariable)
      CVlist <- getColor(allData()[,colorVariable], type=input$scatterColorType, 
                         numCat=input$nscatterColorCat)
      CVlist
  })

  scatterX<-reactive({
      as.numeric(allData()[,which(colnames(allData())==input$variable1)])
  })

  scatterY<-reactive({
      as.numeric(allData()[,which(colnames(allData())==input$variable2)])
  })

  # generate scatter plot
  plotScatter<-function(){
    layout(matrix(c(1,2,1,2),ncol=2), height = c(4,1),width = c(4,4))
    par(mar=c(input$scatterMarBottom,input$scatterMarLeft,
              input$scatterMarTop,input$scatterMarRight))
    colorV <- scatterCVlist()[[1]]
    valueV <- scatterCVlist()[[2]]
    plot(scatterX(), 
         scatterY(), 
         xlab=input$scatterXlab, 
         ylab=input$scatterYlab, 
         main=input$scatterTitle,
         col=colorV, pch="O",
         cex.axis=input$scatterFontSize, cex.main=input$scatterFontSize, 
         cex.lab=input$scatterFontSize, cex=input$scatterPointSize
    )
    plotLegend(colorV, valueV, gradient=(input$scatterColorType=="gradient"), 
               title=input$scatterKeyTitle, min=min(as.numeric(valueV), na.rm=T), 
               max=max(as.numeric(valueV), na.rm=T), cex=input$scatterKeyFontSize,
               keyCol=input$scatterKeyColumns
    )
  }
  
  # save scatter plot
  output$saveScatter <- downloadHandler(
    filename = function() { paste("scatterPlot", fileExtension(), sep=".") },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$scatterFontSize)
        plotScatter()
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        plotScatter()
        dev.off()
      }
    }
  )
  
  # Generate scatter plot of the requested variables
  output$scatterPlot <- renderPlot({
    plotScatter()
  })
  
###################################################################################################
############################################## PCA ################################################
###################################################################################################

  # generate PCA UI
  output$pcaVariableSelection <- renderUI({
    sidebarPanel(
      numericInput("pcX", "Principal component X:", 1),
      numericInput("pcY", "Principal component Y:", 2), 
      selectInput("pcaColorVariable", "Color variable:", choices = colnames(allData())),
      radioButtons("pcaColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Categories" = "category")
      ),
      conditionalPanel(
        condition = 'input.pcaColorType == "category"',
        numericInput("npcaColorCat", "Number of categories:", 4)
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("savePca", "Save Plot"),
      HTML('<br><br>'),
      downloadButton("savePcaEigen", "Save Eigenvectors"),
      HTML('</div>'),
      uiOutput("pcaPlotOptions")
    )
  })

  output$pcaPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("pcaPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.pcaPlotOptions == true",
        sliderInput("pcaFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("pcaPointSize", "Point size", min=0.01, max=3.01, value=1),
        sliderInput("pcaMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("pcaMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("pcaMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("pcaMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("pcaXlab", "X label", value=paste("Principal component ", 
                  input$pcX, " (", pcaPV()[input$pcX], "%)", sep="")),
        textInput("pcaYlab", "Y label", value=paste("Principal component ", 
                  input$pcY, " (", pcaPV()[input$pcY], "%)", sep="")),
        textInput("pcaTitle", "Title", value=paste("Scatter plot of principal components")),
        textInput("pcaKeyTitle", "Legend title", value=input$pcaColorVariable),
        sliderInput("pcaKeyFontSize", "Legend font size", min=0.01, max=3.01, value=1.5),
        sliderInput("pcaKeyColumns", "Number of legend columns", min=1, max=15, value=3)
        
      )
    )
  })

  pcaEigen<-reactive({
    covMat<-cov(microbeData())    # covariance matrix
    eigenMat<-as.matrix(eigen(covMat)$vectors)    # eigenvectors
    row.names(eigenMat)<-row.names(covMat)
    colnames(eigenMat)<-paste("ev", 1:ncol(eigenMat), sep="")
    eigenMat
  })

  pcaObject<-reactive({
    as.matrix(microbeData()) %*% pcaEigen()
  })
  pcX<-reactive({
    pcaObject()[,input$pcX]
  })
  pcY<-reactive({
    pcaObject()[,input$pcY]
  })
  # % variation explained
  pcaPV<-reactive({
    vars<-apply(pcaObject(), 2, sd)^2
    round(vars/sum(vars)*100, digits=2)
  })

  pcaCVlist<-reactive({
     colorVariable<-which(colnames(allData())==input$pcaColorVariable)
     CVlist<-getColor(allData()[,colorVariable], 
                      type=input$pcaColorType, numCat=input$npcaColorCat)
     CVlist
  })

  # generate PCA plot
  plotPca <- function(){
    layout(matrix(c(1,2,1,2),ncol=2), height = c(4,1),width = c(4,4))
    par(mar=c(input$pcaMarBottom,input$pcaMarLeft,input$pcaMarTop,input$pcaMarRight))
    
    colorV <- pcaCVlist()[[1]]
    valueV <- pcaCVlist()[[2]]
    plot(pcX(), pcY(), 
         xlab=input$pcaXlab, 
         ylab=input$pcaYlab, 
         main=input$pcaTitle,
         col=colorV, pch="O",
         cex.axis=input$pcaFontSize, cex.main=input$pcaFontSize, 
         cex.lab=input$pcaFontSize, cex=input$pcaPointSize
    )
    plotLegend(colorV, valueV, gradient=(input$pcaColorType=="gradient"), 
      title=input$pcaKeyTitle, min=min(as.numeric(valueV), na.rm=T), 
      max=max(as.numeric(valueV), na.rm=T), cex=input$pcaKeyFontSize, keyCol=input$pcaKeyColumns
    )
  }

  # save PCA plot
  output$savePca <- downloadHandler(
    filename = function() { paste("PCAplot", fileExtension(), sep=".") },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$pcaFontSize)
        plotPca()
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        plotPca()
        dev.off()
      }
    }
  )

  output$savePcaEigen <- downloadHandler(
    filename = function() { paste("PCAeigenvectors", "csv", sep=".") },
    content = function(filename) {
      write.csv(file=filename, x=pcaEigen(), quote=F)
    }
  )

  # display PCA plot
  output$pcaPlot <- renderPlot({
    plotPca()
  })
  
###################################################################################################
############################################ BAR PLOT #############################################
###################################################################################################

  # dynamic bar plot UI
  output$barVariableSelection <- renderUI({
    sidebarPanel(
      selectInput("barVariable", "Value variable:", choices = colnames(allData())),
      selectInput("categoryVariable", "Bar variable:", choices = colnames(allData())),
      checkboxInput("barCat", "Categorize bar variable", FALSE),
      conditionalPanel(
        condition = "input.barCat == true",
        numericInput("nbarCat", "Number of categories:", 4)
      ),
      helpText("The number of samples in each category is shown in blue."),
      helpText("95% confidence intervals assume normality and independence."),
      helpText("NA values are automatically omitted."),
      HTML('<br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveBar", "Save Plot"),
      HTML('</div>'),
      uiOutput("barPlotOptions")
    )
  })
    
  output$barPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("barPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.barPlotOptions == true",
        sliderInput("barFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("barMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("barMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("barMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("barMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("barXlab", "X label", value=input$categoryVariable),
        textInput("barYlab", "Y label", value=input$barVariable),
        textInput("barTitle", "Title", value="")
      )
    )
  })

  plotBar <- function(){
    barVarCol<-which(colnames(allData())==input$barVariable)
    catVarCol<-which(colnames(allData())==input$categoryVariable)
    groupedValues<-split(allData()[,barVarCol],allData()[,catVarCol])
    if (input$barCat){
      newCats<-cut(allData()[,catVarCol], input$nbarCat)
      barNames<-levels(newCats)
      groupedValues<-split(allData()[,barVarCol], as.numeric(newCats))
    }
    means<-sapply(groupedValues, function(i) mean(i, na.rm=T))
    if (input$barCat){
      names(means)<-barNames
    }
    sds<-sapply(groupedValues, function(i) sd(i, na.rm=T))
    sds[is.na(sds)]<-0
    ns<-sapply(lapply(groupedValues, na.omit), length)
    errors<-qnorm(0.975)*sds/sqrt(ns)
    par(mar=c(input$barMarBottom,input$barMarLeft,input$barMarTop,input$barMarRight))
    x<-barplot(means, 
               xlab=input$barXlab, 
               ylab=input$barYlab, 
               ylim=c(min(means-errors), max(means+errors)),
               main=input$barTitle,
               xpd=F,
               col="#f5f5f5",
               cex.axis=input$barFontSize, cex.names=input$barFontSize, 
               cex.lab=input$barFontSize, cex.main=input$barFontSize
    )
    arrows(x,means+errors,x,means-errors,code=0)
    text(x,min(means-errors)+(max(means+errors)-min(means-errors))/20, 
         ns, col="blue", cex=input$barFontSize       )
  }

  output$saveBar <- downloadHandler(
    filename = function() { paste("barPlot", fileExtension(), sep=".") },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$barFontSize)
        plotBar()
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        plotBar()
        dev.off()
      }
    }
  )
  
  # Generate barplot of the requested variable
  output$barPlot <- renderPlot({
    plotBar()
  })
  
###################################################################################################
############################################# CLUSTER #############################################
###################################################################################################

  # dynamic cluster UI
  output$clusterVariableSelection <- renderUI({
    sidebarPanel(
      selectInput("distMethod", "Distance method:", 
                  choices = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", 
                              "jaccard","gower", "altGower", "morisita", "horn", 
                              "mountford", "raup", "binomial", "chao", "cao")),
      selectInput("hclustMethod", "Cluster method:", 
                  choices = c("ward", "single", "complete", "average", 
                              "mcquitty", "median", "centroid")),
      selectInput("clusterColorVariable", "Color variable:", choices = colnames(allData())),
      radioButtons("clusterColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Categories" = "category")
      ),
      conditionalPanel(
        condition = 'input.clusterColorType == "category"',
        numericInput("nclusterColorCat", "Number of categories:", 4)
      ),
      helpText("The complete tree is cut into subtrees at the red line."),
      sliderInput("clusterCutHeight", "Subtree cut height:", min=0.0, max=1.0, value=0.5),
      numericInput("clusterGroup", "Select subtree", 1),
      
      radioButtons("clusterChoice", "Select features that define samples", 
                   choices=c("Metadata", "Microbe data", "All data", "Custom"), 
                   selected="Microbe data"
      ),
      conditionalPanel(
        condition = "input.clusterChoice == 'Custom'",
        checkboxGroupInput(inputId="customClusterVariables", label="", choices=colnames(allData()))
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveCluster", "Save Plot"),
      HTML('</div>'),
      uiOutput("clusterPlotOptions")
    )
  })
  
  output$clusterPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("clusterPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.clusterPlotOptions == true",
        sliderInput("clusterFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("clusterMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("clusterMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("clusterMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("clusterMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("clusterYlab", "Y label", value="Height"),
        textInput("clusterTitle", "Title", value=""),
        textInput("clusterKeyTitle", "Legend title", value=input$clusterColorVariable),
        sliderInput("clusterKeyFontSize", "Legend font size", min=0.01, max=3.01, value=1.5),
        sliderInput("clusterKeyColumns", "Number of legend columns", min=1, max=15, value=3)
      )
    )
  })

  clusterData <- reactive({
    if (input$clusterChoice == "Metadata"){ data <- metaData() }
    if (input$clusterChoice == "Microbe data"){ data <- microbeData() }
    if (input$clusterChoice == "All data"){ data <- allData()}
    if (input$clusterChoice == "Custom"){
      data <- allData()[,match(input$customClusterVariables, colnames(allData()))]
    }
    # non-numeric values cause cluster plot to fail. This converts them to numbers
    # apply statements fail for unknown reasons
    for (column in 1:ncol(data)){
      data[,column]<-as.numeric(data[,column])
    }

    data
  })
  
  clusterDist <- reactive({
    vegdist(clusterData(), method=input$distMethod, na.rm=T)
  })
  clusterObject <- reactive({
    hclust(clusterDist(), method=input$hclustMethod)
  })
  subtreeGroups <- reactive({
    cutree(clusterObject(), h=input$clusterCutHeight*max(clusterObject()$height))
  })
  subtreeData <- reactive({
    clusterData()[subtreeGroups()==input$clusterGroup,]
  })
  subtreeDist <- reactive({
    vegdist(subtreeData(), method=input$distMethod, na.rm=T)
  })
  subtreeObject <- reactive({
    hclust(subtreeDist(), method=input$hclustMethod)
  })

  plotCompleteTree<-function(){
    colorVariable<-which(colnames(allData())==input$clusterColorVariable)
    CVlist<-getColor(allData()[,colorVariable], type=input$clusterColorType, 
                     numCat=input$nclusterColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    layout(matrix(c(1,2,3,1,2,3),ncol=2), height = c(4,1,1),width = c(4,4))
    par(mar=c(input$clusterMarBottom,input$clusterMarLeft,
              input$clusterMarTop,input$clusterMarRight))
    plotDendroAndColors(clusterObject(), 
                        colors = data.frame(colorV), 
                        dendroLabels = NULL, 
                        abHeight = input$clusterCutHeight*max(clusterObject()$height), 
                        groupLabels = "",
                        main = input$clusterTitle,
                        ylab = input$clusterYlab,
                        setLayout = FALSE, 
                        mar = c(input$clusterMarBottom, input$clusterMarLeft,
                                input$clusterMarTop, input$clusterMarRight),
                        cex.colorLabels = input$clusterFontSize, 
                        cex.dendroLabels = input$clusterFontSize,
                        cex.rowText = input$clusterFontSize, 
                        cex.axis  = input$clusterFontSize, 
                        cex.lab = input$clusterFontSize, cex.main = input$clusterFontSize
    )
    plotLegend(colorV, valueV, gradient=(input$clusterColorType=="gradient"), 
               title=input$clusterKeyTitle, min=min(as.numeric(valueV), na.rm=T), 
               max=max(as.numeric(valueV), na.rm=T), cex=input$clusterKeyFontSize, 
               keyCol=input$clusterKeyColumns
    )
  }

  plotSubTree<-function(){
    colorVariable<-which(colnames(allData())==input$clusterColorVariable)
    CVlist<-getColor(allData()[subtreeGroups()==input$clusterGroup,colorVariable], 
                     type=input$clusterColorType, numCat=input$nclusterColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    layout(matrix(c(1,2,3,1,2,3),ncol=2), height = c(4,1,1),width = c(4,4))
    plotDendroAndColors(subtreeObject(), 
                        colors = data.frame(colorV), 
                        dendroLabels = NULL, 
                        groupLabels = "",
                        main = input$clusterTitle,
                        setLayout = FALSE,
                        mar = c(input$clusterMarBottom, input$clusterMarLeft,
                                input$clusterMarTop, input$clusterMarRight),
                        cex.colorLabels = input$clusterFontSize, 
                        cex.dendroLabels = input$clusterFontSize,
                        cex.rowText=input$clusterFontSize, cex.axis=input$clusterFontSize, 
                        cex.lab=input$clusterFontSize, cex.main=input$clusterFontSize
    )  
    plotLegend(colorV, valueV, gradient=(input$clusterColorType=="gradient"), 
               title=input$clusterKeyTitle, min=min(as.numeric(valueV), na.rm=T), 
               max=max(as.numeric(valueV), na.rm=T),cex=input$clusterKeyFontSize, 
               keyCol=input$clusterKeyColumns
    )
  }

  output$clusterPlot <- renderPlot({
    plotCompleteTree()
  })

  output$clusterGroupPlot <- renderPlot({
    plotSubTree()
  })

  output$saveCluster <- downloadHandler(
    filename = function() { 
      if (input$clusterTab=="complete"){sp<-"complete"}
      if (input$clusterTab=="subtree"){sp<-"subtree"}
      paste("clusterPlot", sp, fileExtension(), sep=".") 
    },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$clusterFontSize)
        if (input$clusterTab=="complete"){plotCompleteTree()}
        if (input$clusterTab=="subtree"){plotSubTree()}
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        if (input$clusterTab=="complete"){plotCompleteTree()}
        if (input$clusterTab=="subtree"){plotSubTree()}
        dev.off()
      }
    }

  )

###################################################################################################
############################################## WGCNA ##############################################
###################################################################################################

  # render sidebars for wgcna plots
  output$wgcnaVariableSelection <- renderUI({
    sidebarPanel(
      helpText("Warning: The Kendall correlation method is very slow."),
      selectInput("corMethod", "Network correlation method:", 
                  choices = c("pearson", "spearman", "kendall")),
      sliderInput("cutLevel", "Cut-off level:", min=0.0, max=1.0, value=0.8),
      numericInput("selectGroup", "Selected group:", 1),
      HTML('<br>'),
      helpText("Note: Metadata correlations use the Pearson method."),
      helpText("NA correlations replaced with 0."),

      HTML('<hr>'),
      HTML('<div align="right">'),
      downloadButton("saveWGCNA", "Save Plot"),
      HTML('</div>')
      )
  })

  cors <- reactive({ 
    pnacors<-cor(microbeData(), method=input$corMethod) 
    pnacors[is.na(pnacors)]<-0
    pnacors
  })
  dADJ <- reactive({ as.dist(1-as.matrix(abs(cors()))) })
  hdADJ <- reactive({ hclust(dADJ(), method="average") })
  
  plotDendrogram <- function(){
    groups<-cutreeStatic(dendro=hdADJ(), minSize=3,cutHeight=input$cutLevel)
    moduleColors<-getColor(groups, "unique")[[1]]
    plotDendroAndColors(hdADJ(), 
                        colors=data.frame(moduleColors), 
                        dendroLabels=F, 
                        abHeight=input$cutLevel, 
                        main="Species dendrogram and module colors",
                        cex.colorLabels=fontSize(), cex.dendroLabels=fontSize(),
                        cex.rowText=fontSize(), cex.axis=fontSize(), cex.lab=fontSize()
    )
  }
  output$dendroPlot <- renderPlot({
    plotDendrogram()
  })
  
  plotHtmp <- function(){
    groups<-cutreeStatic(dendro=hdADJ(), minSize=3,cutHeight=input$cutLevel)
    corlist<-sapply(unique(groups), function(i) cors()[groups==i, groups==i])
    cc<-colorRampPalette(c("white", "blue"))
    heatmap.2(as.matrix(corlist[[input$selectGroup]]), 
              scale="none", 
              trace="none", 
              margins=c(12,12), 
              cexRow=fontSize(), 
              cexCol=fontSize(),
              col=cc)
  }

  output$htmpPlot <- renderPlot({
    plotHtmp()
  })
  
  plotCor <- function(){
    groups<-cutreeStatic(dendro=hdADJ(), minSize=3,cutHeight=input$cutLevel)
    MEs0 = moduleEigengenes(microbeData(), groups+1)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, metaData(), use="p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(metaData()));
    
    textMatrix = paste(signif(moduleTraitCor, 2), 
                       "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar=c(6.1,4.1,0.1,2.1))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(moduleTraitCor),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   setStdMargins = FALSE,
                   zlim = c(-1,1),
                   main = "", 
                   cex.lab.x=fontSize())    
  }

  output$corPlot <- renderPlot({
    plotCor()
  })
  
  output$saveWGCNA <- downloadHandler(
    filename = function() { 
      if (input$wgcnaTab=="ndendrogram"){sp<-"dendrogram"}
      if (input$wgcnaTab=="nheatmap"){sp<-"heatmap"}
      if (input$wgcnaTab=="ncorrelations"){sp<-"correlations"}
      paste("wgcnaPlot", sp, fileExtension(), sep=".") 
    },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$fontSize())
        if (input$wgcnaTab=="ndendrogram"){plotDendrogram()}
        if (input$wgcnaTab=="nheatmap"){plotHtmp()}
        if (input$wgcnaTab=="ncorrelations"){plotCor()}
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        if (input$wgcnaTab=="ndendrogram"){plotDendrogram()}
        if (input$wgcnaTab=="nheatmap"){plotHtmp()}
        if (input$wgcnaTab=="ncorrelations"){plotCor()}
        dev.off()
      }
    }
  )

###################################################################################################
############################################# HEATMAP #############################################
###################################################################################################

  output$heatmapVariableSelection <- renderUI({
    sidebarPanel(
      helpText(
        paste("The heatmap and associated sample clustering are calculated with only a subset of", 
              "taxa. The taxa are ranked by the sum of abundance across samples.",sep=" ")
      ),
      sliderInput("numberHeatmapTaxa", "Number of taxa:", min=3, max=100, value=20),
      selectInput("heatmapSideColorVariable", "Side color variable:", 
                  choices = colnames(allData())),
      radioButtons("heatmapColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Categories" = "category")
      ),
      conditionalPanel(
        condition = 'input.heatmapColorType == "category"',
        numericInput("nheatmapColorCat", "Number of categories:", 4)
      ),
      HTML('<hr>'),
      HTML('<div align="right">'),
      downloadButton("saveHeatmap", "Save Plot"),
      HTML('</div>'),
      uiOutput("heatmapPlotOptions")
    )
  })

  output$heatmapPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("heatmapPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.heatmapPlotOptions == true",
        sliderInput("heatmapFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("heatmapMarCol", "Column margin", min=0.01, max=10.01, value=4.1),
        sliderInput("heatmapMarRow", "Row margin", min=0.01, max=10.01, value=1.1)
      )
    )
  })

  plotHeatmap<-function(){
    colorVariable<-which(colnames(allData())==input$heatmapSideColorVariable)
    CVlist<-getColor(allData()[,colorVariable], type=input$heatmapColorType, 
                     numCat=input$nheatmapColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    heatmapTempData<-(microbeData()[,order(apply(microbeData(), 2, sum),
                      decreasing=T)])[,1:input$numberHeatmapTaxa]
    heatmapData<-t(apply(heatmapTempData, 2, as.numeric))
    colnames(heatmapData)<-row.names(heatmapTempData)
  
    ## This is a terrible way to get colorV in the right order.
    tempclust<-hclust(dist(heatmapTempData))
    colorV<-colorV[match(tempclust$labels, row.names(microbeData()))]
    ##
    
    heatmap.2(heatmapData, scale="none", trace="none",
              lmat = cbind(c(4,2,1),c(5,3,0)), lwid=c(4,1), lhei = c(1,4,0.5), 
              Rowv=NA, dendrogram="column", ColSideColors=colorV, cexRow=input$heatmapFontSize, 
              cexCol=input$heatmapFontSize, margins=c(input$heatmapMarCol, input$heatmapMarRow)
    )
  
  }

  output$heatmapPlot <- renderPlot({
    plotHeatmap()
  })

  output$saveHeatmap <- downloadHandler(
    filename = function() { paste("heatmapPlot", fileExtension(), sep=".") },
    content = function(filename) {
      if (fileExtension()=="png"){
        png(filename, width=2000, height=2000, units="px", pointsize=25*input$heatmapFontSize)
        plotHeatmap()
        dev.off()
      }
      if (fileExtension()=="pdf"){
        pdf(filename, width=10, height=10)
        plotHeatmap()
        dev.off()
      }
    }
  )
  
###################################################################################################
######################################## STACKED BAR PLOT #########################################
###################################################################################################
  
  # dynamically generate stacked barplot UI
  output$stackedbarVariableSelection <- renderUI({  
    # generate sidebar
    sidebarPanel(
      sliderInput("numBars", "Number of taxa:", min=2, max=15, value=5),
      HTML('<br>'),
      selectInput("stackedBarOrderVariable1", "Order samples by:", 
                    choices = c("None", colnames(allData()))          ),
      selectInput("stackedBarOrderVariable2", "Secondary ordering:", 
                    choices = c("None", colnames(allData()))          ),
      selectInput("stackedBarOrderVariable3", "Tertiary ordering:", 
                    choices = c("None", colnames(allData()))          ),
      HTML('<br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveStackedbar", "Save Plot"),
      HTML('</div>'),
      uiOutput("stackedbarPlotOptions")
    )
  })

  output$stackedbarPlotOptions <- renderUI({
    mainPanel(
      checkboxInput("stackedbarPlotOptions", "Show plot options"),
      conditionalPanel(
        condition = "input.stackedbarPlotOptions == true",
        sliderInput("stackedbarFontSize", "Font size", min=0.01, max=3.01, value=1.5),
        sliderInput("stackedbarMarLeft", "Left margin", min=0.01, max=10.01, value=4.1),
        sliderInput("stackedbarMarRight", "Right margin", min=0.01, max=10.01, value=2.1),
        sliderInput("stackedbarMarTop", "Top margin", min=0.01, max=10.01, value=4.1),
        sliderInput("stackedbarMarBottom", "Bottom margin", min=0.01, max=10.01, value=5.1),
        textInput("stackedbarYlab", "Y label", value="Abundance"),
        textInput("stackedbarTitle", "Title", value=""),
        sliderInput("stackedbarKeyFontSize", "Legend font size", min=0.01, max=3.01, value=1.5),
        sliderInput("stackedbarKeyColumns", "Number of legend columns", min=1, max=15, value=3),
        textInput("stackedbarColorOrder", "Color order", 
                  value=paste(1:(input$numBars+1), collapse=","))
        
      )
    )
  })
  
  # data for stacked bar plot
  stackedData <- reactive({
      reorder<-FALSE
      sbov1<-input$stackedBarOrderVariable1
      sbov2<-input$stackedBarOrderVariable2
      sbov3<-input$stackedBarOrderVariable3
      if (sbov1!="None") sampleOrderFeature1<-allData()[,which(colnames(allData())==sbov1)]
      if (sbov2!="None") sampleOrderFeature2<-allData()[,which(colnames(allData())==sbov2)]
      if (sbov3!="None") sampleOrderFeature3<-allData()[,which(colnames(allData())==sbov3)]
      if (sbov1!="None") {sampleOrder<-order(sampleOrderFeature1, decreasing=T); reorder<-TRUE}
      if ((sbov1!="None") & (sbov2!="None")){
        sampleOrder<-order(sampleOrderFeature1, sampleOrderFeature2, decreasing=T)
        reorder<-TRUE
      }
      if ((sbov1!="None") & (sbov2!="None") & (sbov3!="None")) {
        sampleOrder<-order(sampleOrderFeature1, sampleOrderFeature2, 
                           sampleOrderFeature3, decreasing=T)
        reorder<-TRUE
      }

      #### Construct vector of breaks between order variables
      dataLength = dim(allData())[1]
      if (sbov1!="None" && (is.factor(sampleOrderFeature1) || is.character(sampleOrderFeature1))){ 
        orderedFeature1 <- as.numeric(sampleOrderFeature1[sampleOrder]) 
        breakLabels1 = paste("Order1: ",
            paste(
                levels(sampleOrderFeature1)[levels(sampleOrderFeature1)%in%sampleOrderFeature1], 
                collapse=" | " 
            )
        )
      }  
      else {
        orderedFeature1 = rep(0,dataLength)
        breakLabels1 = ""
      }
      breaks1 = orderedFeature1[1:(dataLength-1)]-orderedFeature1[2:dataLength]
     
      if (sbov2!="None" && (is.factor(sampleOrderFeature2) || is.character(sampleOrderFeature2))){ 
        orderedFeature2 <- as.numeric(sampleOrderFeature2[sampleOrder])
        breakLabels2 = paste("Order2: ",
            paste(
                levels(sampleOrderFeature2)[levels(sampleOrderFeature2)%in%sampleOrderFeature2], 
                collapse=" | " 
            )
        )
      } 
      else {
        orderedFeature2 = rep(0,dataLength)  
        breakLabels2 = ""
      }
      breaks2 = orderedFeature2[1:(dataLength-1)]-orderedFeature2[2:dataLength]
      breaks2 = breaks2 * !breaks1
      
      if (sbov3!="None" && (is.factor(sampleOrderFeature3) || is.character(sampleOrderFeature3))){ 
        orderedFeature3 <- as.numeric(sampleOrderFeature3[sampleOrder]) 
        breakLabels3 = paste("Order3: ",
            paste(
                levels(sampleOrderFeature3)[levels(sampleOrderFeature3)%in%sampleOrderFeature3], 
                collapse=" | " 
            )
        )
      }
      else {
        orderedFeature3 = rep(0,dataLength)
        breakLabels3 = ""
      }
      breaks3 = orderedFeature3[1:(dataLength-1)]-orderedFeature3[2:dataLength]
      breaks3 = breaks3 * !breaks2 * breaks1
      breaks = 10 * breaks1 + 4 * breaks2 + 1 * breaks3
      breakLabels = paste(breakLabels1,breakLabels2,breakLabels3, sep="    ")   
    
      topMicrobeCols<-order(apply(microbeData(), 2, sum), decreasing=T)
      topMicrobes<-microbeData()[,topMicrobeCols[1:input$numBars]]
      otherMicrobes<-microbeData()[,topMicrobeCols[(input$numBars+1):length(topMicrobeCols)]]
      Other<-apply(otherMicrobes, 1, sum)
      newData<-cbind(topMicrobes, Other)
      if (reorder) newData <- newData[sampleOrder,]

      # changed return value to include vector of breaks for plot spacing
      list(newData,breaks,breakLabels)  
      
  })
  
  # stacked bar plot
  plotStackedbar <- function(){

    stackedData = stackedData()[[1]]
    breaks = stackedData()[[2]]
    breakLabels = stackedData()[[3]]

    layout(matrix(c(1,2,1,2),ncol=2), height = c(4,1),width = c(4,4))
    par(mar=c(input$stackedbarMarBottom,input$stackedbarMarLeft,
              input$stackedbarMarTop,input$stackedbarMarRight))
    colorV<-discreteColors(input$numBars+1)
    colorV<-colorV[as.numeric(unlist(strsplit(input$stackedbarColorOrder,split=",")))]
    valueV<-colnames(stackedData)
    
    barplot(t(stackedData), 
         beside=F,
         space=c(0,breaks),
         col=colorV,
         ylab=input$stackedbarYlab,
         main=input$stackedbarTitle,
         border=NA, cex.axis=input$stackedbarFontSize, cex.names=input$stackedbarFontSize, 
         cex.lab=input$stackedbarFontSize, cex.main=input$stackedbarFontSize, las=3
    )

    
    plotLegend(colorV, valueV, gradient=F, cex=input$stackedbarKeyFontSize, 
               keyCol=input$stackedbarKeyColumns)
    mtext(breakLabels,side=3,line=37)

  }
  
  # save stacked bar plot
  output$saveStackedbar <- downloadHandler(
    filename = function() { paste("stackedBar", fileExtension(), sep=".") },
    content = function(filename) {
            if (fileExtension()=="png"){
              png(filename, width=2000, height=2000, units="px", 
                  pointsize=25*input$stackedbarFontSize)
              plotStackedbar()
              dev.off()
            }
            if (fileExtension()=="pdf"){
              pdf(filename, width=10, height=10)
              plotStackedbar()
              dev.off()
            }
    }
  )
  
  # display stacked bar plot of the requested variable
  output$stackedbarPlot <- renderPlot(
    plotStackedbar()
  )

###################################################################################################
########################################## PLOT OPTIONS ###########################################
###################################################################################################
  fileExtension<-reactive({
    fileExtension<-input$saveType
    if (length(fileExtension)<1) fileExtension<-"pdf"
    fileExtension
  })

  fontSize<-reactive({
    fontSize<-1.5
    fontSize
  })



})


