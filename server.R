# Created 3/6/2013 by Daniel Beck
# Server file for shiny

library(shiny)
library(vegan)
library(WGCNA)
library(gplots)

shinyServer(function(input, output) {

#####################################################################################################
########################################## FUNCTIONS ################################################
#####################################################################################################

  # functions used by multiple plots

  # color function wrappers
  gradientColors<-function(n){
    colorRampPalette(c("blue", "red"))(n)
  }
  discreteColors<-function(n){
    trim<-as.integer(n*0.1)+3   # prevents color overlap at ends of rainbow
    rainbow(n+trim)[1:n]
  }
  
  # generate color vector for plots. Returns both color vector and value vector (for use in legends).
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
  plotLegend<-function(colorVector, valueVector, gradient=F, title="", min=0, max=1){
    if (gradient){
      par(mar=c(6,4,2,2)+0.1)
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = title)
      legend_image <- as.raster(matrix(gradientColors(100), ncol=100))
      rasterImage(legend_image, 0, 0, 1, 1)
      axis(side=1, at=seq(0,1,l=5), labels=seq(min,max,l=5),col.axis="black")
      mtext("Key", side=2, las=2)
      
    }else{
      par(mar=c(0,2,0,2)+0.1)
      plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = "")
      uniquePairs<-match(levels(as.factor(valueVector)), valueVector)
      values<-levels(as.factor(valueVector))
      colors<-colorVector[uniquePairs]
      legend("center", legend=values, fill=colors, 
             ncol=min(length(colors), 8), 
             box.col="white", title=title)
    }
  }

  # tests to see if file is square (if it has column names for every column)
  isSquare<-function(filename, sep){
    test<-sapply(sapply(readLines(filename,2), function(i) strsplit(i,sep)), length)
    return(test[1]==test[2])
  }
#####################################################################################################
############################################# DATA ##################################################
#####################################################################################################

  # microbeData will contain relative abundances
  microbeData <- reactive({ 
    microbeFile <- input$microbeFilename$datapath
    if (is.null(microbeFile)) return(NULL)
    microbeData <- read.csv(microbeFile, header=input$microbeHeader, sep=input$microbeSep, quote=input$microbeQuote) 
    if (isSquare(microbeFile, input$microbeSep)){
      row.names(microbeData)<-microbeData[,1]
      microbeData<-microbeData[,-1]
    }
    if (input$relativize){
      microbeData <- t(apply(microbeData, 1, function(i) i/sum(i)))
    }
    if (input$presenceAbsence){
      rn<-row.names(microbeData)
      microbeData<-apply(microbeData, 2, function(i) as.numeric(as.logical(i)))
      row.names(microbeData)<-rn
    }
    microbeData
  })
  
  # metaData will contain all sample information other than microbial abundances
  # When reading in metadata, also calculate diversity metrics from microbeData and add to metaData
  metaData <- reactive({ 
    metaFile <- input$metaFilename$datapath
    if (is.null(metaFile)) return(NULL)
    metaData <- read.csv(metaFile, header=input$metaHeader, sep=input$metaSep, quote=input$metaQuote)
    if (isSquare(metaFile, input$metaSep)){
      row.names(metaData)<-metaData[,1]
      metaData<-metaData[,-1]
    }
    if (is.null(input$microbeFilename$datapath)) return(metaData)
    cbind(metaData,
      shannonDiversity=diversity(microbeData(), index="shannon"),
      simpsonDiversity=diversity(microbeData(), index="simpson"),
      inverseSimpsonDiversity=diversity(microbeData(), index="invsimpson")
    ) 
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

#####################################################################################################
########################################### HISTOGRAM ###############################################
#####################################################################################################

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
      HTML('</div>')
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
    hist(as.numeric(allData()[,which(colnames(allData())==input$variable)]), 
         breaks=input$breaks, 
         xlab=input$variable, 
         main=paste("Histogram of", input$variable, sep=" "))
  }
  
  # save histogram plot
  output$saveHist <- downloadHandler(
    filename = function() { paste("hist", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
        plotHistogram()
      dev.off()
    }
  )
  
  # display histogram of the requested variable
  output$histPlot <- renderPlot({
    plotHistogram()
  })
  
#####################################################################################################
########################################### SCATTER PLOT ############################################
#####################################################################################################

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
                       "Cagetories" = "category")
      ),
      conditionalPanel(
        condition = 'input.scatterColorType == "category"',
        numericInput("nscatterColorCat", "Number of categories:", 4)
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveScatter", "Save Plot"),
      HTML('</div>')
    )
  })
  
  # generate scatter plot
  plotScatter<-function(){
    layout(matrix(c(1,2,1,2),ncol=2), height = c(4,1),width = c(4,4))
    
    colorVariable<-which(colnames(allData())==input$scatterColorVariable)
    CVlist <- getColor(allData()[,colorVariable], type=input$scatterColorType, numCat=input$nscatterColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    plot(as.numeric(allData()[,which(colnames(allData())==input$variable1)]), 
         as.numeric(allData()[,which(colnames(allData())==input$variable2)]), 
         xlab=input$variable1, 
         ylab=input$variable2, 
         main=paste("Scatterplot of", input$variable2, "vs.", input$variable1, sep=" "),
         col=colorV
    )
    plotLegend(colorV, valueV, gradient=(input$scatterColorType=="gradient"), 
               title=input$scatterColorVariable, min=min(as.numeric(valueV), na.rm=T), max=max(as.numeric(valueV), na.rm=T))
  }
  
  # save scatter plot
  output$saveScatter <- downloadHandler(
    filename = function() { paste("scatter", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
      plotScatter()
      dev.off()
    }
  )
  
  # Generate scatter plot of the requested variables
  output$scatterPlot <- renderPlot({
    plotScatter()
  })
  
#####################################################################################################
############################################## PCA ##################################################
#####################################################################################################

  # generate PCA UI
  output$pcaVariableSelection <- renderUI({
    sidebarPanel(
      numericInput("pcX", "Principal component X:", 1),
      numericInput("pcY", "Principal component Y:", 2), 
      selectInput("pcaColorVariable", "Color variable:", choices = colnames(allData())),
      radioButtons("pcaColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Cagetories" = "category")
      ),
      conditionalPanel(
        condition = 'input.pcaColorType == "category"',
        numericInput("npcaColorCat", "Number of categories:", 4)
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("savePca", "Save Plot"),
      HTML('</div>')
    )
  })

  # generate PCA plot
  plotPca <- function(){
    layout(matrix(c(1,2,1,2),ncol=2), height = c(4,1),width = c(4,4))
    
    colorVariable<-which(colnames(allData())==input$pcaColorVariable)
    CVlist<-getColor(allData()[,colorVariable], type=input$pcaColorType, numCat=input$npcaColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    principalComponents <- as.matrix(microbeData()) %*% as.matrix(eigen(cov(microbeData()))$vectors)
    plot(principalComponents[,input$pcX], 
         principalComponents[,input$pcY], 
         xlab=paste("Principal component", input$pcX, sep=" "), 
         ylab=paste("Principal component", input$pcY, sep=" "), 
         main=paste("Scatter plot of principal components"),
         col=colorV
    )
    plotLegend(colorV, valueV, gradient=(input$pcaColorType=="gradient"), 
               title=input$pcaColorVariable, min=min(as.numeric(valueV), na.rm=T), max=max(as.numeric(valueV), na.rm=T)
    )
  }

  # sacce PCA plot
  output$savePca <- downloadHandler(
    filename = function() { paste("PCA", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
      plotPca()
      dev.off()
    }
  )
  
  # display PCA plot
  output$pcaPlot <- renderPlot({
    plotPca()
  })
  
#####################################################################################################
############################################ BAR PLOT ###############################################
#####################################################################################################

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
      HTML('</div>')
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
    x<-barplot(means, 
               xlab=colnames(allData())[catVarCol], 
               ylab=colnames(allData())[barVarCol], 
               ylim=c(min(means-errors), max(means+errors)),
               xpd=F,
               col="#f5f5f5")
    arrows(x,means+errors,x,means-errors,code=0)
    text(x,min(means-errors)+(max(means+errors)-min(means-errors))/20, ns, col="blue")
  }

  output$saveBar <- downloadHandler(
    filename = function() { paste("barPlot", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
      plotBar()
      dev.off()
    }
  )
  
  # Generate barplot of the requested variable
  output$barPlot <- renderPlot({
    plotBar()
  })
  
#####################################################################################################
############################################# CLUSTER ###############################################
#####################################################################################################

  # dynamic cluster UI
  output$clusterVariableSelection <- renderUI({
    sidebarPanel(
      selectInput("distMethod", "Distance method:", 
                  choices = c("euclidean", "manhattan", "canberra", "bray", "kulczynski", "jaccard", 
                              "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao")),
      selectInput("hclustMethod", "Cluster method:", 
                  choices = c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")),
      sliderInput("clusterCutHeight", "Subtree cut height:", min=0.0, max=1.0, value=0.5),
      numericInput("clusterGroup", "Select subtree", 1),
      selectInput("clusterColorVariable", "Color variable:", choices = colnames(allData())),
      radioButtons("clusterColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Cagetories" = "category")
      ),
      conditionalPanel(
        condition = 'input.clusterColorType == "category"',
        numericInput("nclusterColorCat", "Number of categories:", 4)
      ),
      radioButtons("clusterChoice", "Select features that define samples", 
                   choices=c("Metadata", "Microbe data", "All data", "Custom"), selected="Microbe data"),
      conditionalPanel(
        condition = "input.clusterChoice == 'Custom'",
        checkboxGroupInput(inputId="customClusterVariables", label="", choices=colnames(allData()))
      ),
      HTML('<br><br><br>'),
      HTML('<div align="right">'),
      downloadButton("saveCluster", "Save Plot"),
      HTML('</div>')
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
    CVlist<-getColor(allData()[,colorVariable], type=input$clusterColorType, numCat=input$nclusterColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    layout(matrix(c(1,2,3,1,2,3),ncol=2), height = c(4,1,1),width = c(4,4))
    plotDendroAndColors(clusterObject(), 
                        colors=data.frame(colorV), 
                        dendroLabels=F, 
                        abHeight=input$clusterCutHeight*max(clusterObject()$height), 
                        groupLabels=input$clusterColorVariable,
                        main="",
                        setLayout=FALSE)
    plotLegend(colorV, valueV, gradient=(input$clusterColorType=="gradient"), 
               title=input$clusterColorVariable, min=min(as.numeric(valueV), na.rm=T), max=max(as.numeric(valueV), na.rm=T))
  }

  plotSubTree<-function(){
    colorVariable<-which(colnames(allData())==input$clusterColorVariable)
    CVlist<-getColor(allData()[subtreeGroups()==input$clusterGroup,colorVariable], 
                     type=input$clusterColorType, numCat=input$nclusterColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    layout(matrix(c(1,2,3,1,2,3),ncol=2), height = c(4,1,1),width = c(4,4))
    plotDendroAndColors(subtreeObject(), 
                        colors=data.frame(colorV), 
                        dendroLabels=F, 
                        groupLabels=input$clusterColorVariable,
                        main="",
                        setLayout=FALSE)  
    plotLegend(colorV, valueV, gradient=(input$clusterColorType=="gradient"), 
               title=input$clusterColorVariable, min=min(as.numeric(valueV), na.rm=T), max=max(as.numeric(valueV), na.rm=T))
  }

  output$clusterPlot <- renderPlot({
    plotCompleteTree()
  })

  output$clusterGroupPlot <- renderPlot({
    plotSubTree()
  })

  output$saveCluster <- downloadHandler(
    filename = function() { paste("clusterPlot", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
      if (input$clusterTab=="complete"){plotCompleteTree()}
      if (input$clusterTab=="subtree"){plotSubTree()}
      dev.off()
    }
  )

#####################################################################################################
############################################## WGCNA ################################################
#####################################################################################################

  # render sidebars for wgcna plots
  output$wgcnaVariableSelection <- renderUI({
    sidebarPanel(
      helpText("Warning: The Kendall correlation method is very slow."),
      selectInput("corMethod", "Network correlation method:", choices = c("pearson", "spearman", "kendall")),
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
    groups<-cutreeDynamic(dendro=hdADJ(), distM=dADJ(), minClusterSize=3,method="tree",cutHeight=input$cutLevel)
    moduleColors<-getColor(groups, "unique")[[1]]
    plotDendroAndColors(hdADJ(), 
                        colors=data.frame(moduleColors), 
                        dendroLabels=F, 
                        abHeight=input$cutLevel, 
                        main="Species dendrogram and module colors")
  }
  output$dendroPlot <- renderPlot({
    plotDendrogram()
  })
  
  plotHtmp <- function(){
    groups<-cutreeDynamic(dendro=hdADJ(), distM=dADJ(), minClusterSize=3,method="tree",cutHeight=input$cutLevel)
    corlist<-sapply(unique(groups), function(i) cors()[groups==i, groups==i])
    cc<-colorRampPalette(c("white", "blue"))
    heatmap.2(as.matrix(corlist[[input$selectGroup]]), 
              scale="none", 
              trace="none", 
              margins=c(12,12), 
              cexRow=0.9, 
              cexCol=0.9,
              col=cc)
  }

  output$htmpPlot <- renderPlot({
    plotHtmp()
  })
  
  plotCor <- function(){
    groups<-cutreeDynamic(dendro=hdADJ(), distM=dADJ(), minClusterSize=3,method="tree",cutHeight=input$cutLevel)
    MEs0 = moduleEigengenes(microbeData(), groups+1)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, metaData(), use="p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(metaData()));
    
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
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
                   cex.lab.x=0.8)    
  }

  output$corPlot <- renderPlot({
    plotCor()
  })
  
  output$saveWGCNA <- downloadHandler(
    filename = function() { paste("wcgnaPlot", '.pdf', sep='') },
    content = function(filename) {
      pdf(filename, width=10, height=10)
        if (input$wgcnaTab=="ndendrogram"){plotDendrogram()}
        if (input$wgcnaTab=="nheatmap"){plotHtmp()}
        if (input$wgcnaTab=="ncorrelations"){plotCor()}
      dev.off()
    }
  )

#####################################################################################################
############################################# HEATMAP ###############################################
#####################################################################################################

output$heatmapVariableSelection <- renderUI({
    sidebarPanel(
      helpText("The heatmap and associated sample clustering are calculated with only a subset of taxa. 
               The taxa are ranked by the sum of abundance across samples."),
      sliderInput("numberHeatmapTaxa", "Number of taxa:", min=3, max=100, value=20),
      selectInput("heatmapSideColorVariable", "Side color variable:", choices = colnames(allData())),
      radioButtons("heatmapColorType", "Color options:", 
                   list("Unique" = "unique",
                        "Gradient" = "gradient",
                        "Cagetories" = "category")
      ),
      conditionalPanel(
        condition = 'input.heatmapColorType == "category"',
        numericInput("nheatmapColorCat", "Number of categories:", 4)
      ),
      HTML('<hr>'),
      HTML('<div align="right">'),
      downloadButton("saveHeatmap", "Save Plot"),
      HTML('</div>')
    )
  })

  plotHeatmap<-function(){
    colorVariable<-which(colnames(allData())==input$heatmapSideColorVariable)
    CVlist<-getColor(allData()[,colorVariable], type=input$heatmapColorType, numCat=input$nheatmapColorCat)
    colorV <- CVlist[[1]]
    valueV <- CVlist[[2]]
    heatmapTempData<-(microbeData()[,order(apply(microbeData(), 2, sum), decreasing=T)])[,1:input$numberHeatmapTaxa]
    heatmapData<-t(apply(heatmapTempData, 2, as.numeric))
    colnames(heatmapData)<-row.names(heatmapTempData)
  
    ## This is a terrible way to get colorV in the right order.
    tempclust<-hclust(dist(heatmapTempData))
    colorV<-colorV[match(tempclust$labels, row.names(microbeData()))]
    ##
    
    heatmap.2(heatmapData, scale="none", trace="none",
              lmat = cbind(c(4,2,1),c(5,3,0)), lwid=c(4,1), lhei = c(1,4,0.5), 
              Rowv=NA, dendrogram="column", ColSideColors=colorV)
  
  }

  output$heatmapPlot <- renderPlot({
    plotHeatmap()
  })

  output$saveHeatmap <- downloadHandler(
    filename = function() { paste("heatmapPlot", ".pdf", sep="")},
    content = function(filename) {
      pdf(filename, width=10, height=10)
      plotHeatmap()
      dev.off()
    }
  )
  


})
