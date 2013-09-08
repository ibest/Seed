# Created 3/6/2013 by Daniel Beck
# Server file for shiny

library(shiny)
library(vegan)
library(WGCNA)
library(gplots)

shinyServer(function(input, output) {

#####################################################################################################
############################################# DATA ##################################################
#####################################################################################################

  # microbeData will contain relative abundances
  microbeData <- reactive({ 
    microbeData <- read.csv(input$microbeFilename$datapath) 
    if (input$relativize){
      microbeData <- t(apply(microbeData, 1, function(i) i/sum(i)))
    }
    microbeData
  })
  
  # metaData will contain all sample information other than microbial abundances
  # When reading in metadata, also calculate diversity metrics from microbeData and add to metaData
  metaData <- reactive({ cbind(read.csv(input$metaFilename$datapath),
                               shannonDiversity=diversity(microbeData(), index="shannon"),
                               simpsonDiversity=diversity(microbeData(), index="simpson"),
                               inverseSimpsonDiversity=diversity(microbeData(), index="invsimpson")) 
                         })
  
  # include all data combined in order to produce comprehensive lists of features
  allData <- reactive({ cbind(metaData(), microbeData()) })

  # extract feature names from allData
  features <- reactive({ colnames(allData) })
  
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
      checkboxInput("scatterGradientColors", "Use color gradient", FALSE),
      checkboxInput("scatterColorCat", "Categorize color variable", FALSE),
      conditionalPanel(
        condition = "input.scatterColorCat == true",
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
    sCV<-which(colnames(allData())==input$scatterColorVariable)
    fsCV<-as.numeric(as.factor(allData()[,sCV]))
    if (input$scatterColorCat){
      fsCV <- as.numeric(cut(allData()[,sCV], input$nscatterColorCat))
    }
    if (input$scatterGradientColors){
      fsCV <- colorRampPalette(c("blue", "black", "red"))(length(unique(fsCV)))[fsCV]
    }
    plot(as.numeric(allData()[,which(colnames(allData())==input$variable1)]), 
         as.numeric(allData()[,which(colnames(allData())==input$variable2)]), 
         xlab=input$variable1, 
         ylab=input$variable2, 
         main=paste("Scatterplot of", input$variable2, "vs.", input$variable1, sep=" "),
         col=fsCV
    )
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
      checkboxInput("pcaGradientColor", "Use color gradient", FALSE),
      checkboxInput("pcaColorCat", "Categorize color variable", FALSE),
      conditionalPanel(
        condition = "input.pcaColorCat == true",
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
    colorV <- as.numeric(as.factor(allData()[,which(colnames(allData())==input$pcaColorVariable)]))
    if (input$pcaColorCat){
      colorV<-as.numeric(cut(allData()[,which(colnames(allData())==input$pcaColorVariable)], input$npcaColorCat))
    }
    if (input$pcaGradientColor){
      colorV<-colorRampPalette(c("blue", "black", "red"))(length(unique(colorV)))[colorV]
    }
    principalComponents <- as.matrix(microbeData()) %*% as.matrix(eigen(cov(microbeData()))$vectors)
    # col<-rainbow(input$npcaCol+10)[1:length(input$npcaCol)]
    plot(principalComponents[,input$pcX], 
         principalComponents[,input$pcY], 
         xlab=paste("Principal component", input$pcX, sep=" "), 
         ylab=paste("Principal component", input$pcY, sep=" "), 
         main=paste("Scatter plot of principal components"),
         col=colorV
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
                  choices = c("euclidean", "maximum", "manhattan", "canberra")),
      selectInput("hclustMethod", "Cluster method:", 
                  choices = c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")),
      sliderInput("clusterCutHeight", "Subtree cut height:", min=0.0, max=1.0, value=0.5),
      numericInput("clusterGroup", "Select subtree", 1),
      selectInput("clusterColorVariable", "Color variable:", choices = colnames(allData())),
      checkboxInput("clusterGradientColor", "Use color gradient", FALSE),
      checkboxInput("clusterColorCat", "Categorize color variable", FALSE),
      conditionalPanel(
        condition = "input.clusterColorCat == true",
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
    data
  })
  
  clusterDist <- reactive({
    dist(clusterData(), method=input$distMethod)
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
    dist(subtreeData(), method=input$distMethod)
  })
  subtreeObject <- reactive({
    hclust(subtreeDist(), method=input$hclustMethod)
  })
  
  output$clusterPlot <- renderPlot({
    colorV <- as.numeric(as.factor(allData()[,which(colnames(allData())==input$clusterColorVariable)]))
    if (input$clusterColorCat){
      colorV<-as.numeric(cut(allData()[,which(colnames(allData())==input$clusterColorVariable)], input$nclusterColorCat))
    }
    if (input$clusterGradientColor){
      colorV<-colorRampPalette(c("blue", "black", "red"))(length(unique(colorV)))[colorV]
    }
    plotDendroAndColors(clusterObject(), 
                        colors=data.frame(colorV), 
                        dendroLabels=F, 
                        abHeight=input$clusterCutHeight*max(clusterObject()$height), 
                        main="Species dendrogram and module colors")  
  })
  
  output$clusterGroupPlot <- renderPlot({
    colorV <- as.numeric(as.factor(allData()[subtreeGroups()==input$clusterGroup,which(colnames(allData())==input$clusterColorVariable)]))
    if (input$clusterColorCat){
      colorV<-as.numeric(cut(allData()[subtreeGroups()==input$clusterGroup,which(colnames(allData())==input$clusterColorVariable)], input$nclusterColorCat))
    }
    if (input$clusterGradientColor){
      colorV<-colorRampPalette(c("blue", "black", "red"))(length(unique(colorV)))[colorV]
    }
    plot(subtreeObject())
  })

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
    ng<-length(unique(groups))
    moduleColors<-(rainbow(ng+3)[1:ng])[groups+1]
    plotDendroAndColors(hdADJ(), 
                        colors=data.frame(moduleColors), 
                        dendroLabels=F, 
                        abHeight=input$cutLevel, 
                        main="Species dendrogram and module colors")
  }
  output$dendroPlot <- renderPlot({
    plotDendrogram()
  })
  
  plotHeatmap <- function(){
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
    plotHeatmap()
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
        if (input$wgcnaTab=="nheatmap"){plotHeatmap()}
        if (input$wgcnaTab=="ncorrelations"){plotCor()}
      dev.off()
    }
  )
})
