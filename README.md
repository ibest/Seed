# Seed
#### Simple Exploration of Ecological Datasets

Seed is nn R/Shiny package for visualizing taxonomic community data. It provides a visual interface for generating a wide variety of plots, including histograms, scatterplots, bar plots, stacked bar plots, PCoA plots, cluster dendrograms, and heatmaps.

## Local Installation

Seed requires R version 2.15 or later. For best results, use the latest version of R.

To install R locally, run the following commands from within R.
```
source("http://bioconductor.org/biocLite.R") 	
biocLite("impute")
biocLite("Heatplus")
install.packages(c("shiny","vegan","WGCNA","gplots","cluster"))
```

To start Seed, begin an R session and run the following commands.
```
library(shiny)
runApp("/path to Seed folder")
```

Shiny can also be downloaded and run automatically 

## Server Setup

Seed is a Shiny application that can be hosted on remote servers using the Shiny Server program. Remote hosting of Seed allows users to access Seed through any web browser. There is no need for local installation of Seed, Shiny, or R. Please find the most recent instructions for installing Shiny Server at https://github.com/rstudio/shiny-server. 


## Documentation and Support

Questions, suggestions, and bug reports are welcome and appreciated. Please contact Daniel Beck at danlbek@gmail.com or Christopher Dennis at christozoan@gmail.com.

