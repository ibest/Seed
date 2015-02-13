# Seed
#### Simple Exploration of Ecological Datasets

Seed is an R/Shiny package for visualizing ecological data. It provides a visual interface for generating a wide variety of plots, including histograms, scatterplots, bar plots, stacked bar plots, PCoA plots, cluster dendrograms, and heatmaps.

## Availability

Seed is available at https://github.com/danlbek/Seed.
Direct download link: https://github.com/danlbek/Seed/archive/master.zip

## Local Installation

Seed requires R version 2.15 or later. For best results, use the latest version of R.

Seed depends on several R packages. To install them, run the following commands from within R.
```r
source("http://bioconductor.org/biocLite.R") 	
biocLite("impute")
biocLite("Heatplus")
biocLite("preprocessCore")
biocLite("GO.db")
install.packages(c("shiny","vegan","WGCNA","gplots","cluster"))
```

To start Seed, begin an R session and run the following commands.
```r
library(shiny)
runApp("/path to Seed folder")
```

After installing Seed dependencies, the latest stable version of Seed can also be downloaded and run simultaneously using the following R commands.
```r
library(shiny)
runGitHub("Seed","danlbek")
```

## Server Setup

Seed is a Shiny application that can be hosted on remote servers using the Shiny Server program. Remote hosting of Seed allows users to access Seed through any web browser. There is no need for local installation of Seed, Shiny, or R. Please find the most recent instructions for installing Shiny Server at https://github.com/rstudio/shiny-server. 


## Documentation and Support
Help files may be accessed from within the Seed interface using the "Help" tab. A detailed walkthrough of an example dataset (walkthrough.pdf) is also included with the Seed source files. 

Questions, suggestions, and bug reports are welcome and appreciated. Please contact Daniel Beck at danlbek@gmail.com or Christopher Dennis at christozoan@gmail.com.

