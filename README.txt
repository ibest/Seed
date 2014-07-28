Seed - An R/Shiny package for visualizing taxonomic community data.

Created by:  Daniel Beck
Contact: danlbek@gmail.com

Maintainer: Christopher Dennis
Contact: christozoan@gmail.com

License: GPLv3

Requirements: R version 2.15 or later

Installation: From within R run the following commands.
	source("http://bioconductor.org/biocLite.R") 
	biocLite("impute")
	biocLite("Heatplus")
	install.packages(c("shiny","vegan","WGCNA","gplots","cluster","devtools"))
	devtools::install_github("shiny-incubator","rstudio")
	
Starting Seed: From within R run the following commands.
	library(shiny)
	runApp("/path to Seed folder")

