Started on 3/6/2013 by Daniel Beck

Project goal:
  Develop an application for viewing microbial community datasets

Strategy:
  My initial approach is to use Shiny as an interface
  The microbial data will be already be processed from sequences into abundances
  Rows are samples, columns are microbes and metadata
  vegan and associated libraries can be used for ecological analyses

Temporary install notes:
  R version 2.15 or later is required
  Install necessary libraries:
              source("http://bioconductor.org/biocLite.R") 
              biocLite("impute")
              install.packages(c("shiny","vegan","WGCNA","qgraph","gplots"))

  To run: 
              library(shiny)
              runApp("/path_to_microbePlot_folder")

As of 10/6/2013, microbePlot can be used without local installation by visiting http://bioinfo-mite.ibest.uidaho.edu:3838/microbePlot/
However, the server is currently only available from within the UI network. 

Notes:
  Input data must be in two files. The metadata file should include sample information. The microbe data file should include the relative abundances of the microbes in each sample.
  The files must be in CSV format with equal numbers of samples.
  Samples must be in rows, with the same order in both files.

TODO:
  Future updates should allow more flexibility for input data.
  Does it make sense to scale any of the data before an analysis? i.e. PCA or clustering. Microbe data as relative abundance is already scaled.
