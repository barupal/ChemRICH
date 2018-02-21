# ChemRICH
R package for Chemical Similarity Enrichment Analysis

ChemRICH installation method

In R, run the following code.
```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")
if (!require("opencpu"))
install.packages('opencpu', repos="http://cran.rstudio.com/")
if (!require("RCurl"))
install.packages('RCurl', repos="http://cran.rstudio.com/")
if (!require("pacman"))
install.packages('pacman', repos="http://cran.rstudio.com/")
library(devtools)
library(RCurl)
library(pacman)
source('https://bioconductor.org/biocLite.R')
pacman::p_load(grid,rcdk, RJSONIO,RCurl, dynamicTreeCut,ape,ggplot2,ggrepel,ReporteRs,XLConnect,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,XLConnect)
install_github('barupal/chemrich')
library(ChemRICH)
library(opencpu)
opencpu::ocpu_start_server()
```
Then go to :
[ChemRICH Local Version](http://localhost:5656/ocpu/library/ChemRICH/www/)

## Online version 

 [ChemRICH Online Version](http://chemrich.fiehnlab.ucdavis.edu)

## Citation

[Barupal, D.K. and Fiehn, O., 2017. Chemical Similarity Enrichment Analysis (ChemRICH) as alternative to biochemical pathway mapping for metabolomic datasets Scientific Report 2017. (https://www.nature.com/articles/s41598-017-15231-w)

## Docker image 

ChemRICH docker image (https://hub.docker.com/r/barupal/chemrich-docker/)
