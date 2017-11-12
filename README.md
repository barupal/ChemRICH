# ChemRICH
R package for Chemical Similarity Enrichment Analysis

ChemRICH installation method

In R, run following code.
```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")
if (!require("opencpu"))
install.packages('opencpu', repos="http://cran.rstudio.com/")
if (!require("RCurl"))
install.packages('RCurl', repos="http://cran.rstudio.com/")
library(devtools)
library(RCurl)
source('https://bioconductor.org/biocLite.R')
install_github('barupal/chemrich')
library(chemrich)
library(opencpu)
opencpu$browse('/library/ChemRICH/www')
```
## Online version 

 [ChemRICH Online Version](http://chemrich.fiehnlab.ucdavis.edu)

## Citation

[Barupal, D.K. and Fiehn, O., 2017. Chemical Similarity Enrichment Analysis (ChemRICH) as alternative to biochemical pathway mapping for metabolomic datasets Scientific Report 2017. (https://www.nature.com/articles/s41598-017-15231-w)

## Docker image 

ChemRICH docker image (https://hub.docker.com/r/barupal/chemrich-docker/)
