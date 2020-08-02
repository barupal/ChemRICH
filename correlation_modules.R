# This script predicts the correlation modules for a metabolomics data-set.
# Author Dinesh Kumar Barupal (dinesh.kumar@mssm.edu) August 2020.

load.ChemRICH.Packages <- function() {
  if (!require("devtools"))
    install.packages('devtools', repos="http://cran.rstudio.com/")
  if (!require("RCurl"))
    install.packages('RCurl', repos="http://cran.rstudio.com/")
  if (!require("pacman"))
    install.packages('pacman', repos="http://cran.rstudio.com/")
  library(devtools)
  library(RCurl)
  library(pacman)
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  pacman::p_load(GGally)
  pacman::p_load(DT)
  pacman::p_load(RCurl)
  pacman::p_load(RJSONIO)
  pacman::p_load(ape)
  pacman::p_load(devEMF)
  pacman::p_load(dynamicTreeCut)
  pacman::p_load(extrafont)
  pacman::p_load(ggplot2)
  pacman::p_load(ggpubr)
  pacman::p_load(ggrepel)
  pacman::p_load(grid)
  pacman::p_load(htmlwidgets)
  pacman::p_load(igraph)
  pacman::p_load(magrittr)
  pacman::p_load(network)
  pacman::p_load(officer)
  pacman::p_load(openxlsx)
  pacman::p_load(phytools)
  pacman::p_load(plotly)
  pacman::p_load(plotrix)
  pacman::p_load(rcdk)
  pacman::p_load(readxl)
  pacman::p_load(rvg)
  pacman::p_load(sna)
  pacman::p_load(visNetwork)
}

chemrich_predict_correlation_modules <- function(inputfile = "nameofthefile") {

  ndf <- data.frame(readxl::read_xlsx(path = inputfile, sheet = 1,col_names = T), stringsAsFactors = F)
  data_matrix.sb <- ndf[,-1]
  data_matrix.sb <- do.call(rbind,lapply(1:nrow(data_matrix.sb),function(x) {as.numeric(data_matrix.sb[x,])}))
  cormat <- cor(t(data_matrix.sb), method = "spearman")
  s <- cormat
  diag(s) <- 0
  s[is.na(s)] <- 0
  hc <- hclust(as.dist(1-s), method="ward.D2") # ward method provide better clusters.
  clust1 <- cutreeDynamic(hc,distM = as.matrix(1-s),deepSplit =4, minClusterSize = 3)
  clusternames <- paste0("clust_",clust1)

  xdf <- data.frame(compound_name = ndf$CompoundID, order = sapply(ndf$CompoundID, function(x) {which(ndf$CompoundID[hc$order]==x)}), pvalue = 0,effect_size=0,set=clusternames, stringsAsFactors = F)

  l <- list("raw" = ndf, "chemrich_input" = xdf )
  openxlsx::write.xlsx(l, file = paste0("correlation_modules_prediction.xlsx"), asTable = TRUE)
  cat(paste0("correlation_modules_prediction.xlsx", " has been saved.\n Voilla!!"))

}





















