# ChemRICH code for Multiple Groups. Author Dinesh Kumar Barupal (dinesh.kumar@mssm.edu)

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

project_name <- "chemrich_1" # Provide this analysis a name. This will be prefixed to all the exported files.
classVariable <- "Class"
inputfile <- "20181015_KOMP_chemrich_input.xlsx"

chemrich_multi_group <- function(inputfile) {
  data_dict <- readxl::read_xlsx(inputfile, sheet="data_dict") # Data Dictionary
  pvalvec <- grep("pvalue",names(data_dict))
  classVec <- names(which(table(data_dict[[classVariable]])>2))
  clusterids <- classVec[which(classVec!="")]
  data_dict$xlogp  =  as.numeric(sapply(data_dict$SMILES, function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))

  l <- list("ChemRICH_Input" = data_dict)

  for(k in pvalvec) {
    df1 <- data.frame(Compound = data_dict$Compound_Name, Class = data_dict[[classVariable]], xlogp = data_dict$xlogp, SMILES = data_dict$SMILES, pvalue= as.numeric(data_dict[[k]]), foldchange= as.numeric(data_dict[[k+1]]), stringsAsFactors = F)

    ## Get CHEMRICH computation
    cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
      cl.member <- which(df1$Class==x)
      if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
        pval.cl.member <- df1$pvalue[cl.member]
        p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
        p.test.results$p.value
      } else {
        1
      }
    })

    cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
    clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)

    clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
      dfx <- df1[which(df1$Class==x),]
      dfx$Compound[which.min(dfx$pvalue)]
    })

    altrat <- sapply(clusterdf$name, function (k) {
      length(which(df1$Class==k & df1$pvalue<0.05))/length(which(df1$Class==k))
    })

    uprat <-sapply(clusterdf$name, function (k) {
      length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange > 1.00000000))/length(which(df1$Class==k & df1$pvalue<0.05))
    })

    clust_s_vec <- sapply(clusterdf$name, function (k) {
      length(which(df1$Class==k))
    })

    clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05))})
    clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange > 1.00000000))})
    clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(df1$Class==k & df1$pvalue<0.05 & df1$foldchange < 1.00000000))})
    clusterdf$upratio <- uprat
    clusterdf$altratio <- altrat
    clusterdf$csize <- clust_s_vec
    clusterdf <- clusterdf[which(clusterdf$csize>2),]
    clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
    clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(df1$xlogp[which(df1$Class==x)]) }))
    #clusterdf
    clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
      dfx <- df1[which(df1$Class==x),]
      paste(dfx$Compound,collapse="<br>")
    }) ## this one is the label on the tooltip of the ggplotly plot.
    clustdf <- clusterdf[which(clusterdf$pvalues!=1),]

    #################################################
    ########## Impact Visualization Graph ###########
    #################################################

    clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]
    clustdf.alt.impact <- clustdf.alt.impact[order(clustdf.alt.impact$xlogp),]
    clustdf.alt.impact$order <- 1:nrow(clustdf.alt.impact) ### Order is decided by the hclust algorithm.
    clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)

    p2 <- ggplot(clustdf.alt.impact,aes(x=xlogp,y=-log(pvalues)))
    p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
      #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
      scale_color_gradient(low = "blue", high = "red", limits=c(0,1))+
      scale_size(range = c(5, 30)) +
      scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
      scale_x_continuous(" Lipophilicity ", limit=c(min(df1$xlogp)-2,max(df1$xlogp)+2)) +
      theme_bw() +
      labs(title = paste0("ChemRICH cluster impact plot for ", gsub("_pvalue","",names(data_dict)[k]))) +
      geom_label_repel(aes(label = name), color = "gray20",family="Arial",data=subset(clustdf.alt.impact, csize>2),force = 5)+
      theme(text=element_text(family="Arial Black"))+
      theme(
        plot.title = element_text(face="bold", size=30,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
        legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        legend.spacing = unit(.05, "cm"),
        axis.text.x = element_text(size=10,angle = 0, hjust = 1),
        axis.text.y = element_text(size=15,angle = 0, hjust = 1)
      )
    read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(dml(ggobj = p2), location = ph_location(type = "body",width=10, height=8,left = 0, top = 0)) %>%
      print(target = paste0(project_name,"_",gsub("_pvalue","",names(data_dict)[k]),"_chemrich_impact_plot.pptx")) %>%
      invisible()

    ggsave(paste0(project_name,"_",gsub("_pvalue","",names(data_dict)[k]),"_chemrich_impact_plot.png"), p2,height = 8, width = 12, dpi=300)
    cat(paste0(project_name,"_",gsub("_pvalue","",names(data_dict)[k]),"_chemrich_impact_plot.pptx"," has been created.\n"))

    clustdf.e <- clusterdf[order(clusterdf$pvalues),]
    clustdf.e$pvalues <- signif(clustdf.e$pvalues, digits = 2)
    clustdf.e$adjustedpvalue <- signif(clustdf.e$adjustedpvalue, digits = 2)
    clustdf.e$upratio <- signif(clustdf.e$upratio, digits = 1)
    clustdf.e$altratio <- signif(clustdf.e$altratio, digits = 1)
    clustdf.e <- clustdf.e[,c("name","csize","pvalues","adjustedpvalue","keycpdname","alteredMetabolites","upcount","downcount","upratio","altratio")]
    names(clustdf.e) <- c("Cluster name","Cluster size","p-values","FDR","Key compound","Altered metabolites","Increased","Decreased","Increased ratio","Altered Ratio")
    l[[gsub("_pvalue","",names(data_dict)[k])]] <- clustdf.e
  }
  openxlsx::write.xlsx(l, file = paste0(project_name, "_ChemRICH_results.xlsx"), asTable = TRUE)
  cat(paste0(project_name, "_ChemRICH_results.xlsx", " has been saved.\n"))
}
