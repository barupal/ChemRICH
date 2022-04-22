# This script runs the ChemRICH analysis for a minimum input.
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

run_chemrich_basic <- function(inputfile = "name_of_the_xlsx_file") {
    ndf <- data.frame(readxl::read_xlsx(path = inputfile, sheet = 1), stringsAsFactors = F)
    ndf$edirection <- "up"
    ndf$efs <- 1

    if(length(which(ndf$effect_size < 0)) >0) { # if regression models
      ndf$edirection[which(ndf$effect_size < 0)] <- "down"
      ndf$edirection[which(ndf$pvalue > 0.1)] <- "no change"
      ndf$efs[which(ndf$effect_size < 0)] <- 1/abs(ndf$effect_size[which(ndf$effect_size < 0)])
      ndf$efs[which(ndf$effect_size > 1)] <- abs(ndf$effect_size[which(ndf$effect_size > 1)])
      ndf$efs [which(ndf$pvalue > 0.1)] <- 1

    } else { # if student test
      ndf$edirection[which(ndf$effect_size < 1)] <- "down"
      ndf$edirection[which(ndf$pvalue > 0.1)] <- "no change"
      ndf$efs[which(ndf$effect_size < 1)] <- 1/ndf$effect_size[which(ndf$effect_size < 1)]
      ndf$efs[which(ndf$effect_size > 1)] <- ndf$effect_size[which(ndf$effect_size > 1)]
      ndf$efs [which(ndf$pvalue > 0.1)] <- 1
    }

    clusterids <- names(which(table(ndf$set)>2))
    clusterids <- clusterids[which(clusterids!="")]
    cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
      cl.member <- which(ndf$set==x)
      if( length(which(ndf$pvalue[cl.member]<.10)) >0 ){
        pval.cl.member <- ndf$pvalue[cl.member]
        p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
        p.test.results$p.value
      } else {
        1
      }
    })

    cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
    #clusterdf <- data.frame(name=clusterids[which(cluster.pvalues!=10)],pvalues=cluster.pvalues[which(cluster.pvalues!=10)], stringsAsFactors = F)
    clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)

    clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
      dfx <- ndf[which(ndf$set==x),]
      dfx$compound_name[which.min(dfx$pvalue)]
    })

    altrat <- sapply(clusterdf$name, function (k) {
      length(which(ndf$set==k & ndf$pvalue<0.1))/length(which(ndf$set==k))
    })

    uprat <-sapply(clusterdf$name, function (k) {
      length(which(ndf$set==k & ndf$pvalue<0.1 & ndf$edirection == "up"))/length(which(ndf$set==k & ndf$pvalue<0.10))
    })

    clust_s_vec <- sapply(clusterdf$name, function (k) {
      length(which(ndf$set==k))
    })

    clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.05))})

    clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.05 & ndf$edirection == "up"))})

    clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(ndf$set==k & ndf$pvalue<0.05 & ndf$edirection == "down"))})

    clusterdf$upratio <- uprat
    clusterdf$altratio <- altrat
    clusterdf$csize <- clust_s_vec
    clusterdf <- clusterdf[which(clusterdf$csize>2),]
    clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")

    clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(ndf$order[which(ndf$set==x)]) }))

    clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
      dfx <- ndf[which(ndf$set==x),]
      paste(dfx$compound_name,collapse="<br>")
    })

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
      scale_x_continuous(" chemical set order (as provided) ") +
      theme_bw() +
      #labs(title = "ChemRICH cluster impact plot") +
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

    p2

    read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(dml(ggobj = p2), location = ph_location(type = "body",width=10, height=8,left = 0, top = 0)) %>%
      print(target = paste0("chemrich_class_impact_plot.pptx")) %>%
      invisible()

    ggsave(paste0("chemrich_class_impact_plot.png"), p2,height = 8, width = 12, dpi=300)

    cat(paste0("chemrich_class_impact_plot.pptx"," has been created.\n"))

    ## Export the result table.
    clustdf.e <- clusterdf[order(clusterdf$pvalues),]
    clustdf.e$pvalues <- signif(clustdf.e$pvalues, digits = 2)
    clustdf.e$adjustedpvalue <- signif(clustdf.e$adjustedpvalue, digits = 2)
    clustdf.e$upratio <- signif(clustdf.e$upratio, digits = 1)
    clustdf.e$altratio <- signif(clustdf.e$altratio, digits = 1)
    clustdf.e <- clustdf.e[,c("name","csize","pvalues","adjustedpvalue","keycpdname","alteredMetabolites","upcount","downcount","upratio","altratio")]
    names(clustdf.e) <- c("Cluster name","Cluster size","p-values","FDR","Key compound","Altered metabolites","Increased","Decreased","Increased ratio","Altered Ratio")
    #df1$TreeLabels <- treeLabels
    ndf$pvalue <- signif(ndf$pvalue, digits = 2)
    ndf$efs <- signif(ndf$efs, digits = 2)
    ndf$FDR <- signif(  p.adjust(ndf$pvalue), digits = 2)
    l <- list("ChemRICH_Results" = clustdf.e, "Compound_ChemRICH" = ndf )
    openxlsx::write.xlsx(l, file = paste0("chemRICH_class_results.xlsx"), asTable = TRUE)
    cat(paste0("chemRICH_class_results.xlsx", " has been saved.\n"))

    ##################################
    #### Interactive Cluster Plot ####
    ##################################

    p2 <- ggplot(clustdf.alt.impact,aes(label=name,label2=pvalues, label3=csize,label4=Compounds))
    p2 <- p2 + geom_point(aes(x=xlogp,y=-log(pvalues),size=csize, color=upratio)) +
      #labs(caption = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
      scale_color_gradient(low = "blue", high = "red", limits=c(0,1))+
      scale_size(range = c(5, 30)) +
      scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+5  )) +
      scale_x_continuous("chemical set order (as provided) ") +
      theme_bw() +
      #labs(title = "ChemRICH cluster impact plot") +
      geom_text(aes(x=xlogp,y=-log(pvalues),label = name), color = "gray20",data=subset(clustdf.alt.impact, csize>2))+
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
        axis.text.x = element_text(size=15,angle = 0, hjust = 1),
        axis.text.y = element_text(size=15,angle = 0, hjust = 1)
      )
    gg <- ggplotly(p2,tooltip = c("label","label2","label4"), width = 1600, height = 1000)
    gg
    saveWidget(gg,file = paste0("chemrich_interactive.html"), selfcontained = T)
    cat(paste0("chemrich_interactive.html", " has been saved.\n"))

    ##############################
    ## Interactive volcano Plot this is the pending code ##
    ##############################

    # # we need to add the interactive volcano plot for the p-values and fold-change sorted by the MeSH tree order.
    # df2 <- chemrich.input.file
    # df2 <- df2[order(sapply(df2$ChemRICHClusters, function(x) { median(df2$xlogp[which(df2$ChemRICHClusters==x)])     }  )),]
    # df2$Changed <- "No Change"
    # df2$Changed[which(df2$pvalue<0.05 & df2$foldchange>1)] <- "UP"
    # df2$Changed[which(df2$pvalue<0.05 & df2$foldchange<1)] <- "DOWN"
    # df2$Changed <- as.factor(df2$Changed)
    # df2$pathway <- "No"
    # df2$pathway[which(df2$PubChemID%in%cid_biosys==TRUE)] <- "yes"
    # df2$pathway <- as.factor(df2$pathway)
    # df2$Compound.Name <- factor(df2$CompoundName, levels =df2$CompoundName)
    # df2$foldchange <- round(sapply(df2$foldchange, function(x) { if(x>1) {x} else {1/x} }), digits = 1)
    # df2$foldchange[ df2$foldchange>5] <- 5
    #
    # p2 <-   ggplot(df2, aes(label=CompoundName,x=CompoundName, y=-log(pvalue,base = 10),colour = Changed,shape=pathway, size=foldchange)) +  scale_size(range = c(1, 10)) +
    #   #geom_line(position=pd, size=2)+
    #   #geom_errorbar(aes(ymin = V2-V3 , ymax=V2+V3), width=.3,size=2,position=pd) +
    #   geom_point(stat = "identity") + # 21 is filled circle
    #   #geom_bar(stat="identity", size=.1,position=position_dodge()) +
    #   scale_y_continuous("pvalue (-log)") +
    #   scale_x_discrete("Metabolites: Red- increased,blue-decreased,yellow-not significant, solid-pathway(s) found ") +
    #   scale_color_manual("Student ttest",values=c("blue", "yellow", "red","white")) +
    #   scale_fill_manual("",values=c("white", "yellow", "red","white")) +
    #   scale_shape_manual("Pathway found",values=c(1,16))+
    #   #scale_shape(solid = FALSE) +
    #   theme_bw() +
    #   labs(title = "Metabolic Dys-regulations sorted by chemical similarity") +
    #   theme(
    #     plot.title = element_text(face="bold", size=30,hjust = 0.5),
    #     axis.title.x = element_text(face="bold", size=20),
    #     axis.title.y = element_text(face="bold", size=30, angle=90),
    #     panel.grid.major = element_blank(), # switch off major gridlines
    #     panel.grid.minor = element_blank(), # switch off minor gridlines
    #     #legend.justification=c(1,0),
    #     #legend.position=c(1,.6),
    #     legend.position = "none",
    #     #legend.title = element_blank(), # switch off the legend title
    #     #legend.text = element_text(size=12),
    #     #legend.key.size = unit(1.5, "lines"),
    #     #legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    #     #legend.spacing = unit(.05, "cm"),
    #     #axis.text.x = element_text(size=15,angle = 45, hjust = 1.0),
    #     axis.text.x= element_blank(),
    #     axis.text.y = element_text(size=15,angle = 0, hjust = 0.5)
    #   )
    # p3 <- ggplotly(p2, width = 1600, height = 1000)
    # saveWidget(gg,file = paste0(project_name, "_chemrich_dysregulation.html"), selfcontained = T)
    # cat("Interactive plots have been generated.")
}
