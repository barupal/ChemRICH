# This script predicts the MeSH chemical class for a given compound list having SMILES codes.
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

getCNames <- function(x) {
  if(length(which(treenames.df$MESHTREE==x))>0){
    as.character(treenames.df$ClassName[which(treenames.df$MESHTREE==x)])
  } else {
    x
  }
}

makeSmiles.clean <- function (smi) {
  smi <- gsub("@","",smi)
  smi <- gsub("O-","O",smi)
  smi <- gsub("H[1-4][+]","",smi)
  smi <- gsub("N[+]","N",smi)
  smi <- gsub("/|[\\]","",smi)
  smi <- gsub("[[]C(@{1,}|)H[]]","C",smi)
  smi <- gsub("[[]CH2[]]","C",smi)
  smi
}



predict_mesh_classes  <- function(inputfile = "nameoftheinputfile") {

  ndf <- data.frame(readxl::read_xlsx(path = inputfile, sheet = 1), stringsAsFactors = F)

  if(!file.exists("cidmesh_smiles_fpbit.RData")) {
    load(url("https://github.com/barupal/ChemRICH/blob/master/cidmesh_smiles_fpbit.RData?raw=true"))
    save(df.mega.mesh, file = "cidmesh_smiles_fpbit.RData")
  }
  if(!file.exists("mesh_bit_loc_list.RData")) {
    load(url("https://github.com/barupal/ChemRICH/blob/master/mesh_bit_loc_list.RData?raw=true"))
    save(bitloclist, file = "mesh_bit_loc_list.RData")
  }
  if(!file.exists("treenames.df.RData")) {
    load(url("https://github.com/barupal/ChemRICH/blob/master/treenames.df.RData?raw=true"))
    save(treenames.df, file = "treenames.df.RData")
  }
  load("mesh_bit_loc_list.RData")
  load("cidmesh_smiles_fpbit.RData")
  load("treenames.df.RData")

  df.mega.mesh$CompoundName <- tolower(df.mega.mesh$CompoundName)
  treenames.df$ClassName <- tolower(treenames.df$ClassName)

  #########################
  ## Fatty acid labels ####
  #########################

  smi.all.fas <- as.character(sapply(ndf$smiles, makeSmiles.clean))
  falabelvec <- sapply(smi.all.fas, function(x) {
    elecount <- table(strsplit(gsub("[0-9]|[)]|[(]|=","",x),"")[[1]])
    falabel <-  ""
    if (length(table(c("c","o")%in%tolower(names(elecount))  )  )==1) {
      if(length(grep("n",x,ignore.case = T))==0) {
        if (elecount['C']>7 & length(grep("CCCC",x))==1 & length(grep("C2",x))!=1  ) { # long carbon but not aromatic or cyclic.
          if (elecount['O']==2) {
            dlen <- length(strsplit(x,"=")[[1]])-2
            falabel <- paste(c("FA",elecount['C'],dlen), collapse="_")
          }
          if(elecount['O']>=3) { ## Put Rules here. How many O and then how many carbon chain. That will make the class.
            if( length(grep("C1",x))==1) {
              if (length(strsplit(x,"C1")[[1]]) ==3 ) {
                dlen <- length(strsplit(x,"=")[[1]])-2
                #falabel <- paste(c("Prostaglandin",elecount['C']), collapse="_")
              } else {
                dlen <- length(strsplit(x,"=")[[1]])-2
                falabel <- paste(c("Epoxy FA",elecount['C']), collapse="_")
              }
            } else {
              if (length(strsplit(x,"=O|CO|OC")[[1]])-2==0){
                dlen <- length(strsplit(x,"=")[[1]])-2
                falabel <- paste(c("OH-FA",elecount['C'],dlen,(elecount['O']-2)), collapse="_")
              } else {
                if (length(strsplit(x,"OC|CO")[[1]]) <3 ) {
                  dlen <- length(strsplit(x,"=")[[1]])-2
                  falabel <- paste(c("O=FA",elecount['C'],dlen), collapse="_")
                }
              }
            }
          }
        }
      }
    }
    falabel
  })

  falabelvec[which(falabelvec=="OH-FA_20_3_2")] <- "DiHETrE"
  falabelvec[which(falabelvec=="OH-FA_20_4_2")] <- "DiHETE"
  falabelvec[which(falabelvec=="O=FA_18_3")] <- "oxo-ODE"
  falabelvec[which(falabelvec=="O=FA_20_5")] <- "oxo-ETE"
  falabelvec[which(falabelvec=="OH-FA_18_1_2")] <- "DiHOME"
  falabelvec[which(falabelvec=="OH-FA_18_1_3")] <- "TriHOME"
  falabelvec[which(falabelvec=="OH-FA_18_2_1")] <- "HODE"
  falabelvec[which(falabelvec=="OH-FA_18_2_2")] <- "DiHODE"
  falabelvec[which(falabelvec=="OH-FA_18_3_1")] <- "HOTrE"
  falabelvec[which(falabelvec=="OH-FA_20_3_1")] <- "HETrE"
  falabelvec[which(falabelvec=="OH-FA_20_4_1")] <- "HETE"
  falabelvec[which(falabelvec=="OH-FA_20_5_1")] <- "HEPE"
  falabelvec[which(falabelvec=="OH-FA_22_5_2")] <- "DiHDPE"
  falabelvec[which(falabelvec=="Epoxy FA_22")] <- "EpDPE"
  falabelvec[which(falabelvec=="Epoxy FA_18")] <- "EpETrE"
  falabelvec[which(falabelvec=="Epoxy FA_20")] <- "EpODE"
  falabelvec[grep("^FA_[0-9]{1,2}_0$", falabelvec)] <- "Saturated FA"
  falabelvec[grep("^FA_[0-9]{1,2}_[1-9]$", falabelvec)] <- "UnSaturated FA"

  cat("Computing sub-structure fingerprint\n")


  fps <- t(sapply(1:nrow(ndf), function(x) {
    #print(x)
    xy <- 0
    xy <- as.character(rcdk::get.fingerprint(rcdk::parse.smiles(ndf$smiles[x])[[1]],type="pubchem"))
    xy
  }))

  cid.mesh.df <- data.frame(CID=df.mega.mesh$CID[df.mega.mesh$CID%in%ndf$pubchem_id],MESHTREE= df.mega.mesh$MESSTREE[df.mega.mesh$CID%in%ndf$pubchem_id], stringsAsFactors = F)

  ## Covered by MeSH ###

  inmesh.vec <- rep("No",nrow(ndf))
  inmesh.vec[ndf$pubchem_id%in%cid.mesh.df$CID] <- "Yes"

  ## chemical similarity matrix.
  df1.bitmat <- do.call(rbind,lapply(fps,function(x) as.integer(strsplit(x,"")[[1]][1:881])))
  df1.bitmat.location <- lapply(1:nrow(df1.bitmat), function(x) { which(df1.bitmat[x,]==1) })
  only_b <- sapply(1:length(bitloclist), function(x) {  length(bitloclist[[x]]) })
  bitmeans <- sapply(1:length(bitloclist), function(x) {  median(bitloclist[[x]]) })
  fpsmeans <- sapply(df1.bitmat.location, function(x){median(x)})

  cat("Obtaining MESH class annotation\n")


  directlabels <- sapply(tolower(ndf$compound_name), function(x) {
    clabel="Not Found"
    findind <- which(df.mega.mesh$CompoundName==x)
    if(length( findind)>0) {
      classvec <- as.character(df.mega.mesh$MESSTREE[findind])
      classvec <- strsplit(classvec[1],";")[[1]]
      if( length(grep("^D01[.]",classvec)) > 0 ) {
        classvec <- classvec[-grep("^D01[.]|^D03[.]",classvec)]
      }
      clabel <- names(which.max(sapply(classvec,nchar)))
    }
    clabel
  })

  labelvec <- sapply(1:nrow(ndf), function(i) {
    clabel <- "Not Found"
    if(falabelvec[i]=="" & inmesh.vec[i]=="No" & directlabels[i]=="Not Found") {
      #if(falabelvec[i]=="") {
      #print(i)
      meanindex <- which(bitmeans < (fpsmeans[i]+5) & bitmeans > (fpsmeans[i]-5))
      bitloclist.sb <- bitloclist[meanindex]
      only_b.sb <- only_b[meanindex]
      overlapvec <- sapply(1:length(bitloclist.sb), function(x) {  length(which(bitloclist.sb[[x]]%in%df1.bitmat.location[[i]]==TRUE)) })
      tmvec <- overlapvec/((length(df1.bitmat.location[[i]])+only_b.sb)-overlapvec)
      if(length(which(tmvec>0.90))>0) {
        if(length(which(tmvec>0.98))>0){
          cidindex <- meanindex[which(tmvec>0.98)]
          if(length(cidindex)==1) {
            clabel <- df.mega.mesh$MESSTREE[cidindex]
          } else {
            clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec>0.98)])],"MeSHUID"])[1:10], function(x) {if(!is.na(x)) { strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID==x)][1],";")  }}))),decreasing = T)
            clabel.table <- which(clabel.table==clabel.table[1])
            clabel.table.sort <- sort(sapply(names(clabel.table),nchar),decreasing = T)
            clabel.table.sort.max <- which.max(clabel.table.sort)
            if(length(clabel.table.sort.max==1)) {
              clabel <- names(clabel.table.sort.max)
            } else {
              clabel <- sort(names(clabel.table.sort.max))[1]
            }
          }
        } else {
          ## CID_MESH
          cidindex <- meanindex[which(tmvec>0.90)]
          if(length(cidindex)==1) {
            clabel <- df.mega.mesh$MESSTREE[cidindex]
          } else {
            clabel.table <- sort(table(unlist(sapply(unique(df.mega.mesh[cidindex[order(tmvec[which(tmvec>0.90)])],"MeSHUID"])[1:10], function(x) {if(!is.na(x)) { strsplit(df.mega.mesh$MESSTREE[which(df.mega.mesh$MeSHUID==x)][1],";")  }}))),decreasing = T)
            clabel.table <- which(clabel.table==clabel.table[1])
            clabel.table.sort <- sort(sapply(names(clabel.table),nchar),decreasing = T)
            clabel.table.sort.max <- which.max(clabel.table.sort)
            if(length(clabel.table.sort.max==1)) {
              clabel <- names(clabel.table.sort.max)
            } else {
              clabel <- sort(names(clabel.table.sort.max))[1]
            }
          }
        }
      }
    }
    clabel
  })

  finalMesh.df <- cid.mesh.df
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%ndf$pubchem_id[which(falabelvec!="")]==FALSE),]
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%ndf$pubchem_id[which(directlabels!="Not Found")]==FALSE),]
  finalMesh.df <- rbind(finalMesh.df,data.frame(CID=ndf$pubchem_id, MESHTREE=labelvec ))
  finalMesh.df$NewMesh <- finalMesh.df$MESHTREE
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!="Not Found"),]
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!=""),]
  finalMesh.df <- finalMesh.df[!duplicated(finalMesh.df),]

  # Calculate the chemical similarity and clusters.

  cat("Computing Chemical Similarity\n")

  m <- df1.bitmat
  mat <- m%*%t(m)
  len <- length(m[,1])
  s <- mat.or.vec(len,len)
  for (i in 1:len) {
    for (j in 1:len){
      s[i,j] <- mat[i,j]/(mat[i,i]+mat[j,j]-mat[i,j])
    }
  }
  diag(s) <- 0
  hc <- hclust(as.dist(1-s), method="ward.D2") # ward method provide better clusters.
  clust1 <- cutreeDynamic(hc,distM = as.matrix(1-s),deepSplit =4, minClusterSize = 3)  # can be used to merge cluster, but it is better to keep them as it is. # Clusters are detected using the average linkage hclust.
  ndf$ClusterNumber <- clust1
  ndf$xlogp <- as.numeric(sapply(ndf$smiles, function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))

  ## creating the final label data frame.

  finalterm.df <- data.frame(CID=ndf$pubchem_id, Clabel=falabelvec,stringsAsFactors = F) # first add the fatty acid labels.
  directlabindex <- as.integer(which(directlabels!="Not Found"))[which(as.integer(which(directlabels!="Not Found"))%in%which(finalterm.df$Clabel=="")==TRUE)] ## then add the direct labels found by names matching
  finalterm.df$Clabel[directlabindex] <- as.character(directlabels[directlabindex])

  for (i in 1:nrow(finalterm.df)) {
    if(finalterm.df$Clabel[i]=="" & length(which(finalMesh.df$CID==ndf$pubchem_id[i]))>0 ) {
      finalterm.df$Clabel[i] <- names(which.max(sapply(unlist(strsplit(finalMesh.df$NewMesh[which(finalMesh.df$CID==ndf$pubchem_id[i])],";")),nchar)))
    }
  }

  ##########################################
  ####  Detect for new compound clusters ###
  ##########################################
  cat("Detecting New Compound Classes\n")

  finalterm.df.2 <- finalterm.df
  newClustVec <- names(which(table(ndf$ClusterNumber[which(finalterm.df.2$Clabel=="")])>3))
  clustMeanvec <- sapply(newClustVec, function(x) {  mean(s[which(ndf$ClusterNumber==x),which(ndf$ClusterNumber==x)])  }  )
  newClustVec <- newClustVec[which(clustMeanvec>0.70)]
  if(length(newClustVec)>0) {
    for(i in which(finalterm.df.2$Clabel=="")) {
      if(ndf$ClusterNumber[i]%in%newClustVec){
        finalterm.df$Clabel[i] <- paste0("NewCluster_",ndf$ClusterNumber[i])
      }
    }
  }

  ##### Map the compounds that have atleast 0.75 similarity to others. Only for compunds that do not have any labels.

  for ( i in which(finalterm.df$Clabel=="")) { ## if there is a metabolite that has score higher than 0.80 then we get the class using that compound.
    if(max(s[i,])>0.75) {
      simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.75)]
      simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if(simorder.class[1]!=""){
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      } else if(length(simorder.class)>1) {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      }
    }
  }

  cat("Detecting non-overlapping class definition\n")


  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))

  exclusionVec <- c("D02","D03.383","D03.633.100","D03.633.300","D03.633.400","D03.633","D03.605","D02.241.081") ## we will have some static list of terms that need to be excluded.
  exclusionVec <- c(exclusionVec, unique(falabelvec)[-1]) ## if we see Fatty acid label, we dont touch them.

  for ( i in which(finalterm.df$gCount<3)) { ## Drop the compound to the neareast one.
    qpat <- gsub("[.][0-9]{2,3}$","",finalterm.df$Clabel[i])
    if(length(grep(qpat,finalterm.df$Clabel))>2 & !qpat%in%exclusionVec){
      finalterm.df$Clabel[i] <- qpat
    }
  }

  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))

  for ( i in which(finalterm.df$gCount<3)) { ## Map to the closest ones.
    if(max(s[i,])>0.85) {
      simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.85)]
      simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if(simorder.class[1]!=""){
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      } else if(length(simorder.class)>1) {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      }
    }
  }

  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))

  # Repeat it one more time.
  for ( i in which(finalterm.df$gCount<3)) { ## Drop the compound to the neareast one.
    qpat <- gsub("[.][0-9]{2,3}$","",finalterm.df$Clabel[i])
    if(length(grep(qpat,finalterm.df$Clabel))>2 & !qpat%in%exclusionVec){
      finalterm.df$Clabel[i] <- qpat
    }
  }
  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))

  for ( i in which(finalterm.df$gCount<3)) { ## Map to the closest ones.
    if(max(s[i,])>0.85) {
      simorder <- order(s[i,],decreasing = T)[which(s[i,][order(s[i,],decreasing = T)]>0.85)]
      simorder.class <- sapply(simorder, function(x) { finalterm.df$Clabel[x]})
      simorder.class <- simorder.class[!is.na(simorder.class)]
      if(simorder.class[1]!=""){
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      } else if(length(simorder.class)>1) {
        finalterm.df$Clabel[i] <- simorder.class[which(simorder.class!="")][1]
      }
    }
  }

  finalterm.df$Count <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(which(finalterm.df$Clabel==x))  }))
  finalterm.df$gCount <- as.integer(sapply(finalterm.df$Clabel, function(x) { length(grep(x,finalterm.df$Clabel))  }))
  finalterm.df$Clabel[which(finalterm.df$Count<3)] <- finalterm.df.2$Clabel[which(finalterm.df$Count<3)] ### We reverse back the original labels as this did not create any higher labels.

  finallabelvec <- finalterm.df$Clabel

  HasSaturatedFats <- names(which(table(finallabelvec[grep("D10|D09.400.410",finallabelvec)[which(sapply(grep("D10|D09.400.410",finallabelvec), function(x)  { length(grep("C=C",ndf$smiles[x]))  })==0)]])>2)) ### we are only selecting lipid classes that has atleast 3 saturated lipids.

  for (i in 1:nrow(finalterm.df)) {
    if(finallabelvec[i]%in%HasSaturatedFats){
      if(length(grep("C=C",ndf$smiles[i]))==0) {
        finallabelvec[i] <- paste0("Saturated_",getCNames(finallabelvec[i]))
      } else {
        finallabelvec[i] <- paste0("Unsaturated_",getCNames(finallabelvec[i]))
      }
    }
  }
  clusterids <- sapply(as.character(finallabelvec),getCNames)
  cat("ChemRICH Class Estimation Finished. ChemRICHClusters column has been added to ndf data.frame\n")
  ndf$MeSH_Class <- clusterids

  l <- list("Results" = ndf )
  openxlsx::write.xlsx(l, file = paste0("MeSh_Prediction_Results.xlsx"), asTable = TRUE)
  cat(paste0("MeSh_Prediction_Results.xlsx", " has been saved.\n Voilla!!"))
}
