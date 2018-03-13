#raw_data <- readChar("./inst/www/chemrich_example_template.txt", nchars = 10000000)
#stat_file = raw_data

load("./data/cidmesh_smiles_fpbit.RData")  ###df.mega.mesh.unique Unique CIDS
load("./data/mesh_bit_loc_list.RData") #bitloclist bit locations for the unique CIDS.
load("./data/treenames.df.RData") ### Name of the TreeBranches. We shall put them on the
load("./data/PubChemMeshCounts.RData") # pubchemdata counts of CIDs and PMIDS for each mesh category.
load("./data/cidsbiosys.RData")

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

arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02,
                          lab.offset=1.06,cex=1,orientation="curved",...){
  ## Credit: this function is adopted from http://blog.phytools.org/2017/03/clade-labels-on-circular-fan-tree.html
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(obj$type!="fan") stop("method works only for type=\"fan\"")
  h<-max(sqrt(obj$xx^2+obj$yy^2))
  if(is.null(tree)){
    tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
               Nnode=obj$Nnode)
    class(tree)<-"phylo"
  }
  d<-getDescendants(tree,node)
  d<-sort(d[d<=Ntip(tree)])
  deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
  ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
  deg[ii]<-360+deg[ii]
  draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),deg2=max(deg),lwd = 2)
  if(orientation=="curved")
    arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex)
  else if(orientation=="horizontal"){
    x0<-lab.offset*cos(median(deg)*pi/180)*h
    y0<-lab.offset*sin(median(deg)*pi/180)*h
    text(x=x0,y=y0,label=text,
         adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
         offset=0)
  }
}

getChemRich_windows <- function (stat_file,cutoff=0.1) {
  letters.x <- c(letters,LETTERS)
  pacman::p_load(grid,rcdk, RJSONIO,RCurl, dynamicTreeCut,ape,ggplot2,ggrepel,ReporteRs,XLConnect,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,XLConnect)
  loadfonts(quiet = T)
  stat_file <- gsub("\r","",stat_file)
  cfile <- strsplit(stat_file,"\n")[[1]]
  df1 <- do.call(rbind, lapply(cfile, function (x) { strsplit(x,"\t")[[1]]  } ))
  #df1 <- df1[c(1:96),]
  colnames(df1) <- sapply(df1[1,],as.character)
  df1 <- df1[-1,]
  df1 <- data.frame(df1,stringsAsFactors = F)
  df1$foldchange <- sapply(df1$foldchange, function(x) { as.numeric(as.character(x))  })
  df1$pvalue <- sapply(df1$pvalue, function(x) { as.numeric(as.character(x))  })
  df1$CID <- as.integer(df1$Pubchem.ID)
  #df1$NewCID <- as.character(df1$CID)

  ### If SMILES code is missing for even one metabolite, we will break the script and ends here.
  if( (length(which(is.na(df1$SMILES)==TRUE)))) { stop("Missing SMILES codes. Please check the input.") }
  if( (length(which(is.na(df1$CID)==TRUE)) & length(which(is.na(df1$SMILES)==TRUE)))) { stop("Missing SMILES codes. Please check the input.") }

  df.mega.mesh$CompoundName <- tolower(df.mega.mesh$CompoundName)
  treenames.df$ClassName <- tolower(treenames.df$ClassName)

  ###########################################
  #### Detection of Fatty Acid Clusters #####
  ###########################################

  smi.all.fas <- as.character(sapply(df1$SMILES, makeSmiles.clean))
  falabelvec <- sapply(smi.all.fas, function(x) {
    elecount <- table(strsplit(gsub("[0-9]|[)]|[(]|=","",x),"")[[1]])
    falabel <-  ""
    if (length(which(names(elecount)%in%c("C","O")==FALSE))==0) {
      if (elecount['C']>7 & length(grep("CCCC",x))==1 & length(grep("C2",x))!=1  ) { # long carbon but not aromatic or cyclic.
        if (elecount['O']==2) {
          dlen <- length(strsplit(x,"=")[[1]])-2
          falabel <- paste(c("FA",elecount['C'],dlen), collapse="_")
        }
        if(elecount['O']>=3) { ## Put Rules here. How many O and then how many carbon chain. That will make the class.
          if( length(grep("C1",x))==1) {
            if (length(strsplit(x,"C1")[[1]]) ==3 ) {
              dlen <- length(strsplit(x,"=")[[1]])-2
              falabel <- paste(c("Prostaglandin",elecount['C']), collapse="_")
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

  exc <- XLConnect::loadWorkbook("ChemRICH_results.xlsx", create = T)

  ### Get the FP annotations from the CouchDB using CIDS
  idlist <- list()
  idlist$keys <- as.integer(df1$CID)[which(is.na(df1$CID)==FALSE)]
  urlres <- getURL("http://chemrich.fiehnlab.ucdavis.edu/db/chemrichdb_cidfp/_design/data/_view/cid_fp",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=RJSONIO::toJSON(idlist))
  urlres.list <- RJSONIO::fromJSON(urlres)
  cid.fp.df <- as.data.frame(do.call(rbind, urlres.list$rows), stringsAsFactors = F)

  ### Get the FP annotations from the CouchDB using SMILES
  idlist <- list()
  idlist$keys <- df1$SMILES
  urlres <- getURL("http://chemrich.fiehnlab.ucdavis.edu/db/chemrichdb_smilesfp/_design/data/_view/smiles_fp",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=RJSONIO::toJSON(idlist))
  urlres.list <- RJSONIO::fromJSON(urlres)
  smiles.fp.df <- as.data.frame(do.call(rbind, urlres.list$rows), stringsAsFactors = F)

  fps <- t(sapply(1:nrow(df1), function(x) {
    #print(x)
    xy <- 0
    if(is.na(df1$CID[x]) | length(which(smiles.fp.df$key==df1$SMILES[x]))!=0  ) {
      whichsmiles <- which(smiles.fp.df$key==df1$SMILES[x])
      if (length(whichsmiles)!=0) {
        xy <- smiles.fp.df$value[whichsmiles][[1]]
      } else {
        xy <- as.character(rcdk::get.fingerprint(rcdk::parse.smiles(df1$SMILES[x])[[1]],type="pubchem"))
      }
    } else {
      whichcid <- which(cid.fp.df$key==df1$CID[x])
      if (length(whichcid)!=0) {
        xy <- cid.fp.df$value[whichcid][[1]]
      } else {
        xy <- as.character(rcdk::get.fingerprint(rcdk::parse.smiles(df1$SMILES[x])[[1]],type="pubchem"))
      }
    }
    xy
  }))

  ### If any of the smiles codes are wrong. Break the code here.
  if( (length(which(fps==0))>0)) { stop("Incorrect SMILES Code provided. Please check the input") }

  misseddf <- data.frame(SMILES=(df1$SMILES[which(df1$SMILES%in%smiles.fp.df$key!=TRUE)]), FP=fps[which(df1$SMILES%in%smiles.fp.df$key!=TRUE)], stringsAsFactors = F)
  if(nrow(misseddf)!=0) {
    elist <- list()
    for (i in 1:nrow(misseddf)) {
      elist[[i]] <- list()
      elist[[i]]$SMILES<- misseddf$SMILES[i]
      elist[[i]]$FP <- misseddf$FP[i]
    }
    elist1<- list()
    elist1$docs <- elist
    urlres <- getURL("http://chemrich.fiehnlab.ucdavis.edu/db/chemrichdb_smilesfp/_bulk_docs",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=RJSONIO::toJSON(elist1))
    gc()
  }

  ################################################################################
  ###### Query Against All the MESH compounds and get back the annotations #######
  ################################################################################

  ## CID_MESH
  idlist <- list()
  idlist$keys <- as.integer(df1$CID)[which(is.na(df1$CID)==FALSE)]
  urlres <- getURL("http://chemrich.fiehnlab.ucdavis.edu/db/chemrichdb_cid/_design/data/_view/cid_tree",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=RJSONIO::toJSON(idlist))
  urlres.list <- RJSONIO::fromJSON(urlres)
  cid.mesh.df <- as.data.frame(do.call(rbind, urlres.list$rows), stringsAsFactors = F)
  cid.mesh.df$CID <- unlist(cid.mesh.df$key)

  ## SMILES_MESH
  idlist <- list()
  idlist$keys <- df1$SMILES
  urlres <- getURL("http://chemrich.fiehnlab.ucdavis.edu/db/chemrichdb_smiles/_design/data/_view/smiles_tree",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=RJSONIO::toJSON(idlist))
  urlres.list <- RJSONIO::fromJSON(urlres)
  smiles.mesh.df <- as.data.frame(do.call(rbind, urlres.list$rows), stringsAsFactors = F)
  smiles.mesh.df$SMILES <- unlist(smiles.mesh.df$key)

  df1$CID[which(is.na(df1$CID)==TRUE)] <- paste0("cid_",which(is.na(df1$CID)==TRUE))

  smiles_cid_map <- as.data.frame(do.call(rbind,lapply(1:nrow(smiles.mesh.df), function(x) {  cbind(smiles.mesh.df[x,],CID=df1$CID[which(df1$SMILES==smiles.mesh.df$SMILES[x])])  })), stringsAsFactors = F)

  ##########################################################################################################
  ######################## Get the MeSH annotation by similarity scoring against mesh compounds. ###########
  ##########################################################################################################

  inmesh.vec <- rep("No",nrow(df1))
  inmesh.vec[unique(c(which(df1$CID%in%cid.mesh.df$CID==TRUE), which(df1$SMILES%in%smiles.mesh.df$SMILES==TRUE)))] <- "Yes"
  df1.bitmat <- do.call(rbind,lapply(fps,function(x) as.integer(strsplit(x,"")[[1]][1:881])))
  df1.bitmat.location <- lapply(1:nrow(df1.bitmat), function(x) { which(df1.bitmat[x,]==1) })
  only_b <- sapply(1:length(bitloclist), function(x) {  length(bitloclist[[x]]) })
  bitmeans <- sapply(1:length(bitloclist), function(x) {  median(bitloclist[[x]]) }) # we decided to use median, as some outliers values shall not affect it.
  fpsmeans <- sapply(df1.bitmat.location, function(x){median(x)})

  ################################################
  ## Direct name searching against MeSH database.#
  ################################################

  directlabels <- sapply(tolower(df1$Compound.Name), function(x) { ## We will ask everyone to use MeSH compound names to enable automated text mining on the metabolomics datasets.
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

  labelvec <- sapply(1:nrow(df1), function(i) {
    clabel <- "Not Found"
    if(falabelvec[i]=="" & inmesh.vec[i]=="No" & directlabels[i]=="Not Found") {
      #if(falabelvec[i]=="") {
      print(i)
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

  xmllist <- which(labelvec!="Not Found")
  if(length(xmllist)!=0) { ## we want to post these high score ones to the couchdb.
    elist <- list()
    for (i in 1:length(xmllist)) {
      elist[[i]] <- list()
      elist[[i]]$SMILES <- df1$SMILES[xmllist[i]]
      elist[[i]]$MeSHTree <- labelvec[xmllist[i]]
    }
    elist1<- list()
    elist1$docs <- elist
    ##urlres <- getURL("http://localhost:5984/chemrichdb_smiles/_bulk_docs",customrequest='POST',httpheader=c('Content-Type'='application/json'),postfields=toJSON(elist1))
    gc()
  }

  smiles_cid_map <- smiles_cid_map[,-4]

  ##########################################
  ### Preparing the final CID_MESH mapping.####
  #########################################

  finalMesh.df <- as.data.frame(rbind(cid.mesh.df, smiles_cid_map))
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%df1$CID[which(falabelvec!="")]==FALSE),] ## remove the compounds that are covered by FA labels.
  finalMesh.df <- finalMesh.df[which(finalMesh.df$CID%in%df1$CID[which(directlabels!="Not Found")]==FALSE),] ## remove the compounds that have been covered by the direct label mapping.

  finalMesh.df <- finalMesh.df[,-1]
  finalMesh.df$value <- unlist(finalMesh.df$value)
  finalMesh.df <- rbind(finalMesh.df,data.frame(key=df1$CID, value=labelvec,CID=df1$CID))
  meshvec <- finalMesh.df$value
  finalMesh.df$NewMesh <- meshvec
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!="Not Found"),]
  finalMesh.df <- finalMesh.df[which(finalMesh.df$NewMesh!=""),]
  finalMesh.df <- finalMesh.df[!duplicated(finalMesh.df),]


  #############################################################
  #### Calculation of the simialrity tree and its clustering ##
  #############################################################

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
  ##

  hc <- hclust(as.dist(1-s), method="average") # ward method provide better clusters.
  clust1 <- cutreeDynamic(hc,distM = as.matrix(1-s),deepSplit =4, minClusterSize = 2)  # can be used to merge cluster, but it is better to keep them as it is. # Clusters are detected using the average linkage hclust.
  glaydf2 <- data.frame(df1,cluster=clust1,stringsAsFactors = F)
  df1.order <- df1[hc$order,]

  #### Clustering order for the tree ordering and the volcano plot ordering.
  s2 <- s[order(df1$Compound.Name),order(df1$Compound.Name)]
  hc2 <- hclust(as.dist(1-s2), method="ward.D") # We will use average to detect cluster on the tree and ward.d to sort the clusters and individual compounds. Apparently ward.d is better in sorting the compounds by their degree of unsaturation. It is also interesting to know that the original order in which the data are had influence on the sorting order.
  df2.order <- df1[order(df1$Compound.Name),][hc2$order,]

  ###########################################
  #### finding the Non-overlapping terms. ###
  ###########################################

  ##
  #### Creation of label dataframe.
  ##

  finalterm.df <- data.frame(CID=df1$CID, Clabel=falabelvec,stringsAsFactors = F) # first add the fatty acid labels.
  directlabindex <- as.integer(which(directlabels!="Not Found"))[which(as.integer(which(directlabels!="Not Found"))%in%which(finalterm.df$Clabel=="")==TRUE)] ## then add the direct labels found by names matching
  finalterm.df$Clabel[directlabindex] <- as.character(directlabels[directlabindex])

  for (i in 1:nrow(finalterm.df)) { # now add the mesh ids obtained from CouchDB. This will include the mesh annotation calcualted previously.
    if(finalterm.df$Clabel[i]=="" & length(which(finalMesh.df$CID==df1$CID[i]))>0 ) {
      print(i)
      finalterm.df$Clabel[i] <- names(which.max(sapply(unlist(strsplit(finalMesh.df$NewMesh[which(finalMesh.df$CID==df1$CID[i])],";")),nchar)))
    }
  }

  ##########################################
  ####  Detect for new compound clusters ###
  ##########################################
  finalterm.df.2 <- finalterm.df
  #finalterm.df.2$Clabel[which(finalterm.df$Clabel=="D09.400.410.420.525.870")] <- "" #we first remove the sphingolyelins to test it.
  #finalterm.df.2$Clabel[which(finalterm.df$Clabel=="D10.351.801")] <- "" #we first remove the sphingolyelins to test it.
  newClustVec <- names(which(table(glaydf2$cluster[which(finalterm.df.2$Clabel=="")])>3))
  clustMeanvec <- sapply(newClustVec, function(x) {  mean(s[which(glaydf2$cluster==x),which(glaydf2$cluster==x)])  }  )
  newClustVec <- newClustVec[which(clustMeanvec>0.70)]
  if(length(newClustVec)>0) {
    for(i in which(finalterm.df.2$Clabel=="")) {
      if(glaydf2$cluster[i]%in%newClustVec){
        finalterm.df$Clabel[i] <- paste0("NewCluster_",glaydf2$cluster[i])
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

  ### Find the saturated and unsaturated compounds in the
  HasSaturatedFats <- names(which(table(finallabelvec[grep("D10|D09.400.410",finallabelvec)[which(sapply(grep("D10|D09.400.410",finallabelvec), function(x)  { length(grep("C=C",df1$SMILES[x]))  })==0)]])>2)) ### we are only selecting lipid classes that has atleast 3 saturated lipids.

  for (i in 1:nrow(finalterm.df)) {
    if(finallabelvec[i]%in%HasSaturatedFats){
      if(length(grep("C=C",df1$SMILES[i]))==0) {
        finallabelvec[i] <- paste0("Saturated_",getCNames(finallabelvec[i]))
      } else {
        finallabelvec[i] <- paste0("Unsaturated_",getCNames(finallabelvec[i]))
      }
    }
  }

  #############################################
  ### Calculation of Enrichment Statistics#####
  #############################################

  clusterids <- names(which(table(sapply(as.character(finallabelvec),getCNames))>2))
  clusterids <- clusterids[which(clusterids!="")]
  df1$ClusterLabel <- as.character(sapply(as.character(finallabelvec),getCNames))

  cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
    cl.member <- which(df1$ClusterLabel==x)
    if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
      pval.cl.member <- df1$pvalue[cl.member]
      p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
      p.test.results$p.value
    } else {
      1
    }
  })

  cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\

  #clusterdf <- data.frame(name=clusterids[which(cluster.pvalues!=10)],pvalues=cluster.pvalues[which(cluster.pvalues!=10)], stringsAsFactors = F)
  clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)

  # clusterdf$keycpd <- sapply(clusterdf$name, function(x) {
  #   dfx <- df1[which(df1$ClusterLabel==x),]
  #   dfx$SMILES[which.min(dfx$pvalue)]
  # })
  clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
    dfx <- df1[which(df1$ClusterLabel==x),]
    dfx$Compound.Name[which.min(dfx$pvalue)]
  })

  altrat <- sapply(clusterdf$name, function (k) {
    length(which(df1$ClusterLabel==k & df1$pvalue<0.10))/length(which(df1$ClusterLabel==k))
  })

  uprat <-sapply(clusterdf$name, function (k) {
    length(which(df1$ClusterLabel==k & df1$pvalue<0.10 & df1$foldchange > 1.00000000))/length(which(df1$ClusterLabel==k))
  })

  clust_s_vec <- sapply(clusterdf$name, function (k) {
    length(which(df1$ClusterLabel==k))
  })

  clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$ClusterLabel==k & df1$pvalue<0.10))})
  clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(df1$ClusterLabel==k & df1$pvalue<0.10 & df1$foldchange > 1.00000000))})
  clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(df1$ClusterLabel==k & df1$pvalue<0.10 & df1$foldchange < 1.00000000))})
  clusterdf$upratio <- uprat
  clusterdf$altratio <- altrat
  clusterdf$csize <- clust_s_vec
  clusterdf <- clusterdf[which(clusterdf$csize>2),]
  clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")
  clustdf <- clusterdf

  clustdf.e <- clustdf[order(clustdf$pvalues),]
  clustdf.e$pvalues <- signif(clustdf.e$pvalues, digits = 2)
  clustdf.e$adjustedpvalue <- signif(clustdf.e$adjustedpvalue, digits = 2)
  clustdf.e$upratio <- signif(clustdf.e$upratio, digits = 1)
  clustdf.e$altratio <- signif(clustdf.e$altratio, digits = 1)

  clustdf.e <- clustdf.e[,c("name","csize","pvalues","adjustedpvalue","keycpdname","alteredMetabolites","upcount","downcount","upratio","altratio")]

  names(clustdf.e) <- c("Cluster name","Cluster size","p-values","FDR","Key compound","Altered metabolites","Increased","Decreased","Increased ratio","Altered Ratio")

  XLConnect::createSheet(exc,'ChemRICH_Results')
  XLConnect::writeWorksheet(exc, clustdf.e, sheet = "ChemRICH_Results", startRow = 1, startCol = 2)


  #write.table(clustdf.e, file="cluster_level_results_altered.txt", col.names = T, row.names = F, quote = F, sep="\t")
  #writeLines(toJSON(clustdf), "clusterJson.json")
  gdf <- datatable(clustdf.e,options = list(pageLength = 10),rownames = F)
  gdf$width  <- "auto"
  gdf$height <- "auto"
  saveWidget(gdf,file="clusterlevel.html",selfcontained = F)

  clusterdf$Compounds <- sapply(clusterdf$name, function(x) {
    dfx <- df1[which(df1$ClusterLabel==x),]
    paste(dfx$Compound.Name,collapse="<br>")
  }) ## this one is the label on the tooltip of the ggplotly plot.

  clustdf <- clusterdf[which(cluster.pvalues!=1),]
  #################################################
  ########## Impact Visualization Graph ###########
  #################################################

  clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]

  clustdf.alt.impact <- clustdf.alt.impact[order(sapply(clustdf.alt.impact$keycpdname, function(x) {which( df2.order$Compound.Name==x)})),]
  clustdf.alt.impact$order <- 1:nrow(clustdf.alt.impact) ### Order is decided by the hclust algorithm.
  clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)

  p2 <- ggplot(clustdf.alt.impact,aes(x=order,y=-log(pvalues)))
  p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
    #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
    scale_color_gradient(low = "blue", high = "red")+
    scale_size(range = c(5, 30)) +
    scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
    scale_x_continuous(" cluster order on the similarity tree ") +
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
  wbp <- pptx(template = system.file("data","chem_rich_temp.pptx", package = "ChemRICH"))
  #wbp <- pptx(template = "./data/chem_rich_temp.pptx")
  wbp <- addSlide( wbp, "lipidClust" )
  wbp <- addPlot( wbp,  function() print(p2), offx =0.0 , offy = 0.0, width = 8, height = 6 , vector.graphic = TRUE )
  writeDoc( wbp, file = "chemrich_impact_plot.pptx" )
  ggsave("chemrich_impact_plot.png", p2,height = 8, width = 12, dpi=300)

  #ggsave("tst.png",height=9,width=12,dpi=72)

  ##########################################################
  ##### Interactive Visualization plot using GGPlotLY ######
  ###########################################################

  p2 <- ggplot(clustdf.alt.impact,aes(label=name,label2=pvalues, label3=csize,label4=Compounds))
  p2 <- p2 + geom_point(aes(x=order,y=-log(pvalues),size=csize, color=upratio)) + scale_color_gradient(low = "blue", high = "red")+ scale_size(range = c(5, 30)) +
    #labs(caption = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
    scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+5  )) +
    scale_x_continuous(" cluster order on the similarity tree ") +
    theme_bw() +
    #labs(title = "ChemRICH cluster impact plot") +
    geom_text(aes(x=order,y=-log(pvalues),label = name), color = "gray20",data=subset(clustdf.alt.impact, csize>2))+
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
  saveWidget(gg,file = "ggplotly.html", selfcontained = F)

  #################################
  ### Interactive volcano plot  ###
  #################################
  # we need to add the interactive volcano plot for the p-values and fold-change sorted by the MeSH tree order.
  df2 <- df2.order
  df2$Changed <- "No Change"
  df2$Changed[which(df2$pvalue<0.05 & df2$foldchange>1)] <- "UP"
  df2$Changed[which(df2$pvalue<0.05 & df2$foldchange<1)] <- "DOWN"
  df2$Changed <- as.factor(df2$Changed)
  df2$pathway <- "No"
  df2$pathway[which(df2$CID%in%cid_biosys==TRUE)] <- "yes"
  df2$pathway <- as.factor(df2$pathway)
  df2$Compound.Name <- factor(df2$Compound.Name, levels =df2$Compound.Name)
  df2$foldchange <- round(sapply(df2$foldchange, function(x) { if(x>1) {x} else {1/x} }), digits = 1)
  df2$foldchange[ df2$foldchange>5] <- 5

  p2 <-   ggplot(df2, aes(label=Compound.Name,x=Compound.Name, y=-log(pvalue,base = 10),colour = Changed,shape=pathway, size=foldchange)) +  scale_size(range = c(1, 10)) +
    #geom_line(position=pd, size=2)+
    #geom_errorbar(aes(ymin = V2-V3 , ymax=V2+V3), width=.3,size=2,position=pd) +
    geom_point(stat = "identity") + # 21 is filled circle
    #geom_bar(stat="identity", size=.1,position=position_dodge()) +
    scale_y_continuous("pvalue (-log)") +
    scale_x_discrete("Metabolites: Red- increased,blue-decreased,yellow-not significant, solid-pathway(s) found ") +
    scale_color_manual("Student ttest",values=c("blue", "yellow", "red","white")) +
    scale_fill_manual("",values=c("white", "yellow", "red","white")) +
    scale_shape_manual("Pathway found",values=c(1,16))+
    #scale_shape(solid = FALSE) +
    theme_bw() +
    labs(title = "Metabolic Dys-regulations sorted by chemical similarity") +
    theme(
      plot.title = element_text(face="bold", size=30,hjust = 0.5),
      axis.title.x = element_text(face="bold", size=20),
      axis.title.y = element_text(face="bold", size=30, angle=90),
      panel.grid.major = element_blank(), # switch off major gridlines
      panel.grid.minor = element_blank(), # switch off minor gridlines
      #legend.justification=c(1,0),
      #legend.position=c(1,.6),
      legend.position = "none",
      #legend.title = element_blank(), # switch off the legend title
      #legend.text = element_text(size=12),
      #legend.key.size = unit(1.5, "lines"),
      #legend.key = element_blank(), # switch off the rectangle around symbols in the legend
      #legend.spacing = unit(.05, "cm"),
      #axis.text.x = element_text(size=15,angle = 45, hjust = 1.0),
      axis.text.x= element_blank(),
      axis.text.y = element_text(size=15,angle = 0, hjust = 0.5)
    )
  p2
  p3 <- ggplotly(p2, width = 1600, height = 1000)
  htmlwidgets::saveWidget(p3, "dysregplot.html", selfcontained = F)

  ######################################################
  ### Visualization of Chemical Tree ###
  ######################################################

  treeLabels <- rep("",nrow(df1))

  plot.fan.chemrich <- function(hc, clus, dirvec,sizevec) {
    nclus <- length(unique(clus))
    palette <- c('black','green','orange','blue','grey','yellow','pink','brown','purple','violet','skyblue','khaki','lavender','magenta',"gold","sienna","tan","seagreen","orchid","linen","skyblue3","wheat","navyblue","moccasin","navy","dodgerblue","deeppink","chocolate",'red','blue','green','orange','maroon2','grey','yellow','pink','brown','purple','violet','skyblue','khaki','lavender','magenta',"gold","sienna","tan","seagreen","orchid","linen","skyblue3","wheat","navyblue","moccasin","navy","dodgerblue","deeppink","chocolate",'black','green','orange','blue','grey','yellow','pink','brown','purple','violet','skyblue','khaki','lavender','magenta',"gold","sienna","tan","seagreen","orchid","linen","skyblue3","wheat","navyblue","moccasin","navy","dodgerblue","deeppink","chocolate",'red','blue','green','orange','maroon2','grey','yellow','pink','brown','purple','violet','skyblue','khaki','lavender','magenta',"gold","sienna","tan","seagreen","orchid","linen","skyblue3","wheat","navyblue","moccasin","navy","dodgerblue","deeppink","chocolate")
    #[1:nclus]
    X <- as.phylo(hc)
    X$tip.label <- df1$Compound.Name
    #X$tip.label <- as.character(clus)
    X$tip.label[which(dirvec=="yellow")] <- ""
    edge.colors <- rep("lightgray",nrow(X$edge))
    for (i in unique(clus)) {
      if(i>0) {
        difvec <- diff(which(X$edge[,2] %in% which(clus==i)))
        if(length(which(difvec>3))==0) {
          edge.colors[min(which(X$edge[,2] %in% which(clus==i))):max(which(X$edge[,2] %in% which(clus==i)))] <- "black"
          edge.colors[min(which(X$edge[,2] %in% which(clus==i)))-1] <- "black"
          #edge.colors[max(which(X$edge[,2] %in% which(clus==i)))+1] <- "black"
          #nodelabels(LETTERS[k], node=   74, adj=c(0.5,0.5), frame = "c", bg = "white", cex = 1.0)
        } else {
          ovec <- which(X$edge[,2] %in% which(clus==i))
          ovec <- ovec[-1]
          edge.colors[min(ovec):max(ovec)] <- "black"
          edge.colors[min(ovec)-1] <- "black"
        }
      }
    }
    XX <- plot(X,type='fan', tip.color="black",edge.color=edge.colors,show.tip.label = F,edge.width = 2, label.offset=.01, cex=0.5,no.margin=T)
    tiplabels(pch = 21,col = dirvec, bg = dirvec,  cex = sizevec)
    #edgelabels()
    k = 1
    for (i in unique(clus)) {
      if(i>0) {
        difvec <- diff(which(X$edge[,2] %in% which(clus==i)))
        if(length(which(difvec>3))==0) {
          nodeiddf <- as.data.frame(X$edge[min(which(X$edge[,2] %in% which(clus==i))):max(which(X$edge[,2] %in% which(clus==i))),])
          #nodelabels(letters[k], node=min(nodeiddf$V1), adj=c(0.5,0.5), frame = "none", bg = "transparent", cex = 2.0,col="green")
          arc.cladelabels(text=letters.x[k],node=min(nodeiddf$V1),cex = 1.0,lab.offset=1.05,ln.offset=1.01)
          k = k+1
        } else {
          ovec <- which(X$edge[,2] %in% which(clus==i))
          ovec <- ovec[-1]
          nodeiddf <- as.data.frame(X$edge[min(ovec):max(ovec),])
          #nodelabels(letters[k], node=min(nodeiddf$V1), adj=c(0.5,0.5), frame = "none", bg = "transparent", cex = 2.0,col="green")
          arc.cladelabels(text=letters.x[k],node=min(nodeiddf$V1),cex = 1.0,lab.offset=1.05,ln.offset=1.01)
          k = k+1
        }
      }
    }
  }
  clus <- as.integer(clust1)
  if(length(which(clus==0))>0) {
    clus[which(clus==0)] <- max(unique(clust1))+1
  }
  sizevec <- rep(1.0,length(clus))
  dirvec <- sapply(glaydf2$foldchange,function(x){ if(x>1) { return("red") } else (return("blue")) })
  dirvec <- sapply(1:length(glaydf2$pvalue), function(x) { if(glaydf2$pvalue[x]>0.05) {"yellow"} else{dirvec[x]}})
  sizevec[which(dirvec=="yellow")] <- 0.2
  clus_pval_list <-  sapply(unique(df1$ClusterLabel), function(x) { sapply(df1$ClusterLabel[which(df1$ClusterLabel==x)], function(y) { clusterdf$pvalues[which(clusterdf$name==y)]  }) })
  alteredClusters <- unique(glaydf2$cluster)[which(unlist(lapply(clus_pval_list, function(xx) { length(which(xx<0.05))/length(xx) }))>0.5)]
  clus[which(!clus%in%alteredClusters)] <- 0 ## All the clusters that are not altered are turned off with this command

  X <- as.phylo(hc)
  k = 1
  for (i in unique(clus)) {
    if(i>0) {
      difvec <- diff(which(X$edge[,2] %in% which(clus==i)))
      if(length(which(difvec>3))==0) {
        nodeiddf <- as.data.frame(X$edge[min(which(X$edge[,2] %in% which(clus==i))):max(which(X$edge[,2] %in% which(clus==i))),])

        treeLabels[which(glaydf2$cluster==i)] <- letters.x[k]
        k = k+1
      } else {
        ovec <- which(X$edge[,2] %in% which(clus==i))
        ovec <- ovec[-1]
        nodeiddf <- as.data.frame(X$edge[min(ovec):max(ovec),])

        treeLabels[which(glaydf2$cluster==i)] <- letters.x[k]
        k = k+1
      }
    }
  }

  wbp <- pptx(template = system.file("data","chem_rich_temp.pptx", package = "ChemRICH" )) # use this one when using the installed packages.
  #wbp <- pptx(template = "./data/chem_rich_temp.pptx")
  wbp <- addSlide( wbp, "lipidClust" )
  wbp <- addPlot( wbp, function() plot.fan.chemrich(hc,clus, dirvec,sizevec), offx =0.1 , offy = -0.1, width = 8, height = 8 , vector.graphic = FALSE )

  #wbp <- addParagraph( wbp, "Compund cluster annotation mapping is provided in the xlsx output file." , offx =0.1 , offy = -0.1, width = 8, height = 8  )

  writeDoc( wbp, file = paste0("chemrich_output_tree.pptx") )
  png(file="chemrich_tree.png",width=12,height=15,units="in",res=300)
  plot.fan.chemrich(hc,clus, dirvec,sizevec)
  text(+0.0, +0.8, "ChemRICH : Chemical Similarity Enrichment Analysis", cex = 2.0)
  dev.off()

  ## Export the final compounds result table.

  df1$TreeLabels <- treeLabels

  df1$pvalue <- signif(df1$pvalue, digits = 2)
  df1$foldchange <- signif(df1$foldchange, digits = 2)
  df1$FDR <- signif(  p.adjust(df1$pvalue), digits = 2)

  gdf <- datatable(df1,options = list(pageLength = 10), rownames = F)

  gdf$width  <- "auto"
  gdf$height <- "auto"
  saveWidget(gdf,file="compoundlevel.html",selfcontained = F)

  XLConnect::createSheet(exc,'Compound_ChemRICH')
  XLConnect::writeWorksheet(exc, df1, sheet = "Compound_ChemRICH", startRow = 1, startCol = 2)

  ### Final Export to XLSx
  XLConnect::saveWorkbook(exc)

}

