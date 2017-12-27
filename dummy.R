# dummy R code

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
