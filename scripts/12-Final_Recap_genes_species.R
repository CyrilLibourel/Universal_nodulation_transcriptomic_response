# Script to count HOGs on species tree
######################
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ComplexHeatmap","circlize"), require, character.only = TRUE)

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")

# load orthogroup file
HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

tables <- paste0(list.files(path="./rnaseq_timepoints/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

######## Multi species ########
Resume_HOGs_Up <- fread(paste0("./results/Nod_final_table_Up.txt"))
Resume_HOGs_Down <- fread(paste0("./results/Nod_final_table_Down.txt"))

for(sp in species){
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  tree <- as(treeObj, "phylo4")
  NODES <- c(names(ancestors(tree, sp)),sp)
  
  # Up 
  subtables <- grep(sp,tables, value=TRUE)
  subtables <- grep("_Up$",subtables, value=TRUE)
  HOG_species <- HOG %>% dplyr::filter(str_detect(Gene_id, paste0("^",sp,"_")))
  
  if(sp=="Medtru"){
    Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
    Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)
    HOG_species <- dplyr::left_join(HOG_species,Medtru_annot, by="Gene_id")
    Resume_HOGs_Up_species <- dplyr::filter(Resume_HOGs_Up,Node_label_Up %in% NODES)
    Resume_HOGs_Up_species <- Resume_HOGs_Up_species %>% dplyr::select(HOG,Category,Description,Node_label_Up,NodFactor,Organogenesis,Release,NFix)
    Resume_HOGs_Up_species <- dplyr::rename(Resume_HOGs_Up_species, NodFactor_Up = NodFactor, Organogenesis_Up = Organogenesis,
                                     Release_Up = Release,NFix_Up=NFix)
  }else{
    Resume_HOGs_Up_species <- dplyr::filter(Resume_HOGs_Up,Node_label_Up %in% NODES)
    Resume_HOGs_Up_species <- Resume_HOGs_Up_species %>% dplyr::select(HOG,Acronym,Category,Description,Node_label_Up,NodFactor,Organogenesis,Release,NFix)
    Resume_HOGs_Up_species <- dplyr::rename(Resume_HOGs_Up_species, NodFactor_Up = NodFactor, Organogenesis_Up = Organogenesis,
                                     Release_Up = Release,NFix_Up=NFix)
  }
  
  HOG_species <- dplyr::left_join(HOG_species,Resume_HOGs_Up_species, by="HOG")
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    DEGs$DEGs <- "Up"
    names(DEGs) <- c("Gene_id",gsub("./rnaseq_timepoints/","",i))
    
    HOG_species <- dplyr::left_join(HOG_species,DEGs, by="Gene_id")
    
  }
  
  # Down 
  subtables <- grep(sp,tables, value=TRUE)
  subtables <- grep("_Down$",subtables, value=TRUE)

  Resume_HOGs_Down_species <- dplyr::filter(Resume_HOGs_Down,Node_label_Down %in% NODES)
  Resume_HOGs_Down_species <- Resume_HOGs_Down_species %>% dplyr::select(HOG,Node_label_Down,NodFactor,Organogenesis,Release,NFix)
  Resume_HOGs_Down_species <- dplyr::rename(Resume_HOGs_Down_species, NodFactor_Down = NodFactor, Organogenesis_Down = Organogenesis,
                                   Release_Down = Release,NFix_Down=NFix)
  
  HOG_species <- dplyr::left_join(HOG_species,Resume_HOGs_Down_species, by="HOG")
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    DEGs$DEGs <- "Down"
    names(DEGs) <- c("Gene_id",gsub("./rnaseq_timepoints/","",i))
    
    HOG_species <- dplyr::left_join(HOG_species,DEGs, by="Gene_id")
    
  }
  
  HOG_species[is.na(HOG_species)] <- ""
  fwrite(HOG_species,paste0("./results/",sp,"_Nod_final_table.txt"), sep="\t",quote=FALSE)
  
}


