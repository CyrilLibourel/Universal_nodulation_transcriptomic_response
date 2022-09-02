######################
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ComplexHeatmap","circlize"), require, character.only = TRUE)

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")

##### Up ancestral nodes ##### 
Nodes_Up <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZonesUp.tsv"))
Nodes_Up <- Nodes_Up %>% tidyr::unite("NodFactor",c("Mimpud_NodFactor","Lotjap_NodFactor","Medtru_NodFactor"), sep="", remove=FALSE)
Nodes_Up$NodFactor[Nodes_Up$NodFactor=="UpUpUp"|Nodes_Up$NodFactor=="UpUp"] <- "NodFactor"
Nodes_Up$NodFactor[Nodes_Up$NodFactor!="NodFactor"] <- ""
Nodes_Up$Organogenesis[Nodes_Up$Mimpud_Organogenesis=="Up"]<- "Organogenesis"

Nodes_Up <- Nodes_Up %>% tidyr::unite("NFix",c("Mimpud_NFix","ZIII"), sep="", remove=FALSE)
Nodes_Up$NFix[Nodes_Up$NFix=="UpUp"] <- "NFix"
Nodes_Up$NFix[Nodes_Up$NFix!="NFix"] <- ""

Nodes_Up <- Nodes_Up %>% tidyr::unite("MimpudIntraCell",c("Mimpud_Release","Mimpud_Persist"), sep="", remove=FALSE)
Nodes_Up$MimpudIntraCell[Nodes_Up$MimpudIntraCell=="UpUp"] <- "Up"
Nodes_Up <- Nodes_Up %>% tidyr::unite("Release",c("MimpudIntraCell","FIId"), sep="", remove=FALSE)
Nodes_Up$Release[Nodes_Up$Release=="UpUp"] <- "Release"
Nodes_Up$Release[Nodes_Up$Release!="Release"] <- ""

Nodes_Up[is.na(Nodes_Up)] <- ""

Nodes_Up <- Nodes_Up %>% dplyr::select("HOG","Acronym","Category","Description","Mean_seq_number_HOG","Mean_seq_number_NodspecieswithDEGs",
                                       "Percentage_of_Nodspecies_in_HOG","Percentage_of_nonNodspecies_in_HOG","Node_label_Up","Ranking_Up",
                                       "Aeseve_Nod_FDR005_logFC0.5_Up","Arahyp_Nod_FDR005_logFC1.5_Up","Datglo_Nod_FDR005_logFC1_Up",
                                       "Hiprha_Nod_FDR005_logFC0.5_Up","Lotjap_Nod_FDR005_logFC1.5_Up","Lupalb_Nod_FDR005_logFC1_Up",
                                       "Medtru_Nod_FDR005_logFC1.5_Up","Mimpud_Nod_FDR005_logFC1.5_Up","Parand_Nod_FDR005_logFC1_Up",
                                       "Sampling_percentage","Mimpud_NF_inhib_def","Mimpud_NodFactor","Mimpud_Organogenesis","MimpudIntraCell","Mimpud_Release",
                                       "Mimpud_Persist","Mimpud_NFix","Lotjap_NodFactor","Medtru_NodFactor","FI","FIIp","FIId","IZ","ZIII",
                                       "NodFactor","Organogenesis","Release","NFix")

fwrite(Nodes_Up,paste0("./results/Nod_final_table_Up.txt"),sep="\t")

NODES <- c("N1","N2","N7","N12","N21","N22")
Classes <- c("NodFactor","Organogenesis","Release","NFix")

for(i in Classes){
  specific_class <- Nodes_Up[Nodes_Up$Node_label_Up=="N1",]
  # specific_class <- specific_class %>% arrange(Acronym)
  specific_class <- eval(parse(text = paste0("specific_class %>% dplyr::filter(",i,"=='",i,"')")))

  subset_species <- dplyr::select(specific_class,HOG,ends_with("_Up"))
  subset_species <- dplyr::select(subset_species,-Node_label_Up,-Ranking_Up)
  
  subset_species <- subset_species %>% tidyr::unite("PASTE",2:ncol(subset_species), sep="", remove=FALSE)
  subset_species$PASTE <- gsub("ns","",subset_species$PASTE)
  
  subset_species <- arrange(subset_species, desc(PASTE), HOG)
  
  Gene_nickname <- dplyr::left_join(subset_species,specific_class, by="HOG") %>% dplyr::select(HOG,Acronym, Category)
  subset_species <- subset_species %>% dplyr::select(-HOG,-PASTE)
  subset_species[subset_species == ""] <- "missing"
  
  SPECIES <- data.frame(class=names(subset_species))
  SPECIES <- tidyr::separate(SPECIES,class,into=c("Species","Nod","FDR","logFC","UpDown"),sep = "_")
  SPECIES <- unique(SPECIES$Species)
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  HOGs_tree <- keep.tip(treeObj,SPECIES)
  names(subset_species) <- SPECIES
  subset_species <- subset_species %>% dplyr::select(all_of(HOGs_tree$tip.label))
  
  
  dat <- as.matrix(subset_species)
  rownames(dat) <- Gene_nickname$Acronym
  OG_colors <- c("missing" = "white", "ns" = "cornsilk2", "Up" = "firebrick")
  
  ht1 <- Heatmap(dat, name = "HOG", cluster_rows = TRUE,rect_gp = gpar(col = "white", lwd = 0.2),
                 row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
                 cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
                 row_names_side = "left", col = OG_colors , width = unit(4, "cm"), height = unit(nrow(subset_species)/3, "cm"),row_title_rot = 0)

  pdffile = paste0("./results/Heatmap_NFN_Node_",i,"_Up.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=nrow(subset_species)/3,width=15)
  # par(mar=c(5,5,5,5))
  draw(ht1)
  dev.off()
  
}


##### Down ancestral nodes ##### 
Nodes_Down <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZonesDown.tsv"))
Nodes_Down <- Nodes_Down %>% tidyr::unite("NodFactor",c("Mimpud_NodFactor","Lotjap_NodFactor","Medtru_NodFactor"), sep="", remove=FALSE)
Nodes_Down$NodFactor[Nodes_Down$NodFactor=="DownDownDown"|Nodes_Down$NodFactor=="DownDown"] <- "NodFactor"
Nodes_Down$NodFactor[Nodes_Down$NodFactor!="NodFactor"] <- ""
Nodes_Down$Organogenesis[Nodes_Down$Mimpud_Organogenesis=="Down"]<- "Organogenesis"

Nodes_Down <- Nodes_Down %>% tidyr::unite("NFix",c("Mimpud_NFix","ZIII"), sep="", remove=FALSE)
Nodes_Down$NFix[Nodes_Down$NFix=="DownDown"] <- "NFix"
Nodes_Down$NFix[Nodes_Down$NFix!="NFix"] <- ""

Nodes_Down <- Nodes_Down %>% tidyr::unite("MimpudIntraCell",c("Mimpud_Release","Mimpud_Persist"), sep="", remove=FALSE)
Nodes_Down$MimpudIntraCell[Nodes_Down$MimpudIntraCell=="DownDown"] <- "Down"
Nodes_Down <- Nodes_Down %>% tidyr::unite("Release",c("MimpudIntraCell","FIId"), sep="", remove=FALSE)
Nodes_Down$Release[Nodes_Down$Release=="DownDown"] <- "Release"
Nodes_Down$Release[Nodes_Down$Release!="Release"] <- ""

Nodes_Down[is.na(Nodes_Down)] <- ""
Nodes_Down <- Nodes_Down %>% dplyr::select("HOG","Acronym","Category","Description","Mean_seq_number_HOG","Mean_seq_number_NodspecieswithDEGs",
                                       "Percentage_of_Nodspecies_in_HOG","Percentage_of_nonNodspecies_in_HOG","Node_label_Down","Ranking_Down",
                                       "Aeseve_Nod_FDR005_logFC0.5_Down","Arahyp_Nod_FDR005_logFC1.5_Down","Datglo_Nod_FDR005_logFC1_Down",
                                       "Hiprha_Nod_FDR005_logFC0.5_Down","Lotjap_Nod_FDR005_logFC1.5_Down","Lupalb_Nod_FDR005_logFC1_Down",
                                       "Medtru_Nod_FDR005_logFC1.5_Down","Mimpud_Nod_FDR005_logFC1.5_Down","Parand_Nod_FDR005_logFC1_Down",
                                       "Sampling_percentage","Mimpud_NF_inhib_def","Mimpud_NodFactor","Mimpud_Organogenesis","MimpudIntraCell","Mimpud_Release",
                                       "Mimpud_Persist","Mimpud_NFix","Lotjap_NodFactor","Medtru_NodFactor","FI","FIIp","FIId","IZ","ZIII",
                                       "NodFactor","Organogenesis","Release","NFix")

fwrite(Nodes_Down,paste0("./results/Nod_final_table_Down.txt"),sep="\t")


NODES <- c("N1","N2","N7","N12","N21","N22")
Classes <- c("NodFactor","Organogenesis","Release","NFix")

for(i in Classes){
  specific_class <- Nodes_Down[Nodes_Down$Node_label_Down=="N1",]
  # specific_class <- specific_class %>% arrange(Acronym)
  specific_class <- eval(parse(text = paste0("specific_class %>% dplyr::filter(",i,"=='",i,"')")))
  
  subset_species <- dplyr::select(specific_class,HOG,ends_with("_Down"))
  subset_species <- dplyr::select(subset_species,-Node_label_Down,-Ranking_Down)
  
  subset_species <- subset_species %>% tidyr::unite("PASTE",2:ncol(subset_species), sep="", remove=FALSE)
  subset_species$PASTE <- gsub("ns","",subset_species$PASTE)
  
  subset_species <- arrange(subset_species, desc(PASTE), HOG)
  
  Gene_nickname <- dplyr::left_join(subset_species,specific_class, by="HOG") %>% dplyr::select(HOG,Acronym, Category)
  subset_species <- subset_species %>% dplyr::select(-HOG,-PASTE)
  subset_species[subset_species == ""] <- "missing"
  
  SPECIES <- data.frame(class=names(subset_species))
  SPECIES <- tidyr::separate(SPECIES,class,into=c("Species","Nod","FDR","logFC","UpDown"),sep = "_")
  SPECIES <- unique(SPECIES$Species)
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  HOGs_tree <- keep.tip(treeObj,SPECIES)
  names(subset_species) <- SPECIES
  subset_species <- subset_species %>% dplyr::select(all_of(HOGs_tree$tip.label))
  
  
  dat <- as.matrix(subset_species)
  rownames(dat) <- Gene_nickname$Acronym
  OG_colors <- c("missing" = "white", "ns" = "cornsilk2", "Down" = "royalblue")

  ht1 <- Heatmap(dat, name = "HOG", cluster_rows = TRUE,rect_gp = gpar(col = "white", lwd = 0.2),
                 row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
                 cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
                 row_names_side = "left", col = OG_colors , width = unit(4, "cm"), height = unit(nrow(subset_species)/3, "cm"),row_title_rot = 0)


  pdffile = paste0("./results/Heatmap_NFN_Node_",i,"_Down.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=nrow(subset_species)/3,width=15)
  # par(mar=c(5,5,5,5))
  draw(ht1)
  dev.off()

}

