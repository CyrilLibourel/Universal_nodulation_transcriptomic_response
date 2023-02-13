######################
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","DBI","grid","mixOmics","DECIPHER","Biostrings","gridExtra",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ComplexHeatmap","circlize","dbplyr","wesanderson"), library, character.only = TRUE)

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")
sqlite_db_out <- "C:/Users/cyril.libourel/OneDrive/mimosaProject/work/Sql_db.sqlite"
outcon <- dbConnect(RSQLite::SQLite(), dbname=sqlite_db_out)
src_dbi(outcon)

##### Write ancestral nodes recruited in SQL db ##### 
# Nodes_Up <- dbGetQuery(outcon, paste0("select * from ","HOG_by_node_asrmkER_Up"))
# Nodes_Down <- dbGetQuery(outcon, paste0("select * from ","HOG_by_node_asrmkER_Down"))

##### Write species log fold changes in SQL db ##### 
Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

tables <- data.frame(Files=paste0(list.files(path="./Intraspe_tables/",pattern=".txt$",full.names=FALSE)))
tables$Files <- gsub(".txt$","",tables$Files)
tables$Files <- gsub("DEG_","",tables$Files)
tables = tidyr::separate(tables,Files,c("species","condition"),remove=FALSE,sep="_")

Species <- unique(tables$species)

HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))
HOG <- tidyr::separate(HOG,Gene_id,c("species"),remove=FALSE,sep="_")

Full_HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected"))

for(sp in Species){
  species_table <- tables %>% dplyr::filter(species==sp)
  upsetdf <- data.frame(HOG %>% dplyr::filter(species==sp) %>% dplyr::select(Gene_id,HOG))
  
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  tree <- as(treeObj, "phylo4")
  NODES <- c(names(ancestors(tree, sp)),sp)
  
  if(sp=="Medtru"){
    Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
    Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)
    upsetdf <- dplyr::left_join(upsetdf,Medtru_annot, by="Gene_id")
    Resume_HOGs_species <- Full_HOG
    Resume_HOGs_species <- Resume_HOGs_species %>% dplyr::select(HOG,Category,Description,starts_with(paste0(sp,"_WSR_")),
                                                                       Node_label_Up,Ranking_Up,Sampling_percentage_Up,
                                                                       Node_label_Down,Ranking_Down,Sampling_percentage_Down)
  }else{
    Resume_HOGs_species <- Full_HOG
    Resume_HOGs_species <- Resume_HOGs_species %>% dplyr::select(HOG,Acronym,Category,Description,starts_with(paste0(sp,"_WSR_")),
                                                                       Node_label_Up,Ranking_Up,Sampling_percentage_Up,
                                                                       Node_label_Down,Ranking_Down,Sampling_percentage_Down)
  }
  
  upsetdf <- dplyr::left_join(upsetdf,Resume_HOGs_species,by="HOG")
  # upsetdf <- dplyr::filter(upsetdf,Node_label_Up %in% NODES | Node_label_Down %in% NODES)
  upsetdf <- distinct(upsetdf)
  
  for(back in unique(species_table$condition)){
    spfile <- species_table %>% dplyr::filter(condition==back)
    sp_data <- fread(paste0("./Intraspe_tables/DEG_",spfile$Files,".txt"),h=TRUE)
    sp_data <- sp_data %>% dplyr::select(Gene_id,starts_with("logFC_"),starts_with("FDR"))
    
    upsetdf <- full_join(upsetdf,sp_data,by="Gene_id")
    
  }

  ips_tsv <- fread(paste0("../interproscan/IPR/IPR_aggregate/Table_Gene_id_IPR_aggregate_",sp,".txt"), h=TRUE, sep="\t", fill=TRUE, quote="")
  ips_tsv$Gene_id <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",ips_tsv$Gene_id)

  upsetdf <- left_join(upsetdf,ips_tsv,by="Gene_id")
  upsetdf$IPR[is.na(upsetdf$IPR)] <- ""
  
  upsetdf <- upsetdf %>% dplyr::select(Gene_id,HOG, Acronym,Category, Description,starts_with(paste0(sp,"_WSR_")),Node_label_Up,Node_label_Down,starts_with("logFC_"),starts_with("FDR"), IPR)

  dbWriteTable(outcon, paste0(sp,"_intra_DEGs"), upsetdf, overwrite=TRUE)
  message(sp," added in the SQL database")
}

Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))
HOG <- tidyr::separate(HOG,Gene_id,c("species"),remove=FALSE,sep="_")

for(i in species){
  species_DEGs <- dbGetQuery(outcon, paste0("select * from ",i,"_intra_DEGs"))
  
  fwrite(species_DEGs,paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/",i,"_intra_DEGs.txt"))
}

#### Add Up/Down data ####
Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

tables <- paste0(list.files(path="./rnaseq_timepoints/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))
HOG <- tidyr::separate(HOG,Gene_id,c("species"),remove=FALSE,sep="_")

Full_HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected"))


for(sp in Species){
  subtables <- grep(sp,tables, value=TRUE)
  subtables <- grep("_Up$",subtables, value=TRUE)
  
  upsetdf <- data.frame(HOG %>% dplyr::filter(species==sp) %>% dplyr::select(Gene_id,HOG))
  
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  tree <- as(treeObj, "phylo4")
  NODES <- c(names(ancestors(tree, sp)),sp)
  
  if(sp=="Medtru"){
    Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
    Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)
    upsetdf <- dplyr::left_join(upsetdf,Medtru_annot, by="Gene_id")
    Resume_HOGs_species <- Full_HOG
    Resume_HOGs_species <- Resume_HOGs_species %>% dplyr::select(HOG,Category,Description,
                                                                 Node_label_Up,Ranking_Up,Sampling_percentage_Up,
                                                                 Node_label_Down,Ranking_Down,Sampling_percentage_Down)
  }else{
    Resume_HOGs_species <- Full_HOG
    Resume_HOGs_species <- Resume_HOGs_species %>% dplyr::select(HOG,Acronym,Category,Description,
                                                                 Node_label_Up,Ranking_Up,Sampling_percentage_Up,
                                                                 Node_label_Down,Ranking_Down,Sampling_percentage_Down)
  }
  
  upsetdf <- dplyr::left_join(upsetdf,Resume_HOGs_species,by="HOG")
  upsetdf <- distinct(upsetdf)
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    DEGs$DEGs <- "Up"
    names(DEGs) <- c("Gene_id",gsub("./rnaseq_timepoints/","",i))
    
    upsetdf <- dplyr::left_join(upsetdf,DEGs, by="Gene_id")
    
  }
  
  #Down
  subtables <- grep(sp,tables, value=TRUE)
  subtables <- grep("_Down$",subtables, value=TRUE)
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    DEGs$DEGs <- "Down"
    names(DEGs) <- c("Gene_id",gsub("./rnaseq_timepoints/","",i))
    
    upsetdf <- dplyr::left_join(upsetdf,DEGs, by="Gene_id")
    
  }
  
  ips_tsv <- fread(paste0("../interproscan/IPR/IPR_aggregate/Table_Gene_id_IPR_aggregate_",sp,".txt"), h=TRUE, sep="\t", fill=TRUE, quote="")
  ips_tsv$Gene_id <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",ips_tsv$Gene_id)
  
  upsetdf <- left_join(upsetdf,ips_tsv,by="Gene_id")
  upsetdf[is.na(upsetdf)] <- ""

  
  SP_HOG <- data.frame(HOG %>% dplyr::filter(species==sp) %>% dplyr::select(Gene_id,HOG)); SP_HOG <- unique(SP_HOG$HOG)
  HOG_species <- data.frame(HOG %>% dplyr::filter(species!=sp) %>% dplyr::select(Gene_id,HOG)); HOG_species <- unique(HOG_species$HOG)
  SP_specific_HOG <- data.frame(HOG=unique(setdiff(SP_HOG,HOG_species)),specific=paste0(sp,"_specific"))
  
  upsetdf$New_Node_label_Up <- upsetdf[,grep(paste0("_WSR_Nod_Up"),names(upsetdf))]
  upsetdf$New_Node_label_Down <- upsetdf[,grep(paste0("_WSR_Nod_Down"),names(upsetdf))]

  upsetdf$Node_label_Up[upsetdf$HOG=="" & upsetdf$New_Node_label_Up=="Up"] <- paste0(sp,"_specific")
  upsetdf$Node_label_Down[upsetdf$HOG=="" & upsetdf$New_Node_label_Down=="Down"] <- paste0(sp,"_specific")
  
  for(hogs in SP_specific_HOG$HOG){
    upsetdf$Node_label_Up[upsetdf$HOG==hogs & upsetdf$Node_label_Up==sp] <- paste0(sp,"_specific")
    upsetdf$Node_label_Down[upsetdf$HOG==hogs & upsetdf$Node_label_Down==sp] <- paste0(sp,"_specific")
    
  }
  
  PresAbs_HOG <- dbGetQuery(outcon, paste0("select * from ","PresenceAbsence_HOG_node_asrmkER"))
  PresAbs_HOG <- PresAbs_HOG %>% dplyr::select(HOG,Node_label)
  PresAbs_HOG <- dplyr::filter(PresAbs_HOG,Node_label %in% NODES);names(PresAbs_HOG) <- c("HOG","HOG_node_label")
  PresAbs_HOG <- distinct(PresAbs_HOG)
  
  upsetdf <- dplyr::left_join(upsetdf,PresAbs_HOG,by="HOG")
  upsetdf <- upsetdf %>% dplyr::select(-New_Node_label_Up,-New_Node_label_Down)

  dbWriteTable(outcon, paste0(sp,"_UpDown"), upsetdf, overwrite=TRUE)
  message(sp," added in the SQL database")
  
}


Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

for(i in Species){
  species_DEGs <- dbGetQuery(outcon, paste0("select * from ",i,"_UpDown"))
  fwrite(species_DEGs,paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/Nature_Plants/",i,"_intra_UpDown.txt"),quote=FALSE, sep="\t")
}


#### Count number of genes per species ####

Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

species_node_summary_up <- data.frame(Node_label=c("N1","N2","N3","N7","N12","N15","N21","N22",Species,paste0(Species,"_specific")))
species_node_summary_down <- data.frame(Node_label=c("N1","N2","N3","N7","N12","N15","N21","N22",Species,paste0(Species,"_specific")))

for(i in Species){
  species_DEGs <- dbGetQuery(outcon, paste0("select * from ",i,"_UpDown"))
  NODES <- c(names(ancestors(tree, i)),i, paste0(i,"_specific"))
  
  species_DEGs_Up <- species_DEGs %>% dplyr::filter(Node_label_Up %in% NODES) %>% dplyr::select(Gene_id,Node_label_Up, starts_with(paste0(i,"_WSR_Nod_Up"))) %>%
    dplyr::filter(.[[3]]=="Up")
  summary_species_DEGs_Up <- species_DEGs_Up %>% group_by(Node_label_Up) %>% summarise(Count = length(unique(Gene_id)))
  names(summary_species_DEGs_Up) <- c("Node_label",i)
  
  species_node_summary_up <- dplyr::left_join(species_node_summary_up,summary_species_DEGs_Up,by="Node_label")
  
  species_DEGs_Down <- species_DEGs %>% dplyr::filter(Node_label_Down %in% NODES) %>% dplyr::select(Gene_id,Node_label_Down, starts_with(paste0(i,"_WSR_Nod_Down"))) %>%
    dplyr::filter(.[[3]]=="Down")
  species_DEGs_Down <- species_DEGs_Down %>% group_by(Node_label_Down) %>% summarise(Count = length(unique(Gene_id)))
  names(species_DEGs_Down) <- c("Node_label",i)
  species_node_summary_down <- dplyr::left_join(species_node_summary_down,species_DEGs_Down,by="Node_label")

}
species_node_summary_up[is.na(species_node_summary_up)] <- 0
species_node_summary_down[is.na(species_node_summary_down)] <- 0

dbWriteTable(outcon, paste0("Table_count_genes_by_node_species_asrmkER_Up"), species_node_summary_up, overwrite=TRUE)
dbWriteTable(outcon, paste0("Table_count_genes_by_node_species_asrmkER_Down"), species_node_summary_down, overwrite=TRUE)

#### Count number of genes per traits ####

Species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

species_node_summary_up <- data.frame(Node_label=c("N1","N2","N3","N7","N12","N15","N21","N22",Species,paste0(Species,"_specific")))
species_node_summary_down <- data.frame(Node_label=c("N1","N2","N3","N7","N12","N15","N21","N22",Species,paste0(Species,"_specific")))

for(i in Species){
  species_DEGs <- dbGetQuery(outcon, paste0("select * from ",i,"_UpDown"))
  NODES <- c(names(ancestors(tree, i)),i, paste0(i,"_specific"))
  COL <- grep(i,names(species_DEGs),value=TRUE)
  COL_Up <- grep("Up",COL,value=TRUE)
  for(j in COL_Up){
  species_DEGs_Up <- species_DEGs %>% dplyr::filter(Node_label_Up %in% NODES) %>% dplyr::select(Gene_id,Node_label_Up, starts_with(paste0(j))) %>%
    dplyr::filter(.[[3]]=="Up")
  summary_species_DEGs_Up <- species_DEGs_Up %>% group_by(Node_label_Up) %>% summarise(Count = length(unique(Gene_id)))
  names(summary_species_DEGs_Up) <- c("Node_label",j)
  
  species_node_summary_up <- dplyr::left_join(species_node_summary_up,summary_species_DEGs_Up,by="Node_label")
  }
  
  COL_Down <- grep("Down",COL,value=TRUE)
  for(j in COL_Down){
  species_DEGs_Down <- species_DEGs %>% dplyr::filter(Node_label_Down %in% NODES) %>% dplyr::select(Gene_id,Node_label_Down, starts_with(paste0(j))) %>%
    dplyr::filter(.[[3]]=="Down")
  species_DEGs_Down <- species_DEGs_Down %>% group_by(Node_label_Down) %>% summarise(Count = length(unique(Gene_id)))
  names(species_DEGs_Down) <- c("Node_label",j)
  species_node_summary_down <- dplyr::left_join(species_node_summary_down,species_DEGs_Down,by="Node_label")
  }
  
}
species_node_summary_up[is.na(species_node_summary_up)] <- 0
species_node_summary_down[is.na(species_node_summary_down)] <- 0

dbWriteTable(outcon, paste0("Table_count_genes_by_conditions_Up"), species_node_summary_up, overwrite=TRUE)
dbWriteTable(outcon, paste0("Table_count_genes_by_conditions_Down"), species_node_summary_down, overwrite=TRUE)



#### Add Up/Down traits data to N0_corrected ####
tables <- paste0(list.files(path="./rnaseq_other/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

Full_HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected"))


  subtables <- grep("_Up$",tables, value=TRUE)
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    names(DEGs) <- c("Gene_id")
    subupsetdf <- dplyr::left_join(DEGs,HOG, by="Gene_id") %>% dplyr::filter(HOG!="")
    subupsetdf <- data.frame(HOG=unique(subupsetdf$HOG),"V2"="Up")
    names(subupsetdf) <- c("HOG",gsub("./rnaseq_other/","",i))
    Full_HOG <- dplyr::left_join(Full_HOG,subupsetdf, by="HOG")
    
  }
  
  Full_HOG[is.na(Full_HOG)] <- ""
  
  subtables <- grep("_Down$",tables, value=TRUE)
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    DEGs <- DEGs %>% distinct()
    names(DEGs) <- c("Gene_id")
    subupsetdf <- dplyr::left_join(DEGs,HOG, by="Gene_id") %>% dplyr::filter(HOG!="")
    subupsetdf <- data.frame(HOG=unique(subupsetdf$HOG),"V2"="Down")
    names(subupsetdf) <- c("HOG",gsub("./rnaseq_other/","",i))
    Full_HOG <- dplyr::left_join(Full_HOG,subupsetdf, by="HOG")
    
  }
  
  Full_HOG[is.na(Full_HOG)] <- ""
  
dbWriteTable(outcon, paste0("N0_corrected_with_traits"), Full_HOG, overwrite=TRUE)

#### Figure 2 AMS Lateral Ancestral Node ####
Full_HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_with_traits"))

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

Subset_Full_HOG <- Full_HOG %>% dplyr::select(HOG,Acronym,Category,Node_label_Up,Lotjap_Myc_Up,Medtru_Myc_Up,Medtru_Lateralroot_Up)
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::filter(Acronym!="")
Subset_Full_HOG$Node_label_Up[Subset_Full_HOG$Node_label_Up!="N1"] <- ""
Subset_Full_HOG$Node_label_Up[Subset_Full_HOG$Node_label_Up=="N1"] <- "Up"
Myc <- Subset_Full_HOG %>% tidyr::unite("Myc",Lotjap_Myc_Up,Medtru_Myc_Up, sep="", remove=FALSE) %>% dplyr::select(Myc)
Subset_Full_HOG$Myc <- Myc$Myc
Subset_Full_HOG$Myc[Subset_Full_HOG$Myc!=""] <- "Up"
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::select(-Lotjap_Myc_Up,-Medtru_Myc_Up)
Subset_Full_HOG <- dplyr::distinct(Subset_Full_HOG)
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::filter(Node_label_Up=="Up" | Medtru_Lateralroot_Up=="Up" | Myc=="Up" )
Subset_Full_HOG[Subset_Full_HOG==""] <- "ns"
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::arrange(Acronym)

fwrite(Subset_Full_HOG,"C:/Users/cyril.libourel/OneDrive/mimosaProject/work/Heatmap_Fig2_RNS_AMS_Latroot_Up.txt",quote=FALSE, sep="\t")

dat <- as.matrix(Subset_Full_HOG %>% dplyr::select(Node_label_Up,Myc,Medtru_Lateralroot_Up))
rownames(dat) <- Subset_Full_HOG$Acronym
OG_colors <- c("ns" = "cornsilk2", "Up" = "firebrick")

split <- Subset_Full_HOG %>% tidyr::unite("PASTE",Node_label_Up,Myc,Medtru_Lateralroot_Up, sep="", remove=FALSE) %>% dplyr::select(PASTE)

ht1 <- Heatmap(dat, name = "HOG", cluster_rows = TRUE,rect_gp = gpar(col = "white", lwd = 0.2),
               row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
               cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
               row_names_side = "left", col = OG_colors , width = unit(2, "cm"), height = unit(nrow(Subset_Full_HOG)/3, "cm"),row_title_rot = 0,split = split)

pdffile = paste0("./Results_paper/Heatmap_Fig2_RNS_AMS_Latroot_Up.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=nrow(Subset_Full_HOG)/2,width=10)
# par(mar=c(5,5,5,5))
draw(ht1)
dev.off()

Full_HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_with_traits"))

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

Subset_Full_HOG <- Full_HOG %>% dplyr::select(HOG,Acronym,Category,Node_label_Down,Lotjap_Myc_Down,Medtru_Myc_Down,Medtru_Lateralroot_Down)
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::filter(Acronym!="")
Subset_Full_HOG$Node_label_Down[Subset_Full_HOG$Node_label_Down!="N1"] <- ""
Subset_Full_HOG$Node_label_Down[Subset_Full_HOG$Node_label_Down=="N1"] <- "Down"
Myc <- Subset_Full_HOG %>% tidyr::unite("Myc",Lotjap_Myc_Down,Medtru_Myc_Down, sep="", remove=FALSE) %>% dplyr::select(Myc)
Subset_Full_HOG$Myc <- Myc$Myc
Subset_Full_HOG$Myc[Subset_Full_HOG$Myc!=""] <- "Down"
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::select(-Lotjap_Myc_Down,-Medtru_Myc_Down)
Subset_Full_HOG <- dplyr::distinct(Subset_Full_HOG)
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::filter(Node_label_Down=="Down" | Medtru_Lateralroot_Down=="Down" | Myc=="Down" )
Subset_Full_HOG[Subset_Full_HOG==""] <- "ns"
Subset_Full_HOG <- Subset_Full_HOG %>% dplyr::arrange(Acronym)

fwrite(Subset_Full_HOG,"C:/Users/cyril.libourel/OneDrive/mimosaProject/work/Heatmap_Fig2_RNS_AMS_Latroot_Down.txt",quote=FALSE, sep="\t")

dat <- as.matrix(Subset_Full_HOG %>% dplyr::select(Node_label_Down,Myc,Medtru_Lateralroot_Down))
rownames(dat) <- Subset_Full_HOG$Acronym
OG_colors <- c("ns" = "cornsilk2", "Down" = "royalblue")

split <- Subset_Full_HOG %>% tidyr::unite("PASTE",Node_label_Down,Myc,Medtru_Lateralroot_Down, sep="", remove=FALSE) %>% dplyr::select(PASTE)

ht1 <- Heatmap(dat, name = "HOG", cluster_rows = TRUE,rect_gp = gpar(col = "white", lwd = 0.2),
               row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
               cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
               row_names_side = "left", col = OG_colors , width = unit(2, "cm"), height = unit(nrow(Subset_Full_HOG)/3, "cm"),row_title_rot = 0,split = split)

pdffile = paste0("./Results_paper/Heatmap_Fig2_RNS_AMS_Latroot_Down.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=nrow(Subset_Full_HOG)/2,width=10)
# par(mar=c(5,5,5,5))
draw(ht1)
dev.off()

#### Analysis traits multi species ####
#### Analysis of mutants vs 5dpi_niroot Figure 3 ####
treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")
NODES <- c(names(ancestors(tree, "Mimpud")),"Mimpud")

# Up 
Mimpud_traits <- dbGetQuery(outcon, paste0("select * from ","Mimpud_UpDown"))
Mimpud_traits <- Mimpud_traits %>% dplyr::select(Gene_id, ends_with("_Up")) 
Mimpud_traits_shared <- Mimpud_traits %>% dplyr::filter(Mimpud_WSR_Nod_Up=="Up") %>% dplyr::filter(Node_label_Up %in% NODES | Node_label_Up=="Mimpud_specific")

Full_summary <- data.frame()
Clones <- c("RsolGMI1000","pRalta","hrcV","hrpG","efpR","nifH","WSR")
All_DEGs <- data.frame()

for(i in Clones){
  subclones <- Mimpud_traits_shared %>% dplyr::select(Gene_id,starts_with(paste0("Mimpud_",i)),starts_with(paste0("Mimpud_WSR_")),Node_label_Up) %>% 
    dplyr::filter_at(2,all_vars(.=="Up")) %>%
    group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))
  names(subclones) <- c("Node_label","Overlap")
  subclones$Sum <- sum(subclones$Overlap);subclones$Clones <- i
  All <- data.frame(Node_label_Up="All",Overlap = Mimpud_traits %>% dplyr::select(Gene_id,starts_with(paste0("Mimpud_",i))) %>% 
                      dplyr::filter_at(2,all_vars(.=="Up")) %>% dplyr::summarise(Overlap=length(unique(Gene_id))))
  subclones$Totalshared[subclones$Node_label=="N1"] <- All$Overlap
  Full_summary <- rbind(Full_summary,subclones)
}

Full_summary$Up_Down <- "Up"

# Down 
Mimpud_traits <- dbGetQuery(outcon, paste0("select * from ","Mimpud_UpDown"))
Mimpud_traits <- Mimpud_traits %>% dplyr::select(Gene_id, ends_with("_Down")) 
Mimpud_traits_shared <- Mimpud_traits %>% dplyr::filter(Mimpud_WSR_Nod_Down=="Down") %>% dplyr::filter(Node_label_Down %in% NODES | Node_label_Down=="Mimpud_specific")
Full_summary_Down <- data.frame()
Clones <- c("RsolGMI1000","pRalta","hrcV","hrpG","efpR","nifH","WSR")
All_DEGs_Down <- data.frame()

for(i in Clones){
  subclones <- Mimpud_traits_shared %>% dplyr::select(Gene_id,starts_with(paste0("Mimpud_",i)),starts_with(paste0("Mimpud_WSR_")),Node_label_Down) %>% 
    dplyr::filter_at(2,all_vars(.=="Down")) %>%
    group_by(Node_label_Down) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))
  names(subclones) <- c("Node_label","Overlap")
  subclones$Sum <- sum(subclones$Overlap);subclones$Clones <- i
  All <- data.frame(Node_label_Down="All",Overlap = Mimpud_traits %>% dplyr::select(Gene_id,starts_with(paste0("Mimpud_",i))) %>% 
                      dplyr::filter_at(2,all_vars(.=="Down")) %>% dplyr::summarise(Overlap=length(unique(Gene_id))))
  subclones$Totalshared[subclones$Node_label=="N1"] <- All$Overlap
  Full_summary_Down <- rbind(Full_summary_Down,subclones)

}

Full_summary_Down$Up_Down <- "Down"

Full_summary <- rbind(Full_summary,Full_summary_Down)

Clones <- c("RsolGMI1000","pRalta","hrcV","hrpG","efpR","nifH","WSR")

Nodes <- c("N1","N2","Mimpud","Mimpud_specific")
for(UpDown in c("Up","Down")){
for(Node in Nodes){
  for(i in Clones){
    CloneNode <- Full_summary$Overlap[Full_summary$Clones==i & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    CtaiNode <- Full_summary$Overlap[Full_summary$Clones=="WSR" & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])-CloneNode,
                                                  unique(Full_summary$Sum[Full_summary$Clones=="WSR" & Full_summary$Up_Down==UpDown])-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- test$p.value
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
    Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$estimate
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
  }
}
}
Full_summary$Total[Full_summary$Clones=="WSR" & Full_summary$Node_label=="N1"] <- unique(Full_summary$Sum[Full_summary$Clones=="WSR"])

{Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"

Full_summary$Node_label[Full_summary$Node_label=="Mimpud"] <- "Mimpud_DEOGs"
Full_summary$Node_label[Full_summary$Node_label=="Mimpud_specific"] <- "Mimpud_specific_DEOGs"
Full_summary$Clones[Full_summary$Clones=="WSR"] <- "C.tai"
Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","Mimpud_DEOGs","Mimpud_specific_DEOGs"))
Full_summary$Clones <- factor(Full_summary$Clones,levels = c("RsolGMI1000","pRalta","hrcV","hrpG","efpR","nifH","C.tai"))

Full_summary$Up_Down <- factor(Full_summary$Up_Down, levels=c("Up","Down"))
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than C.tai"}

LABEL_NODE <- c(`Up` = "Up",
                `Down` = "Down")

Full_summary <- Full_summary %>% tidyr::unite("Clone_Up_Down",c("Clones","Up_Down"), sep="_", remove=FALSE)

mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]


pdf(file = paste0("./Results_paper/Barplot_count_shared_with_Ctai_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
ggplot(data=Full_summary, aes(x=Clones, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="stack",  width = 0.65) +
  scale_fill_manual(values = mycol,labels = c("NFN node", "Fabales node", "Mimpud_DEOGs","Mimpud_specific_DEOGs")) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
  labs(y= "Number of genes shared with C.tai", x = "Strains") +
  geom_text(data = Full_summary, aes(y = max(Total*1.05,na.rm=TRUE), label = Total)) +
  theme_classic(base_size = 16) +
  geom_point(aes(x=Clones, y=Totalshared, col="red")) +
  geom_text(data = Full_summary, aes(y = Totalshared, label = Totalshared))
dev.off()

pdf(file = paste0("./Results_paper/Barplot_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
ggplot(data=Full_summary, aes(x=Clones, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = mycol,labels = c("NFN node", "Fabales node", "Mimpud_DEOGs","Mimpud_specific_DEOGs")) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
  labs(y= "Number of genes shared with C.tai", x = "Strains") +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.5))+
  scale_colour_manual(values=c("coral2", "black")) +
  theme_classic(base_size = 16) + labs(fill = "Node") 
dev.off()


#### Figure traits (Figure 4) ####
treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")
MimpudNODES <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
LotjapNODES <- c(names(ancestors(tree, "Lotjap")),"Lotjap")
MedtruNODES <- c(names(ancestors(tree, "Medtru")),"Medtru")

#### Mimpud traits #### 
Mimpud_traits <- dbGetQuery(outcon, paste0("select * from ","Mimpud_UpDown"))
Mimpud_traits <- Mimpud_traits %>% dplyr::select(Gene_id, ends_with("_Up")) %>% dplyr::filter(Mimpud_WSR_Nod_Up=="Up") %>% dplyr::filter(Node_label_Up %in% MimpudNODES | Node_label_Up=="Mimpud_specific")

Full_summary <- data.frame()
Clones <- c("NF","Organogenesis","Release","Persist","NFix","WSR")

for(i in Clones){
  subclones <- Mimpud_traits %>% dplyr::select(Gene_id,starts_with(paste0("Mimpud_",i)),Node_label_Up) %>% 
    dplyr::filter_at(2,all_vars(.=="Up")) %>%
    group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))
  names(subclones) <- c("Node_label","Overlap")
  subclones$Sum <- sum(subclones$Overlap);subclones$Clones <- i
  Full_summary <- rbind(Full_summary,subclones)
  
}

Full_summary$Up_Down <- "Up"
Nodes <- c("N1","N2","Mimpud","Mimpud_specific")
UpDown<- "Up"
  for(Node in Nodes){
    for(i in Clones){
      CloneNode <- Full_summary$Overlap[Full_summary$Clones==i & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
      CtaiNode <- Full_summary$Overlap[Full_summary$Clones=="WSR" & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
      M <- as.table(rbind(c(CloneNode, CtaiNode), c(unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])-CloneNode,
                                                    unique(Full_summary$Sum[Full_summary$Clones=="WSR" & Full_summary$Up_Down==UpDown])-CtaiNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- test$p.value
      Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
      Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$estimate
      Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
    }
  }

Full_summary$Total[Full_summary$Clones=="WSR" & Full_summary$Node_label=="N1"] <- unique(Full_summary$Sum[Full_summary$Clones=="WSR"])

{Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  
  Full_summary$Node_label[Full_summary$Node_label=="Mimpud"] <- "Mimpud_DEOGs"
  Full_summary$Node_label[Full_summary$Node_label=="Mimpud_specific"] <- "Mimpud_specific_DEOGs"
  Full_summary$Clones[Full_summary$Clones=="WSR"] <- "C.tai"
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","Mimpud_DEOGs","Mimpud_specific_DEOGs"))
  Full_summary$Clones <- factor(Full_summary$Clones,levels = c("NF","Organogenesis","Release","Persist","NFix","C.tai"))
  
  Full_summary$Up_Down <- factor(Full_summary$Up_Down, levels=c("Up","Down"))
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than C.tai"}

LABEL_NODE <- c(`Up` = "Up",
                `Down` = "Down")

Full_summary <- Full_summary %>% tidyr::unite("Clone_Up_Down",c("Clones","Up_Down"), sep="_", remove=FALSE)
Mimpud_summary <- Full_summary

dbWriteTable(outcon, paste0("Mimpud_traits_table"), Mimpud_summary, overwrite=TRUE)


mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]



#### Medtru traits #### 
Medtru_traits <- dbGetQuery(outcon, paste0("select * from ","Medtru_UpDown"))
Medtru_traits <- Medtru_traits %>% dplyr::select(Gene_id, ends_with("_Up")) %>% dplyr::filter(Medtru_WSR_Nod_Up=="Up") %>% dplyr::filter(Node_label_Up %in% MedtruNODES | Node_label_Up=="Medtru_specific")
WLRR <- Medtru_traits %>% tidyr::unite("Medtru_WLRR",starts_with("Medtru_latroot"), sep="", remove=FALSE) %>% dplyr::select(Medtru_WLRR)
Medtru_traits$Medtru_WLRR <- WLRR$Medtru_WLRR
Medtru_traits$Medtru_WLRR[Medtru_traits$Medtru_WLRR!=""] <- "Up"
WMycR <- Medtru_traits %>% tidyr::unite("Medtru_WMycR",starts_with("Medtru_Myc"), sep="", remove=FALSE) %>% dplyr::select(Medtru_WMycR)
Medtru_traits$Medtru_WMycR <- WMycR$Medtru_WMycR
Medtru_traits$Medtru_WMycR[Medtru_traits$Medtru_WMycR!=""] <- "Up"

Full_summary <- data.frame()
Clones <- c("NF","FI","FIId","FIIp","IZ","ZIII","WLRR","WMycR","WSR")

for(i in Clones){
  subclones <- Medtru_traits %>% dplyr::select(Gene_id,starts_with(paste0("Medtru_",i)),Node_label_Up) %>% 
    dplyr::filter_at(2,all_vars(.=="Up")) %>%
    group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))
  names(subclones) <- c("Node_label","Overlap")
  subclones$Sum <- sum(subclones$Overlap);subclones$Clones <- i
  Full_summary <- rbind(Full_summary,subclones)
  
}

Full_summary$Up_Down <- "Up"
Nodes <- unique(Full_summary$Node_label)
UpDown<- "Up"
for(Node in Nodes){
  for(i in Clones){
    CloneNode <- Full_summary$Overlap[Full_summary$Clones==i & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    CtaiNode <- Full_summary$Overlap[Full_summary$Clones=="WSR" & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])-CloneNode,
                                                  unique(Full_summary$Sum[Full_summary$Clones=="WSR" & Full_summary$Up_Down==UpDown])-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- test$p.value
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
    Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$estimate
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
  }
}

Full_summary$Total[Full_summary$Clones=="WSR" & Full_summary$Node_label=="N1"] <- unique(Full_summary$Sum[Full_summary$Clones=="WSR"])

{Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  
  Full_summary$Node_label[Full_summary$Node_label=="Medtru"] <- "Medtru_DEOGs"
  Full_summary$Node_label[Full_summary$Node_label=="Medtru_specific"] <- "Medtru_specific_DEOGs"
  Full_summary$Clones[Full_summary$Clones=="WSR"] <- "Medtru"
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N7","N12","N21","Medtru_DEOGs","Medtru_specific_DEOGs"))
  Full_summary$Clones <- factor(Full_summary$Clones,levels = c("NF","FI","FIId","FIIp","IZ","ZIII","WLRR","WMycR","Medtru"))
  
  Full_summary$Up_Down <- factor(Full_summary$Up_Down, levels=c("Up","Down"))
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than Medtru"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than Medtru"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than Medtru"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than Medtru"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than Medtru"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than Medtru"}

LABEL_NODE <- c(`Up` = "Up",
                `Down` = "Down")

Full_summary <- Full_summary %>% tidyr::unite("Clone_Up_Down",c("Clones","Up_Down"), sep="_", remove=FALSE)
Medtru_summary <- Full_summary

dbWriteTable(outcon, paste0("Medtru_traits_table"), Medtru_summary, overwrite=TRUE)

#### Lotjap traits #### 
Lotjap_traits <- dbGetQuery(outcon, paste0("select * from ","Lotjap_UpDown"))
Lotjap_traits <- Lotjap_traits %>% dplyr::select(Gene_id, ends_with("_Up")) %>% dplyr::filter(Lotjap_WSR_Nod_Up=="Up") %>% dplyr::filter(Node_label_Up %in% LotjapNODES | Node_label_Up=="Lotjap_specific")


Full_summary <- data.frame()
Clones <- c("NF","Myc","WSR")

for(i in Clones){
  subclones <- Lotjap_traits %>% dplyr::select(Gene_id,starts_with(paste0("Lotjap_",i)),Node_label_Up) %>% 
    dplyr::filter_at(2,all_vars(.=="Up")) %>%
    group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))
  names(subclones) <- c("Node_label","Overlap")
  subclones$Sum <- sum(subclones$Overlap);subclones$Clones <- i
  Full_summary <- rbind(Full_summary,subclones)
  
}

Full_summary$Up_Down <- "Up"
Nodes <- unique(Full_summary$Node_label)
UpDown<- "Up"
for(Node in Nodes){
  for(i in Clones){
    CloneNode <- Full_summary$Overlap[Full_summary$Clones==i & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    CtaiNode <- Full_summary$Overlap[Full_summary$Clones=="WSR" & Full_summary$Node_label==Node & Full_summary$Up_Down==UpDown]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])-CloneNode,
                                                  unique(Full_summary$Sum[Full_summary$Clones=="WSR" & Full_summary$Up_Down==UpDown])-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- test$p.value
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
    Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- fishertest$estimate
    Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Clones==i & Full_summary$Up_Down==UpDown] <- unique(Full_summary$Sum[Full_summary$Clones==i & Full_summary$Up_Down==UpDown])
  }
}

Full_summary$Total[Full_summary$Clones=="WSR" & Full_summary$Node_label=="N1"] <- unique(Full_summary$Sum[Full_summary$Clones=="WSR"])

{Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  
  Full_summary$Node_label[Full_summary$Node_label=="Lotjap"] <- "Lotjap_DEOGs"
  Full_summary$Node_label[Full_summary$Node_label=="Lotjap_specific"] <- "Lotjap_specific_DEOGs"
  Full_summary$Clones[Full_summary$Clones=="WSR"] <- "Lotjap"
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N7","N12","N21","Lotjap_DEOGs","Lotjap_specific_DEOGs"))
  Full_summary$Clones <- factor(Full_summary$Clones,levels = c("NF","Myc","Lotjap"))
  
  Full_summary$Up_Down <- factor(Full_summary$Up_Down, levels=c("Up","Down"))
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than Lotjap"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than Lotjap"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than Lotjap"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than Lotjap"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than Lotjap"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than Lotjap"}

LABEL_NODE <- c(`Up` = "Up",
                `Down` = "Down")

Full_summary <- Full_summary %>% tidyr::unite("Clone_Up_Down",c("Clones","Up_Down"), sep="_", remove=FALSE)
Lotjap_summary <- Full_summary

dbWriteTable(outcon, paste0("Lotjap_traits_table"), Lotjap_summary, overwrite=TRUE)

#### Merging Tables Fig 4 ####
Mimpud_summary <- dbGetQuery(outcon, paste0("select * from ","Mimpud_traits_table"))
Mimpud_summary$Clones <- paste0("Mimpud_",Mimpud_summary$Clones)
Medtru_summary <- dbGetQuery(outcon, paste0("select * from ","Medtru_traits_table"))
Medtru_summary$Clones <- paste0("Medtru_",Medtru_summary$Clones)
Lotjap_summary <- dbGetQuery(outcon, paste0("select * from ","Lotjap_traits_table"))
Lotjap_summary$Clones <- paste0("Lotjap_",Lotjap_summary$Clones)

Full_summary <- rbind(Mimpud_summary,Medtru_summary)
Full_summary <- rbind(Full_summary,Lotjap_summary)
Full_summary[is.na(Full_summary)] <- ""

Full_summary$Clones[Full_summary$Clones=="Mimpud_C.tai"] <- "Mimpud"
Full_summary$Node_label[Full_summary$Node_label=="Mimpud_DEOGs"] <- "Species_DEOGs"
Full_summary$Node_label[Full_summary$Node_label=="Medtru_DEOGs"] <- "Species_DEOGs"
Full_summary$Node_label[Full_summary$Node_label=="Lotjap_DEOGs"] <- "Species_DEOGs"
Full_summary$Node_label[Full_summary$Node_label=="Mimpud_specific_DEOGs"] <- "Species_specific"
Full_summary$Node_label[Full_summary$Node_label=="Medtru_specific_DEOGs"] <- "Species_specific"
Full_summary$Node_label[Full_summary$Node_label=="Lotjap_specific_DEOGs"] <- "Species_specific"
Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N7","N12","N21","Species_DEOGs","Species_specific"))
Full_summary$Clones <- factor(Full_summary$Clones,
                             levels = c("Mimpud","Mimpud_NF","Mimpud_Organogenesis","Mimpud_Release","Mimpud_Persist","Mimpud_NFix",
                                        "Medtru_Medtru","Medtru_NF","Medtru_FI","Medtru_FIId","Medtru_FIIp","Medtru_IZ","Medtru_ZIII","Medtru_WLRR","Medtru_WMycR",
                                        "Lotjap_Lotjap","Lotjap_NF","Lotjap_Myc"))
Full_summary$Color <- "black"
Full_summary$Color[grep("More",Full_summary$Odds)] <- "white"

mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4", "Species_DEOGs" = "#E1AF00","Species_specific"="#F21A00")

pdf(file = paste0("./Results_paper/Barplot_Proportion_shared_with_Multispecies_release_multi_Nodes_Zones_Traits_Up.pdf"), width=20, height=10)
ggplot(data=Full_summary, aes(x=Clones, y=Overlap, fill=Node_label)) +
        geom_bar(stat="identity", position="fill",  width = 0.65) +
        scale_fill_manual(values = mycol,labels = c("NFN clade", "Fabales clade", "Papilionoideae clade",
                                                    "Dalbergioid+Hologalegina clade","Hologalegina clade", "Species_DEOGs","Species_specific")) + 
        labs(y= "Number of genes shared with species transcriptomic response", x = "Categories") +
        geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
        geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher), col=Full_summary$Color,size=8,position = position_fill(vjust = 0.45)) + 
        theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle = 45))
dev.off()

# Create Heatmap # 
Heatmap_summary <- Full_summary %>% dplyr::select(Node_label,Clones,oddsratio)

test <- tidyr::spread(Heatmap_summary, Clones, oddsratio)
test[is.na(test)] <- 1
M <- as.matrix(test %>% dplyr::select(-Node_label))
row.names(M) <- test$Node_label

inc <- 0.005
lowc <- c(seq(0,0.9-inc,by=inc))
colors <- c(seq(0.9,1.1,by=inc))
highc <- c(seq(1.1+inc,2,by=inc))
brs <- c(lowc,colors,highc)
my_palette <- c(colorRampPalette(colors=c("darkblue","cadetblue1"))(n=length(lowc)),colorRampPalette(colors = c("cadetblue1", "white", "yellow"))(n = length(colors)), colorRampPalette(colors = c("yellow", "red"))(n = length(highc)))

col_fun = colorRamp2(c(brs), c(my_palette))

ht1 <- Heatmap(M, name = "oddsRatio",rect_gp = gpar(col = "white", lwd = 0.3), cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
               column_labels = c("M. pudica","NF","Organogenesis","Release","Persist","Nfix",
                                 "M. truncatula","NF","FI","FIId","FIIp","IZ","ZIII","Lateral_root","Myc",
                                 "L. japonicus","NF","Myc"),
               cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 45, column_names_side = "top",
               row_names_side = "left", col = col_fun , width = unit(15, "cm"), height = unit(8, "cm"),row_gap = unit(25, "mm"))
ht1

pdffile = paste0("./Results_paper/Heatmap_oddsratio_Multispecies_release_multi_Nodes_Zones_Traits_Up.pdf")
pdf(pdffile,height=10,width=12)
par(mar=c(1,1,1,1))
draw(ht1)
dev.off()

Heatmap_summary <- Full_summary %>% dplyr::select(Node_label,Clones,oddsratio,FishTest)
Heatmap_summary$oddsratio[Heatmap_summary$FishTest>0.05] <- 1

test <- tidyr::spread(Heatmap_summary %>% dplyr::select(-FishTest), Clones, oddsratio)
test[is.na(test)] <- 1
M <- as.matrix(test %>% dplyr::select(-Node_label))
row.names(M) <- test$Node_label

inc <- 0.005
lowc <- c(seq(0,0.9-inc,by=inc))
colors <- c(seq(0.9,1.1,by=inc))
highc <- c(seq(1.1+inc,2,by=inc))
brs <- c(lowc,colors,highc)
my_palette <- c(colorRampPalette(colors=c("darkblue","cadetblue1"))(n=length(lowc)),colorRampPalette(colors = c("cadetblue1", "white", "yellow"))(n = length(colors)), colorRampPalette(colors = c("yellow", "red"))(n = length(highc)))

col_fun = colorRamp2(c(brs), c(my_palette))

ht1 <- Heatmap(M, name = "oddsRatio",rect_gp = gpar(col = "white", lwd = 0.3), cluster_rows = FALSE,
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
               column_labels = c("M. pudica","NF","Organogenesis","Release","Persist","Nfix",
                                 "M. truncatula","NF","FI","FIId","FIIp","IZ","ZIII","Lateral_root","Myc",
                                 "L. japonicus","NF","Myc"),
               cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 45, column_names_side = "top",
               row_names_side = "left", col = col_fun , width = unit(15, "cm"), height = unit(8, "cm"),row_gap = unit(25, "mm"))
ht1

pdffile = paste0("./Results_paper/Heatmap_oddsratiosignificant_Multispecies_release_multi_Nodes_Zones_Traits_Up.pdf")
pdf(pdffile,height=10,width=12)
par(mar=c(1,1,1,1))
draw(ht1)
dev.off()


#### Release and FIId specificity ####
treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")
MimpudNODES <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
MedtruNODES <- c(names(ancestors(tree, "Medtru")),"Medtru")

Mimpud_traits <- dbGetQuery(outcon, paste0("select * from ","Mimpud_UpDown"))
Mimpud_traits <- Mimpud_traits %>% dplyr::select(Gene_id,HOG,Mimpud_Release_Up, Mimpud_WSR_Nod_Up,Node_label_Up) %>%
  dplyr::filter(Node_label_Up=="Mimpud_specific") %>% dplyr::filter(Mimpud_WSR_Nod_Up=="Up")

# Test signalP
Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Mimpud_prot <- dplyr::select(Mimpud_prot,-Acronym)
Mimpud_prot$Gene_id <- row.names(Mimpud_prot)
Mimpud_prot <- Mimpud_prot %>% dplyr::filter(Gene_id %in% unique(Mimpud_traits$Gene_id))
Mimpud_prot$SignalP[Mimpud_prot$SignalP=="signal_peptide"] <- 1
Mimpud_prot$SignalP[Mimpud_prot$SignalP==""] <- 0; Mimpud_prot$SignalP <- as.numeric(Mimpud_prot$SignalP)

Mimpud_prot <- distinct(dplyr::left_join(Mimpud_traits,Mimpud_prot,by="Gene_id"))
Mimpud_prot <- Mimpud_prot %>% tidyr::unite("PASTE",all_of(c("Mimpud_Release_Up","Node_label_Up")), sep="_", remove=FALSE)
Mimpud_prot$PASTE[Mimpud_prot$PASTE=="_Mimpud_specific"] <- "noRelease"
Mimpud_prot$PASTE[Mimpud_prot$PASTE=="Up_Mimpud_specific"] <- "Release"

Mimpud_prot <- Mimpud_prot %>% group_by(PASTE,SignalP) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))

test <- Mimpud_prot %>% tidyr::spread(PASTE, Overlap) %>% dplyr::select(-SignalP)
M <- as.matrix(test)

fishertest <- fisher.test(M)
fishertest
fishertable <- data.frame(Test="Fisher's Exact Test for Count Data",pval=fishertest$p.value,oddsRatio=fishertest$estimate)

dbWriteTable(outcon, paste0("Mimpud_Release_signalP_stats"), fishertable, overwrite=TRUE)

#### Medtru FIId ####
Medtru_traits <- dbGetQuery(outcon, paste0("select * from ","Medtru_UpDown"))
Medtru_traits <- Medtru_traits %>% dplyr::select(Gene_id,HOG,Medtru_FIId_Up, Medtru_WSR_Nod_Up,Node_label_Up) %>%
  dplyr::filter(Node_label_Up=="Medtru_specific") %>% dplyr::filter(Medtru_WSR_Nod_Up=="Up")

# Test signalP
Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Medtru_prot <- dplyr::select(Medtru_prot,-Acronym)
Medtru_prot$Gene_id <- row.names(Medtru_prot)
Medtru_prot <- Medtru_prot %>% dplyr::filter(Gene_id %in% unique(Medtru_traits$Gene_id))
Medtru_prot$SignalP[Medtru_prot$SignalP=="signal_peptide"] <- 1
Medtru_prot$SignalP[Medtru_prot$SignalP==""] <- 0; Medtru_prot$SignalP <- as.numeric(Medtru_prot$SignalP)

Medtru_prot <- distinct(dplyr::left_join(Medtru_traits,Medtru_prot,by="Gene_id"))
Medtru_prot <- Medtru_prot %>% tidyr::unite("PASTE",all_of(c("Medtru_FIId_Up","Node_label_Up")), sep="_", remove=FALSE)
Medtru_prot$PASTE[Medtru_prot$PASTE!="Up_Medtru_specific"] <- "noRelease"
Medtru_prot$PASTE[Medtru_prot$PASTE=="Up_Medtru_specific"] <- "Release"

Medtru_prot <- Medtru_prot %>% group_by(PASTE,SignalP) %>% dplyr::summarise(Overlap=length(unique(Gene_id)))

test <- Medtru_prot %>% tidyr::spread(PASTE, Overlap) %>% dplyr::select(-SignalP)
M <- as.matrix(test)

fishertest <- fisher.test(M)
fishertest
fishertable <- data.frame(Test="Fisher's Exact Test for Count Data",pval=fishertest$p.value,oddsRatio=fishertest$estimate)

dbWriteTable(outcon, paste0("Medtru_FIId_signalP_stats"), fishertable, overwrite=TRUE)

#### PLS-DA on Mimosa signalP proteins ####
Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Mimpud_prot <- dplyr::select(Mimpud_prot,-Acronym)
Mimpud_prot$Gene_id <- row.names(Mimpud_prot)
Mimpud_prot <- Mimpud_prot %>% dplyr::filter(Gene_id %in% unique(Mimpud_traits$Gene_id)) %>% dplyr::filter(SignalP=="signal_peptide") %>% dplyr::select(-Gene_id,-SignalP) 

X<-data.matrix(t(Mimpud_prot), rownames.force = NA)
colors<-rep(c("black","blue1"))
colors2<-rep(c("black","blue1"), c(10,28))
Y = data.frame(Gene_id=row.names(Mimpud_prot))
TRAIT_df <- Mimpud_traits %>% dplyr::filter(Mimpud_Release_Up=="Up") %>% dplyr::select(Gene_id,Mimpud_Release_Up,Node_label_Up)
TRAIT_df <- TRAIT_df %>% dplyr::filter(Node_label_Up=="Mimpud_specific") %>% dplyr::select(Gene_id,Mimpud_Release_Up)
Y <- dplyr::left_join(Y,TRAIT_df,by="Gene_id")
Y[is.na(Y)] <- "noRelease"; Y$Mimpud_Release_Up[Y$Mimpud_Release_Up=="Up"] <- "Release"
Y = as.factor(Y$Mimpud_Release_Up) 

# PLS DA
res_plsda <- plsda(t(X), Y=Y, ncomp=4, near.zero.var=FALSE, scale=TRUE) 
# plotIndiv(res_plsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
plotVar(res_plsda, comp=c(1,2),cex=3)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Release species specific") # graphe de contribution des IPR à la séparation des LFA et nLFA 

pdf(file = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Loadings_Mimpud_Release_Up_aminoacids_on_signalP.pdf"), width=12, height=8)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Release species specific")
dev.off()

#### PLS-DA on Mimosa signalP proteins ####
Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Medtru_prot <- dplyr::select(Medtru_prot,-Acronym)
Medtru_prot$Gene_id <- row.names(Medtru_prot)
Medtru_prot <- Medtru_prot %>% dplyr::filter(Gene_id %in% unique(Medtru_traits$Gene_id)) %>% dplyr::filter(SignalP=="signal_peptide") %>% dplyr::select(-Gene_id,-SignalP) 

X<-data.matrix(t(Medtru_prot), rownames.force = NA)
colors<-rep(c("black","blue1"))
colors2<-rep(c("black","blue1"), c(10,28))
Y = data.frame(Gene_id=row.names(Medtru_prot))
TRAIT_df <- Medtru_traits %>% dplyr::filter(Medtru_FIId_Up=="Up") %>% dplyr::select(Gene_id,Medtru_FIId_Up,Node_label_Up)
TRAIT_df <- TRAIT_df %>% dplyr::filter(Node_label_Up=="Medtru_specific") %>% dplyr::select(Gene_id,Medtru_FIId_Up)
Y <- dplyr::left_join(Y,TRAIT_df,by="Gene_id")
Y[is.na(Y)] <- "noRelease"; Y$Medtru_FIId_Up[Y$Medtru_FIId_Up=="Up"] <- "Release"
Y = as.factor(Y$Medtru_FIId_Up) 

# PLS DA
res_plsda <- plsda(t(X), Y=Y, ncomp=4, near.zero.var=FALSE, scale=TRUE) 
# plotIndiv(res_plsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
plotVar(res_plsda, comp=c(1,2),cex=3)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Medtru FIId species specific") # graphe de contribution des IPR à la séparation des LFA et nLFA 

pdf(file = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Loadings_Medtru_FIId_Up_aminoacids_on_signalP.pdf"), width=12, height=8)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Medtru FIId species specific")
dev.off()

##### Mimpud Release Stats on protein characteristics ####
Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Mimpud_prot <- dplyr::select(Mimpud_prot,-Acronym)
Mimpud_prot$Gene_id <- row.names(Mimpud_prot)
Mimpud_prot <- Mimpud_prot %>% dplyr::filter(Gene_id %in% unique(Mimpud_traits$Gene_id))
Mimpud_prot$SignalP[Mimpud_prot$SignalP=="signal_peptide"] <- 1
Mimpud_prot$SignalP[Mimpud_prot$SignalP==""] <- 0; Mimpud_prot$SignalP <- as.numeric(Mimpud_prot$SignalP)

Mimpud_prot <- distinct(dplyr::left_join(Mimpud_traits,Mimpud_prot,by="Gene_id"))
Mimpud_prot <- Mimpud_prot %>% tidyr::unite("PASTE",all_of(c("Mimpud_Release_Up","Node_label_Up")), sep="_", remove=FALSE)
Mimpud_prot$PASTE[Mimpud_prot$PASTE!="Up_Mimpud_specific"] <- "noRelease"
Mimpud_prot$PASTE[Mimpud_prot$PASTE=="Up_Mimpud_specific"] <- "Release"
Mimpud_prot <- Mimpud_prot %>% dplyr::filter(SignalP==1)

categories <- names(Mimpud_prot %>% dplyr::select(Pep_length,ends_with("_prop")))
Stats_summary=data.frame(Category=categories)
for(i in categories){
  sub_medtru <- as.vector(as.data.frame(Mimpud_prot %>% dplyr::filter(PASTE=="noRelease") %>% dplyr::select(paste0(i))))
  sub_trait <- as.vector(as.data.frame(Mimpud_prot %>% dplyr::filter(PASTE=="Release") %>% dplyr::select(paste0(i))))
  t_test <- t.test(as.numeric(sub_trait[,1]),as.numeric(sub_medtru[,1]))
  Stats_summary$ttest_p.val[Stats_summary$Category==i] <- t_test$p.value
  Stats_summary$Mean_Release[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))
  Stats_summary$Mean_Mimpud[Stats_summary$Category==i] <- mean(as.numeric(sub_medtru[,1]))
  Stats_summary$Ratio_means[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))/mean(as.numeric(sub_medtru[,1]))
}

Stats_summary <- arrange(Stats_summary, Ratio_means)
Stats_summary$ttest_FDR <- p.adjust(Stats_summary$ttest_p.val, method="fdr")

dbWriteTable(outcon, paste0("Mimpud_Release_signalP_ProtStats"), Stats_summary, overwrite=TRUE)

##### Mimpud scatterplot ####
##### HighLight Proline rich short proteins ####
Biplot_data_signalP <- Mimpud_prot %>% dplyr::select(Pep_length,P_prop,H_prop,PASTE)

# Signal P Proline proportion and size
scatterPlot <- ggplot(Biplot_data_signalP,aes(P_prop, Pep_length, color=PASTE)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 10)+
  geom_hline(yintercept = 200)+
  theme_classic()+theme(legend.position=c(0.9,1), legend.justification=c(1,1),axis.text.y = element_text(angle=90,hjust=0.5))

xdensity <- ggplot(Biplot_data_signalP, aes(P_prop, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_vline(xintercept = 10) +
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5))

ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_hline(yintercept = 200) +
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5)) 

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdffile = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Mimpud_Release_signalP_Proline_length.pdf")
pdf(pdffile,height=7,width=7)
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

##### HighLight Cysteine rich short proteins ####
Biplot_data_signalP <- Mimpud_prot %>% dplyr::select(Pep_length,C_prop,H_prop,PASTE)

# Signal C Cysteine proportion and size
scatterPlot <- ggplot(Biplot_data_signalP,aes(C_prop, Pep_length, color=PASTE)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_hline(yintercept = 100)+
  theme_classic()+theme(legend.position=c(0.9,1), legend.justification=c(1,1),axis.text.y = element_text(angle=90,hjust=0.5))

xdensity <- ggplot(Biplot_data_signalP, aes(C_prop, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_vline(xintercept = 6) +
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5))

ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_hline(yintercept = 100) +
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5)) 

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdffile = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Mimpud_Release_signalP_Cysteine_length.pdf")
pdf(pdffile,height=7,width=7)
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()



##### Medtru Release Stats on protein characteristics ####
Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
Medtru_prot <- dplyr::select(Medtru_prot,-Acronym)
Medtru_prot$Gene_id <- row.names(Medtru_prot)
Medtru_prot <- Medtru_prot %>% dplyr::filter(Gene_id %in% unique(Medtru_traits$Gene_id))
Medtru_prot$SignalP[Medtru_prot$SignalP=="signal_peptide"] <- 1
Medtru_prot$SignalP[Medtru_prot$SignalP==""] <- 0; Medtru_prot$SignalP <- as.numeric(Medtru_prot$SignalP)

Medtru_prot <- distinct(dplyr::left_join(Medtru_traits,Medtru_prot,by="Gene_id"))
Medtru_prot <- Medtru_prot %>% tidyr::unite("PASTE",all_of(c("Medtru_FIId_Up","Node_label_Up")), sep="_", remove=FALSE)
Medtru_prot$PASTE[Medtru_prot$PASTE!="Up_Medtru_specific"] <- "noFIId"
Medtru_prot$PASTE[Medtru_prot$PASTE=="Up_Medtru_specific"] <- "FIId"
Medtru_prot <- Medtru_prot %>% dplyr::filter(SignalP==1)

categories <- names(Medtru_prot %>% dplyr::select(Pep_length,ends_with("_prop")))
Stats_summary=data.frame(Category=categories)
for(i in categories){
  sub_medtru <- as.vector(as.data.frame(Medtru_prot %>% dplyr::filter(PASTE=="noFIId") %>% dplyr::select(paste0(i))))
  sub_trait <- as.vector(as.data.frame(Medtru_prot %>% dplyr::filter(PASTE=="FIId") %>% dplyr::select(paste0(i))))
  t_test <- t.test(as.numeric(sub_trait[,1]),as.numeric(sub_medtru[,1]))
  Stats_summary$ttest_p.val[Stats_summary$Category==i] <- t_test$p.value
  Stats_summary$Mean_FIId[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))
  Stats_summary$Mean_Medtru[Stats_summary$Category==i] <- mean(as.numeric(sub_medtru[,1]))
  Stats_summary$Ratio_means[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))/mean(as.numeric(sub_medtru[,1]))
}

Stats_summary <- arrange(Stats_summary, Ratio_means)
Stats_summary$ttest_FDR <- p.adjust(Stats_summary$ttest_p.val, method="fdr")

dbWriteTable(outcon, paste0("Medtru_FIId_signalP_ProtStats"), Stats_summary, overwrite=TRUE)

#### Medtru scatterplot ####
##### HighLight Proline rich short proteins ####
Biplot_data_signalP <- Medtru_prot %>% dplyr::select(Pep_length,P_prop,H_prop,PASTE)

# Signal P Proline proportion and size
scatterPlot <- ggplot(Biplot_data_signalP,aes(P_prop, Pep_length, color=PASTE)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#E69F00','#999999')) + 
  geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 10)+
  geom_hline(yintercept = 200)+
  theme_classic()+theme(legend.position=c(0.9,1), legend.justification=c(1,1),axis.text.y = element_text(angle=90,hjust=0.5))

xdensity <- ggplot(Biplot_data_signalP, aes(P_prop, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_vline(xintercept = 10) +
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5))

ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_hline(yintercept = 200) +
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5)) 

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdffile = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Medtru_FIId_signalP_Proline_length.pdf")
pdf(pdffile,height=7,width=7)
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

##### HighLight Cysteine rich short proteins ####
Biplot_data_signalP <- Medtru_prot %>% dplyr::select(Pep_length,C_prop,H_prop,PASTE)

# Signal C Cysteine proportion and size
scatterPlot <- ggplot(Biplot_data_signalP,aes(C_prop, Pep_length, color=PASTE)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#E69F00','#999999')) + 
  geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_hline(yintercept = 100)+
  theme_classic()+theme(legend.position=c(0.9,1), legend.justification=c(1,1),axis.text.y = element_text(angle=90,hjust=0.5))

xdensity <- ggplot(Biplot_data_signalP, aes(C_prop, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_vline(xintercept = 6) +
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5))

ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=PASTE)) + 
  geom_density(alpha=.5) + 
  geom_hline(yintercept = 100) +
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme_classic() + theme(legend.position = "none",axis.text.y = element_text(angle=90,hjust=0.5)) 

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdffile = paste0("C:/Users/cyril.libourel/OneDrive/mimosaProject/figures/Medtru_FIId_signalP_Cysteine_length.pdf")
pdf(pdffile,height=7,width=7)
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

