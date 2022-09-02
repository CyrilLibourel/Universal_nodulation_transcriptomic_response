######################
lapply(c("data.table","dplyr","stringr"), require, character.only = TRUE)


setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")

##### Write files with HOG and gene_id correspondance HOG #####
HOG <- as.data.frame(fread(paste0("../Results_Apr14/Phylogenetic_Hierarchical_Orthogroups/N0_corrected.tsv"),h=T))
HOG_NAMES <- names(HOG)[c(4:ncol(HOG))]

All_df <- data.table()

for(i in HOG_NAMES){
  test_HOG <- dplyr:::select(HOG,HOG,paste(i))
  test_list <- strsplit(test_HOG[,2],", ")
  names(test_list) <- HOG$HOG
  df <- reshape2::melt(test_list)
  names(df) <- c("Gene_id","HOG")
  df$Gene_id <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",df$Gene_id)

  All_df <- rbind(All_df,df)
  message(i," Done")
}

fwrite(All_df,paste0("./results/N0_corrected_Gene_id_HOGs_correspondance.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

All_HOGs <- data.table(HOG=unique(HOG$HOG))

for(i in tables){
  DEGs <- fread(paste0(i,".txt"),h=FALSE)
  names(DEGs) <- "Gene_id"
  temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
  names(temp) <- c(gsub("./rnaseq/","",i),"HOG")
  temp[is.na(temp)] <- ""
  temp[temp[,1] != "",1] <- unlist(strsplit(i,split="_"))[5]
  temp <- temp[!duplicated(temp[,c('HOG')]),]
  All_HOGs <- dplyr::left_join(All_HOGs,temp)
  All_HOGs[is.na(All_HOGs)] <- ""
}

fwrite(All_HOGs,paste0("./results/N0_corrected_nodMimosa_project.tsv"),sep="\t",quote=FALSE)

##### Add Medtr annotation and sequence number by HOG #####

HOG <- as.data.frame(fread(paste0("../Results_Apr14/Phylogenetic_Hierarchical_Orthogroups/N0_corrected.tsv"),h=T))
HOG_NAMES <- names(HOG)[c(4:ncol(HOG))]

HOG_count <- HOG[,c(4:length(HOG))]
HOG_count <- HOG_count %>% dplyr::mutate_each(funs(unlist(lapply(strsplit(.,","),length))))
rownames(HOG_count) <- HOG[,1]

HOG_species <- HOG_count
HOG_species[HOG_species>0] <- 1

HOG_species_Nod <- HOG_species %>% rowwise() %>% 
  transmute(Percentage_of_Nodspecies_in_HOG = sum(c_across(all_of(c("Aeseve","Alnglu","Arahyp","Casgla","Datglo","Distri","Drydru","Glymax","Hiprha","Lotjap","Lupalb","Medtru","Mimpud","Parand"))),na.rm=T)/14*100)
HOG_species_Nod$HOG <- HOG[,1]

HOG_species_noNod <- HOG_species %>% rowwise() %>% 
  transmute(Percentage_of_nonNodspecies_in_HOG = sum(c_across(all_of(c("Begfuc","Betpen","Carfan","Casaus","Cucmax","Fraves","Jugreg","Lagsic","Mornot","Nissch","Treori"))),na.rm=T)/11*100)
HOG_species_noNod$HOG <- HOG[,1]

is.na(HOG_count) <- HOG_count==0

HOG_count <- HOG_count %>% rowwise() %>% 
  transmute(Mean_seq_number_HOG = mean(c_across(all_of(HOG_NAMES)),na.rm=T),
            Mean_seq_number_NodspecieswithDEGs = mean(c_across(all_of(c("Aeseve","Arahyp","Datglo","Hiprha","Lotjap","Medtru","Mimpud","Parand","Lupalb"))),na.rm=T))
HOG_count$HOG <- HOG[,1]

HOG <- fread(paste0("./results/N0_corrected_Gene_id_HOGs_correspondance.txt"),h=T)
HOG <- HOG %>% filter(str_detect(Gene_id, "^Medtru_"))
Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
HOG <- as.data.frame(left_join(HOG,Medtru_annot))
HOG[is.na(HOG)] <- ""

Concat_HOG <- HOG %>% group_by(HOG) %>% transmute(Acronym = paste0(sort(unique(unlist(strsplit(Acronym,",")))),collapse=","),
                                                  Category = paste0(sort(unique(Category)),collapse="|"),
                                                  Description = paste0(unique(Description),collapse=","),
                                                  IPR = paste0(sort(unique(unlist(strsplit(IPR,",")))),collapse=","))
Concat_HOG <- distinct(Concat_HOG)

All_HOGs <- fread(paste0("./results/N0_corrected_nodMimosa_project.tsv"))
All_HOGs <- left_join(All_HOGs,Concat_HOG)
All_HOGs <- left_join(All_HOGs,HOG_count)
All_HOGs <- left_join(All_HOGs,HOG_species_Nod)
All_HOGs <- left_join(All_HOGs,HOG_species_noNod)
All_HOGs[is.na(All_HOGs)] <- ""

fwrite(All_HOGs,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot.tsv"),sep="\t",quote=FALSE)

All_HOGs <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot.tsv"))
All_HOGs[is.na(All_HOGs)] <- ""
HOG <- fread(paste0("./results/N0_corrected_Gene_id_HOGs_correspondance.txt"),h=T)
HOG_NAMES <- names(All_HOGs)[c(2:(ncol(All_HOGs)-6))]

Resume_HOGs <- dplyr::select(All_HOGs,HOG,Acronym, Category,Description,IPR,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                             Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)

for(i in HOG_NAMES){
  species_HOG <- unique(HOG %>% filter(str_detect(Gene_id, paste0("^",str_split(i,pattern="_")[[1]][1]))) %>% dplyr::select(HOG))$HOG
  
  subset_HOG_present <- All_HOGs %>% dplyr::select(HOG, paste0(i)) %>% filter(HOG %in% species_HOG)
  names(subset_HOG_present) <- c("HOG","Genes")
  subset_HOG_present$Genes[subset_HOG_present$Genes==""] <- "ns"
  
  subset_HOG_absent <- All_HOGs %>% dplyr::select(HOG, paste0(i)) %>% filter(!HOG %in% species_HOG)
  names(subset_HOG_absent) <- c("HOG","Genes")
  subset_HOG <- rbind(subset_HOG_present,subset_HOG_absent)
  names(subset_HOG) <- c("HOG",i)
  Resume_HOGs <- left_join(Resume_HOGs,subset_HOG)
}

Resume_HOGs$Description <- gsub("^,","",Resume_HOGs$Description)
Resume_HOGs$Description <- gsub(",$","",Resume_HOGs$Description)
Resume_HOGs$Category <- gsub("^\\|","",Resume_HOGs$Category)
Resume_HOGs$Category <- gsub("\\|$","",Resume_HOGs$Category)
fwrite(Resume_HOGs,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes.tsv"),sep="\t",quote=FALSE)

#### Ranking of HOG ####
Resume_HOGs <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes.tsv"))
HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)
All_HOGs <- data.table(HOG=unique(HOG$HOG))

UPDOWN <- c("Down","Up")

for(updown in UPDOWN){
  subtables <- grep(updown,tables, value=TRUE)
  All_ranking <- Resume_HOGs %>% dplyr::select(HOG)
for(i in subtables){
  DEGs <- fread(paste0(i,".txt"),h=FALSE)
  names(DEGs) <- "Gene_id"
  temp_DEGs <- as.data.frame(dplyr::left_join(DEGs,HOG))
  names(temp_DEGs) <- c(gsub("./rnaseq/","",i),"HOG")
  temp_DEGs[is.na(temp_DEGs)] <- ""
  sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
  count_DEGs <- temp_DEGs %>% group_by(HOG) %>% dplyr::count(HOG);names(count_DEGs) <- c("HOG","DEGs")
  count_genesHOG <- HOG %>% dplyr::filter(grepl(sp,Gene_id)) %>% group_by(HOG) %>% dplyr::count(HOG) ;names(count_genesHOG) <- c("HOG","seqcount")
  
  HOG_sp <- Resume_HOGs %>% dplyr::select(HOG)
  HOG_sp <- left_join(HOG_sp,count_genesHOG,by="HOG");HOG_sp <- left_join(HOG_sp,count_DEGs,by="HOG")
  HOG_sp[is.na(HOG_sp)] <- 0
  HOG_sp$Ranking <- round((HOG_sp$DEGs/HOG_sp$seqcount)*100,2)
  HOG_sp <- dplyr::select(HOG_sp,HOG,Ranking)
  names(HOG_sp) <- c("HOG",paste0(sp,"_ranking"))
  
  All_ranking <- left_join(All_ranking,HOG_sp,by="HOG")
}
  All_ranking$Ranking <- round(rowMeans(dplyr::select(All_ranking,ends_with("_ranking")), na.rm=TRUE),2)
  All_ranking <- All_ranking %>% dplyr::select(HOG,Ranking);names(All_ranking) <- c("HOG",paste0("Ranking_",updown))
  Resume_HOGs <- left_join(Resume_HOGs,All_ranking)

}

All_HOGs <- dplyr::select(Resume_HOGs,HOG,Acronym, Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,Ranking_Down,
                          Ranking_Up,Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
All_HOGs2 <- dplyr::select(Resume_HOGs,-Acronym, -Category,-Description,-IPR,-Mean_seq_number_HOG,-Mean_seq_number_NodspecieswithDEGs,
                           -Ranking_Down,-Ranking_Up,-Percentage_of_Nodspecies_in_HOG,-Percentage_of_nonNodspecies_in_HOG)
All_HOGs <- left_join(All_HOGs,All_HOGs2)
fwrite(All_HOGs,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking.tsv"),sep="\t",quote=FALSE)

