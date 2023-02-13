######################
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","DBI","grid","mixOmics","DECIPHER","Biostrings","gridExtra",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ComplexHeatmap","circlize","dbplyr","wesanderson"), library, character.only = TRUE)



setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")
sqlite_db_out <- "C:/Users/cyril.libourel/OneDrive/mimosaProject/work/Sql_db.sqlite"
outcon <- dbConnect(RSQLite::SQLite(), dbname=sqlite_db_out)

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

dbWriteTable(outcon, paste0("N0_corrected_Gene_id_HOGs_correspondance"), All_df, overwrite=TRUE)
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")
for(sp in species){
  cds <- Biostrings::readDNAStringSet(paste0("./cds_DB/",sp,"_longest_isoform.cds"))
  Gene_id_df <- data.frame(Gene_id=gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds)))
  HOG <- full_join(HOG,Gene_id_df, by="Gene_id")
}

HOG[is.na(HOG)] <- ""
dbWriteTable(outcon, paste0("N0_corrected_Gene_id_HOGs_correspondance"), HOG, overwrite=TRUE)
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

# Up
All_HOGs <- data.frame()
tables_up <- grep("Up",tables,value=TRUE)

for(i in tables_up){
  sp <- gsub("./rnaseq/","",unlist(strsplit(i,split="_"))[1])
  DEGs <- fread(paste0(i,".txt"),h=FALSE)
  DEGs <- distinct(DEGs)
  DEGs$V2 <- unlist(strsplit(i,split="_"))[5]
  names(DEGs) <- c("Gene_id",paste0("Nod_WSR_",unlist(strsplit(i,split="_"))[5]))
  sp_HOGs <- HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
  sp_HOGs <- dplyr::left_join(sp_HOGs,DEGs,by="Gene_id")
  All_HOGs <- rbind(All_HOGs,sp_HOGs)
}

All_HOGs[is.na(All_HOGs)] <- ""
dbWriteTable(outcon, paste0("N0_DE_Up"), All_HOGs, overwrite=TRUE)

# Down
All_HOGs <- data.frame()
tables_Down <- grep("Down",tables,value=TRUE)

for(i in tables_Down){
  sp <- gsub("./rnaseq/","",unlist(strsplit(i,split="_"))[1])
  DEGs <- fread(paste0(i,".txt"),h=FALSE)
  DEGs <- distinct(DEGs)
  DEGs$V2 <- unlist(strsplit(i,split="_"))[5]
  names(DEGs) <- c("Gene_id",paste0("Nod_WSR_",unlist(strsplit(i,split="_"))[5]))
  sp_HOGs <- HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
  sp_HOGs <- dplyr::left_join(sp_HOGs,DEGs,by="Gene_id")
  All_HOGs <- rbind(All_HOGs,sp_HOGs)
}

All_HOGs[is.na(All_HOGs)] <- ""
dbWriteTable(outcon, paste0("N0_DE_Down"), All_HOGs, overwrite=TRUE)

##### DEOG based on HOG #####
# Up
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_DE_Up"))

All_HOG <- data.frame(HOG=unique(HOG$HOG))
species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

for(sp in species){
  sub_HOG <- HOG %>% dplyr::filter(Nod_WSR_Up=="Up") %>% dplyr::filter(grepl(sp,Gene_id)) %>% group_by(HOG) %>% dplyr::count(HOG)
  
  sub_HOG$HOG[sub_HOG$HOG==""] <-"Species_specific"
  names(sub_HOG) <- c("HOG",paste0(sp,"_WSR_Up"))
  
  All_HOG <- dplyr::left_join(All_HOG,sub_HOG,by="HOG")

}

All_HOG[is.na(All_HOG)] <- 0;All_HOG$HOG <- unique(HOG$HOG)
dbWriteTable(outcon, paste0("N0_DEOG_Up"), All_HOG, overwrite=TRUE)

All_HOG[is.na(All_HOG)] <- 0;All_HOG[All_HOG>0] <- 1;All_HOG$HOG <- unique(HOG$HOG)

dbWriteTable(outcon, paste0("N0_DEOG_Up_binary"), All_HOG, overwrite=TRUE)


# Down
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_DE_Down"))

All_HOG <- data.frame(HOG=unique(HOG$HOG))
species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")

for(sp in species){
  sub_HOG <- HOG %>% dplyr::filter(Nod_WSR_Down=="Down") %>% dplyr::filter(grepl(sp,Gene_id)) %>% group_by(HOG) %>% dplyr::count(HOG)
  
  sub_HOG$HOG[sub_HOG$HOG==""] <-"Species_specific"
  names(sub_HOG) <- c("HOG",paste0(sp,"_WSR_Down"))
  
  All_HOG <- dplyr::left_join(All_HOG,sub_HOG,by="HOG")
  
}

All_HOG[is.na(All_HOG)] <- 0;All_HOG$HOG <- unique(HOG$HOG)
dbWriteTable(outcon, paste0("N0_DEOG_Down"), All_HOG, overwrite=TRUE)
All_HOG[is.na(All_HOG)] <- 0;All_HOG[All_HOG>0] <- 1;All_HOG$HOG <- unique(HOG$HOG)

dbWriteTable(outcon, paste0("N0_DEOG_Down_binary"), All_HOG, overwrite=TRUE)

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

# HOG <- fread(paste0("./results/N0_corrected_Gene_id_HOGs_correspondance.txt"),h=T)
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

HOG <- HOG %>% filter(str_detect(Gene_id, "^Medtru_"))
Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
HOG <- as.data.frame(left_join(HOG,Medtru_annot))
HOG[is.na(HOG)] <- ""

Concat_HOG <- HOG %>% group_by(HOG) %>% transmute(Acronym = paste0(sort(unique(unlist(strsplit(Acronym,",")))),collapse=","),
                                                  Category = paste0(sort(unique(Category)),collapse="|"),
                                                  Description = paste0(unique(Description),collapse=","),
                                                  IPR = paste0(sort(unique(unlist(strsplit(IPR,",")))),collapse=","))
Concat_HOG <- distinct(Concat_HOG)

All_HOGs_Up <- dbGetQuery(outcon, paste0("select * from ","N0_DEOG_Up"))
All_HOGs_Down <- dbGetQuery(outcon, paste0("select * from ","N0_DEOG_Down"))
All_HOGs <- left_join(All_HOGs_Up,All_HOGs_Down,by="HOG")
# All_HOGs <- fread(paste0("./results/N0_corrected_nodMimosa_project.tsv"))
All_HOGs <- left_join(All_HOGs,Concat_HOG)
All_HOGs <- left_join(All_HOGs,HOG_count)
All_HOGs <- left_join(All_HOGs,HOG_species_Nod)
All_HOGs <- left_join(All_HOGs,HOG_species_noNod)
All_HOGs[is.na(All_HOGs)] <- ""
All_HOGs$Percentage_of_nonNodspecies_in_HOG <- round(as.numeric(All_HOGs$Percentage_of_nonNodspecies_in_HOG),2)
All_HOGs$Percentage_of_Nodspecies_in_HOG <- round(as.numeric(All_HOGs$Percentage_of_Nodspecies_in_HOG),2)
All_HOGs$Mean_seq_number_NodspecieswithDEGs <- round(as.numeric(All_HOGs$Mean_seq_number_NodspecieswithDEGs),2)
All_HOGs$Mean_seq_number_HOG <- round(as.numeric(All_HOGs$Mean_seq_number_HOG),2)

# fwrite(All_HOGs,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot.tsv"),sep="\t",quote=FALSE)

dbWriteTable(outcon, paste0("N0_corrected"), All_HOGs, overwrite=TRUE)

All_HOGs <- dbGetQuery(outcon, paste0("select * from ","N0_corrected"))
All_HOGs[is.na(All_HOGs)] <- ""
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

Resume_HOGs <- dplyr::select(All_HOGs,HOG,Acronym, Category,Description,IPR,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                             Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)

HOG_NAMES <- setdiff(names(All_HOGs),names(Resume_HOGs))

for(i in HOG_NAMES){
  species_HOG <- unique(HOG %>% filter(str_detect(Gene_id, paste0("^",str_split(i,pattern="_")[[1]][1]))) %>% dplyr::select(HOG))$HOG
  
  subset_HOG_present <- All_HOGs %>% dplyr::select(HOG, paste0(i)) %>% filter(HOG %in% species_HOG)
  names(subset_HOG_present) <- c("HOG","Genes")
  subset_HOG_present$Genes[subset_HOG_present$Genes==0] <- "ns"
  
  subset_HOG_absent <- All_HOGs %>% dplyr::select(HOG, paste0(i)) %>% filter(!HOG %in% species_HOG)
  names(subset_HOG_absent) <- c("HOG","Genes")
  subset_HOG <- rbind(subset_HOG_present,subset_HOG_absent)
  names(subset_HOG) <- c("HOG",i)
  Resume_HOGs <- dplyr::left_join(Resume_HOGs,subset_HOG,by="HOG")
}

Resume_HOGs$Description <- gsub("^,","",Resume_HOGs$Description)
Resume_HOGs$Description <- gsub(",$","",Resume_HOGs$Description)
Resume_HOGs$Category <- gsub("^\\|","",Resume_HOGs$Category)
Resume_HOGs$Category <- gsub("\\|$","",Resume_HOGs$Category)

dbWriteTable(outcon, paste0("N0_corrected"), Resume_HOGs, overwrite=TRUE)

#### Ranking of HOG ####
Resume_HOGs <- dbGetQuery(outcon, paste0("select * from ","N0_corrected"))
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

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
  temp_DEGs <- as.data.frame(dplyr::left_join(DEGs,HOG, by="Gene_id"))
  names(temp_DEGs) <- c(gsub("./rnaseq/","",i),"HOG")
  temp_DEGs[is.na(temp_DEGs)] <- ""
  sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
  count_DEGs <- temp_DEGs %>% group_by(HOG) %>% dplyr::count(HOG);names(count_DEGs) <- c("HOG","DEGs")
  count_genesHOG <- HOG %>% dplyr::filter(grepl(sp,Gene_id)) %>% group_by(HOG) %>% dplyr::count(HOG) ;names(count_genesHOG) <- c("HOG","seqcount")
  
  HOG_sp <- Resume_HOGs %>% dplyr::select(HOG)
  HOG_sp <- dplyr::left_join(HOG_sp,count_genesHOG,by="HOG");HOG_sp <- dplyr::left_join(HOG_sp,count_DEGs,by="HOG")
  HOG_sp[is.na(HOG_sp)] <- 0
  HOG_sp$Ranking <- round((HOG_sp$DEGs/HOG_sp$seqcount)*100,2)
  HOG_sp <- dplyr::select(HOG_sp,HOG,Ranking)
  names(HOG_sp) <- c("HOG",paste0(sp,"_ranking"))
  
  All_ranking <- dplyr::left_join(All_ranking,HOG_sp,by="HOG")
}
  All_ranking$Ranking <- round(rowMeans(dplyr::select(All_ranking,ends_with("_ranking")), na.rm=TRUE),2)
  All_ranking <- All_ranking %>% dplyr::select(HOG,Ranking);names(All_ranking) <- c("HOG",paste0("Ranking_",updown))
  Resume_HOGs <- dplyr::left_join(Resume_HOGs,All_ranking, by="HOG")

}

All_HOGs <- dplyr::select(Resume_HOGs,HOG,Acronym, Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,Ranking_Down,
                          Ranking_Up,Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
All_HOGs2 <- dplyr::select(Resume_HOGs,-Acronym, -Category,-Description,-IPR,-Mean_seq_number_HOG,-Mean_seq_number_NodspecieswithDEGs,
                           -Ranking_Down,-Ranking_Up,-Percentage_of_Nodspecies_in_HOG,-Percentage_of_nonNodspecies_in_HOG)
All_HOGs <- left_join(All_HOGs,All_HOGs2,by="HOG")
dbWriteTable(outcon, paste0("N0_corrected"), All_HOGs, overwrite=TRUE)

