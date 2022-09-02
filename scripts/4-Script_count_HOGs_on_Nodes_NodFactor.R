# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","reshape2","ComplexHeatmap","ggVennDiagram",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ggpubr","wesanderson","DECIPHER","Biostrings","gridExtra","ggvenn"), require, character.only = TRUE)

source(file="//194.199.55.66/evo/commun/projects/nodMimosa/Ranking_HOG_nodMimosa/enrichment_analysis.R")
# path to input and output folder
slashtype = substr(unique(gsub("[^\\/]", "", getwd())),1,1)
PATH = paste0(getwd(),slashtype,"res",slashtype)

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")
out_path <- "//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3/results/"

# create output folder if necessary
if(dir.exists(out_path)==FALSE){dir.create(out_path,recursive = T)}
# load orthogroup file
HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

######## Up ########
UP_or_DOWN <- c("Up","Down")

for(UP_DOWN in UP_or_DOWN){
Medtru_Correspond <- unique(fread(input = paste0("//194.199.55.66/evo/commun/projects/nodMimosa/Ranking_HOG_nodMimosa/results/Medtru_Gene_id_correspondance_HOG_WT_UpDown_Nod.txt"),header=TRUE))
Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits",UP_DOWN,".tsv"))
# Recap_HOG <- Recap_HOG %>% dplyr::select(HOG,starts_with("Node_label_"),Acronym,Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
#                                          Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
# Medtru
Medtru_NF_Up <- unique(fread(input = paste0("../Nod_factor/Medtru_NF_logFC0_",UP_DOWN,".txt"),header=FALSE))
names(Medtru_NF_Up) <- "Gene_id"
Medtru_NF_Up$Gene_id <- unique(Medtru_NF_Up$Gene_id)
Medtru_NF_Up$NF <- UP_DOWN

Medtru_NF_Up <- left_join(Medtru_NF_Up,HOG,by="Gene_id")

Medtru_Up <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0("Medtru_.*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
names(Medtru_Up) <- "Gene_id"
Medtru_Up$Gene_id <- unique(Medtru_Up$Gene_id)
Medtru_Up$Medtru <- UP_DOWN
Medtru_NF_Up <- left_join(Medtru_NF_Up,Medtru_Up,by="Gene_id")

Medtru_NF_Up[is.na(Medtru_NF_Up)] <- ""
Medtru_NF_Up <- Medtru_NF_Up %>% dplyr::filter(NF==UP_DOWN, Medtru==UP_DOWN, HOG!="")

Medtru_NF_Up_Recap <- left_join(Medtru_NF_Up,Recap_HOG,by="HOG")

# Lotjap
Lotjap_NF_Up <- unique(fread(input = paste0("../Nod_factor/Lotjap_NF_logFC0_",UP_DOWN,".txt"),header=FALSE))
names(Lotjap_NF_Up) <- "Gene_id"
Lotjap_NF_Up$Gene_id <- unique(Lotjap_NF_Up$Gene_id)
Lotjap_NF_Up$NF <- UP_DOWN

Lotjap_NF_Up <- left_join(Lotjap_NF_Up,HOG,by="Gene_id")

Lotjap_Up <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0("Lotjap_.*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
names(Lotjap_Up) <- "Gene_id"
Lotjap_Up$Gene_id <- unique(Lotjap_Up$Gene_id)
Lotjap_Up$Lotjap <- UP_DOWN
Lotjap_NF_Up <- left_join(Lotjap_NF_Up,Lotjap_Up,by="Gene_id")

Lotjap_NF_Up[is.na(Lotjap_NF_Up)] <- ""
Lotjap_NF_Up <- Lotjap_NF_Up %>% dplyr::filter(NF==UP_DOWN, Lotjap==UP_DOWN, HOG!="")

Lotjap_NF_Up_Recap <- left_join(Lotjap_NF_Up,Recap_HOG,by="HOG")

# Mimpud
Mimpud_NF_Up <- unique(fread(input = paste0("./trait_analysis/NF.txt"),header=TRUE))
names(Mimpud_NF_Up) <- c("Gene_id","Mimpud")
Mimpud_NF_Up <- Mimpud_NF_Up %>% dplyr::filter(Mimpud==UP_DOWN) %>% dplyr::select(Gene_id)
Mimpud_NF_Up$Gene_id <- unique(Mimpud_NF_Up$Gene_id)
Mimpud_NF_Up$Mimpud <- UP_DOWN

Mimpud_NF_Up <- left_join(Mimpud_NF_Up,HOG,by="Gene_id")

Mimpud_NF_Up[is.na(Mimpud_NF_Up)] <- ""
Mimpud_NF_Up <- Mimpud_NF_Up %>% dplyr::filter(Mimpud==UP_DOWN, HOG!="")

Mimpud_NF_Up_Recap <- left_join(Mimpud_NF_Up,Recap_HOG,by="HOG")

NF_data <- dplyr::full_join(Mimpud_NF_Up %>% dplyr::select(HOG, Mimpud),Medtru_NF_Up %>% dplyr::select(HOG, Medtru), by="HOG")
NF_data <- dplyr::full_join(NF_data,Lotjap_NF_Up %>% dplyr::select(HOG, Lotjap), by="HOG")
NF_data[is.na(NF_data)] <- ""

HOG_all_df <- left_join(Recap_HOG,NF_data,by="HOG")
HOG_all_df[is.na(HOG_all_df)] <- ""
HOG_all_df <- as_tibble(HOG_all_df)
HOG_all_df <- dplyr::rename(HOG_all_df, Mimpud_NodFactor = Mimpud)
HOG_all_df <- dplyr::rename(HOG_all_df, Medtru_NodFactor = Medtru)
HOG_all_df <- dplyr::rename(HOG_all_df, Lotjap_NodFactor = Lotjap)
HOG_all_df <- HOG_all_df %>% distinct()

fwrite(HOG_all_df,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF",UP_DOWN,".tsv"), quote=FALSE,sep="\t")


x <- list(
  M.truncatula = NF_data$HOG[NF_data$Medtru==UP_DOWN], 
  L.japonicus = NF_data$HOG[NF_data$Lotjap==UP_DOWN],
  M.pudica = NF_data$HOG[NF_data$Mimpud==UP_DOWN]
)

pdffile = paste0("./results/Jvenn_NF_",UP_DOWN,".pdf")
pdf(pdffile,height=8,width=8)
ggvenn(x,fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4)
dev.off()


Recap_HOG <- Recap_HOG %>% dplyr::select(HOG,starts_with("Node_label_"),Acronym,Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                                         Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
NF_data_HOG <- dplyr::left_join(NF_data,Recap_HOG, by="HOG")

Medtru_NF_Up_Recap <- Medtru_NF_Up_Recap
Lotjap_NF_Up_Recap <- Lotjap_NF_Up_Recap
Mimpud_NF_Up_Recap <- Mimpud_NF_Up_Recap

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

Summary_count <- data.frame()
i="Medtru"
  DEGs <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0(i,".*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
  names(DEGs) <- "Gene_id"
  temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
  # names(temp) <- c(gsub("./rnaseq/","",i),"HOG")
  names(temp) <- c(i,"HOG")
  temp[is.na(temp)] <- ""
  # temp <- temp[!duplicated(temp[,c('HOG')]),]
  cor_Genes_HOG <- dplyr::left_join(temp,Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN)))
  
  names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
  NODES <- c(names(ancestors(tree, i)),i)
  cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
  cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
  cor_Genes_HOG$species <- i
  cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
  cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
  Summary_count <- rbind(Summary_count,cor_Genes_HOG)
  
  cor_Genes_HOG <- Medtru_NF_Up_Recap%>% dplyr::select(Gene_id,HOG,paste0("Node_label_",UP_DOWN))
  
  names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
  cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
  cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
  cor_Genes_HOG$species <- "NF"
  cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
  cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
  Summary_count <- rbind(Summary_count,cor_Genes_HOG)

  
  Full_summary <- Summary_count
  
  Clones_sum <- data.frame(species=unique(Full_summary$species),Sum=0)
  for(sp in unique(Full_summary$species)){
    Clones_sum$Sum[Clones_sum$species==sp] <- sum(Full_summary$Overlap[Full_summary$species==sp])
    
  }

  Categories <- unique(Full_summary$Node_label)
  for(Node in Categories){
      CloneNode <- Full_summary$Overlap[Full_summary$species=="NF" & Full_summary$Node_label==Node]
      CtaiNode <- Full_summary$Overlap[Full_summary$species==i & Full_summary$Node_label==Node]
      M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$species=="NF"]-CloneNode,
                                                    Clones_sum$Sum[Clones_sum$species==i]-CtaiNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- test$p.value
      Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- fishertest$estimate
      Full_summary$Total[Full_summary$Node_label=="Mimpud_HOG" & Full_summary$species=="NF"] <- Clones_sum$Sum[Clones_sum$species=="NF"]
  }
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than WT"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than WT"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than WT"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than WT"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than WT"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than WT"
  
  fwrite(Full_summary,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Medtru_NodFactor_multi_Nodes_",UP_DOWN,".txt"),sep="\t")
  
  
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% Full_summary$Node_label)

p <- ggplot(data=Full_summary, aes(x=species, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "Number of genes shared with M. truncatula", x = "M. truncatula Nod factor response") +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 18)


pdf(file = paste0(out_path,"../results/Barplot_Medtru_multi_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_Medtru_multi_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()

Summary_count <- data.frame()
i="Mimpud"
DEGs <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0(i,".*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
names(DEGs) <- "Gene_id"
temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
# names(temp) <- c(gsub("./rnaseq/","",i),"HOG")
names(temp) <- c(i,"HOG")
temp[is.na(temp)] <- ""
# temp <- temp[!duplicated(temp[,c('HOG')]),]
cor_Genes_HOG <- dplyr::left_join(temp,Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN)))

names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
NODES <- c(names(ancestors(tree, i)),i)
cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
cor_Genes_HOG$species <- i
cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
Summary_count <- rbind(Summary_count,cor_Genes_HOG)

cor_Genes_HOG <- Mimpud_NF_Up_Recap%>% dplyr::select(Gene_id,HOG,paste0("Node_label_",UP_DOWN))

names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
cor_Genes_HOG$species <- "NF"
cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
Summary_count <- rbind(Summary_count,cor_Genes_HOG)


Summary_count$Node_label <- factor(Summary_count$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% Summary_count$Node_label)

p <- ggplot(data=Summary_count, aes(x=species, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "Number of genes shared with M. pudica", x = "M. pudica Nod factor response") +
  geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
  theme_classic(base_size = 18)


pdf(file = paste0(out_path,"../results/Barplot_Mimpud_multi_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_Mimpud_multi_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()

Summary_count <- data.frame()
i="Lotjap"
DEGs <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0(i,".*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
names(DEGs) <- "Gene_id"
temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
# names(temp) <- c(gsub("./rnaseq/","",i),"HOG")
names(temp) <- c(i,"HOG")
temp[is.na(temp)] <- ""
# temp <- temp[!duplicated(temp[,c('HOG')]),]
cor_Genes_HOG <- dplyr::left_join(temp,Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN)))

names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
NODES <- c(names(ancestors(tree, i)),i)
cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
cor_Genes_HOG$species <- i
cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
Summary_count <- rbind(Summary_count,cor_Genes_HOG)

cor_Genes_HOG <- Lotjap_NF_Up_Recap%>% dplyr::select(Gene_id,HOG,paste0("Node_label_",UP_DOWN))

names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
cor_Genes_HOG$species <- "NF"
cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
Summary_count <- rbind(Summary_count,cor_Genes_HOG)

Full_summary <- Summary_count

Clones_sum <- data.frame(species=unique(Full_summary$species),Sum=0)
for(sp in unique(Full_summary$species)){
  Clones_sum$Sum[Clones_sum$species==sp] <- sum(Full_summary$Overlap[Full_summary$species==sp])
  
}

Categories <- unique(Full_summary$Node_label)
for(Node in Categories){
  CloneNode <- Full_summary$Overlap[Full_summary$species=="NF" & Full_summary$Node_label==Node]
  CtaiNode <- Full_summary$Overlap[Full_summary$species==i & Full_summary$Node_label==Node]
  M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$species=="NF"]-CloneNode,
                                                Clones_sum$Sum[Clones_sum$species==i]-CtaiNode)))
  test <- chisq.test(M)
  fishertest <- fisher.test(M)
  Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- test$p.value
  Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- fishertest$p.value
  Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$species=="NF"] <- fishertest$estimate
  Full_summary$Total[Full_summary$Node_label=="Mimpud_HOG" & Full_summary$species=="NF"] <- Clones_sum$Sum[Clones_sum$species=="NF"]
}

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than WT"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than WT"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than WT"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than WT"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than WT"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than WT"

fwrite(Full_summary,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Lotjap_NodFactor_multi_Nodes_",UP_DOWN,".txt"),sep="\t")


Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% Full_summary$Node_label)

p <- ggplot(data=Full_summary, aes(x=species, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "Number of genes shared with L. japonicus", x = "L. japonicus Nod factor response") +
  geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  theme_classic(base_size = 18)

pdf(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()

}


