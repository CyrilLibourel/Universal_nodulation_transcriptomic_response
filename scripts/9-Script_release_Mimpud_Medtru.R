# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","reshape2","ComplexHeatmap","ggVennDiagram","mixOmics",
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
Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZones",UP_DOWN,".tsv"))

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


# Medtru Mimpud 
Multi_sp_traits <- Recap_HOG %>% dplyr::filter(FI==UP_DOWN | FIIp==UP_DOWN | FIId==UP_DOWN | IZ==UP_DOWN | ZIII==UP_DOWN |
                                                 Mimpud_Organogenesis==UP_DOWN | Mimpud_Release==UP_DOWN | Mimpud_Persist==UP_DOWN |
                                                 Mimpud_NFix==UP_DOWN) %>%
  dplyr::select(HOG, starts_with("Node_label"),FI, FIId, FIIp, IZ, ZIII,Mimpud_Organogenesis, Mimpud_Release, Mimpud_Persist, Mimpud_NFix)
names(Multi_sp_traits) <- c("HOG","Node_label","FI","FIId","FIIp","IZ","ZIII","Mimpud_Organogenesis","Mimpud_Release","Mimpud_Persist","Mimpud_NFix")

Node="N1"
N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)

temp <- N1_data_HOG
temp[temp==UP_DOWN] <- 1
temp[temp==""] <- 0
dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
rownames(dat) <- N1_data_HOG$HOG

upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                   order.by = "freq",keep.order = TRUE, 
                   mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
dev.off()

pdf(file=paste0("./results/upsetR_",Node,"_Node_Mimpud_traits_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
print(upsetPlot)
grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
dev.off()


# Medtru Mimpud Release
Multi_sp_traits <- Recap_HOG %>% dplyr::filter(FIId==UP_DOWN | Mimpud_Release==UP_DOWN | Mimpud_Persist==UP_DOWN) %>%
  dplyr::select(HOG, starts_with("Node_label"), FIId, Mimpud_Release, Mimpud_Persist)
names(Multi_sp_traits) <- c("HOG","Node_label","FIId","Mimpud_Release","Mimpud_Persist")

Node="N1"
N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)

temp <- N1_data_HOG
temp[temp==UP_DOWN] <- 1
temp[temp==""] <- 0
dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
rownames(dat) <- N1_data_HOG$HOG

upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                   order.by = "freq",keep.order = TRUE, 
                   mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
dev.off()

pdf(file=paste0("./results/upsetR_",Node,"_Node_Medru_Mimpud_intracel_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
print(upsetPlot)
grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
dev.off()


treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

Mimpud_nodes <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
Medtru_nodes <- c(names(ancestors(tree, "Medtru")),"Medtru")
Lotjap_nodes <- c(names(ancestors(tree, "Lotjap")),"Lotjap")


# 2 sp
MpMt <- intersect(Medtru_nodes,Mimpud_nodes)
data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MpMt)
data_HOG[data_HOG=="Up"]<- 1
data_HOG[data_HOG==""]  <- 0
temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
  rowwise() %>%
  mutate(SumMt = sum(c(FIId)),
         SumMp = sum(c(Mimpud_Release,Mimpud_Persist)))
temp$SumMt[temp$SumMt==2] <- 1
temp$SumMp[temp$SumMp==2] <- 1
temp <- temp %>% mutate_if(is.character, as.numeric) %>%
  rowwise() %>%
  mutate(Sum = sum(c(SumMt,SumMp))) %>% dplyr::select(-SumMt,-SumMp)
temp$Node_label <- data_HOG$Node_label
temp <- dplyr::filter(temp, Sum>1)

summary_2sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
summary_2sp$Group <- "2species"

summary <- summary_2sp


summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)

p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "", x = "Intracellular response") +
  geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 16)


pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_Release_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_shared_Nodes_Release_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()


# Medtru Mimpud Nod Factor
Multi_sp_traits <- Recap_HOG %>% dplyr::filter(Mimpud_NodFactor==UP_DOWN | Lotjap_NodFactor==UP_DOWN | Medtru_NodFactor==UP_DOWN) %>%
  dplyr::select(HOG, starts_with("Node_label"),Mimpud_NodFactor, Lotjap_NodFactor, Medtru_NodFactor)
names(Multi_sp_traits) <- c("HOG","Node_label","Mimpud_NodFactor","Lotjap_NodFactor","Medtru_NodFactor")

# Node="N1"
# N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)
# 
# temp <- N1_data_HOG
# temp[temp==UP_DOWN] <- 1
# temp[temp==""] <- 0
# dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
# rownames(dat) <- N1_data_HOG$HOG
# 
# upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
#                    order.by = "freq",keep.order = TRUE, 
#                    mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
# dev.off()
# 
# pdf(file=paste0("./results/upsetR_",Node,"_Node_Medru_Mimpud_NodFactor_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
# print(upsetPlot)
# grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
# dev.off()

treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

Mimpud_nodes <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
Medtru_nodes <- c(names(ancestors(tree, "Medtru")),"Medtru")
Lotjap_nodes <- c(names(ancestors(tree, "Lotjap")),"Lotjap")

# 3 sp
MpMtLj <- intersect(intersect(Mimpud_nodes,Medtru_nodes),Lotjap_nodes)
MtLj <- intersect(Medtru_nodes,Lotjap_nodes)
# data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MpMtLj)

threespnodes_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MtLj)
threespnodes_data_HOG$Mimpud_NodFactor[threespnodes_data_HOG$Node_label %in% setdiff(MtLj,MpMtLj)] <- ""
fwrite(threespnodes_data_HOG,"./results/Test_NF.txt")

temp <- threespnodes_data_HOG
temp[temp==UP_DOWN] <- 1
temp[temp==""] <- 0
# temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
rownames(dat) <- threespnodes_data_HOG$HOG

upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                   order.by = "freq",keep.order = TRUE, 
                   mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
dev.off()

pdf(file=paste0("./results/upsetR_Multi_Nodes_Medru_Mimpud_Lotjap_NodFactor_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
print(upsetPlot)
grid.text(paste0("Overlaping NF ",UP_DOWN," HOG among species."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
dev.off()


temp <- temp  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
  rowwise() %>%
  mutate(Sum = sum(c(Mimpud_NodFactor,Lotjap_NodFactor,Medtru_NodFactor)))
temp$Node_label <- theespnodes_data_HOG$Node_label
# temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
temp <- dplyr::filter(temp, Sum>1)
summary_3sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
summary_3sp$Group <- "3species"

summary <- summary_3sp


summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)

p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "", x = "Nod factor response") +
  geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 16)



pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()


data_HOG[data_HOG=="Up"]<- 1
data_HOG[data_HOG==""]  <- 0
temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
  rowwise() %>%
  mutate(Sum = sum(c(Mimpud_NodFactor,Lotjap_NodFactor,Medtru_NodFactor)))
temp$Node_label <- data_HOG$Node_label
temp <- dplyr::filter(temp, Sum>1)

summary_3sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))

# 2 sp
MtLj <- intersect(Medtru_nodes,Lotjap_nodes)
data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MtLj)
data_HOG[data_HOG=="Up"]<- 1
data_HOG[data_HOG==""]  <- 0
temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label,-Mimpud_NodFactor) %>% mutate_if(is.character, as.numeric) %>%
  rowwise() %>%
  mutate(Sum = sum(c(Lotjap_NodFactor,Medtru_NodFactor)))
temp$Node_label <- data_HOG$Node_label
temp <- dplyr::filter(temp, Sum>1)

temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
summary_2sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))


summary_3sp$Group <- "3species"
summary_2sp$Group <- "2species"

summary <- rbind(summary_3sp,summary_2sp)

summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                   "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                      Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                      mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                              "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)

p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
  labs(y= "", x = "Nod factor response") +
  geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 16)



pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()


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
  
  cor_Genes_HOG <- Medtru_NF_Up_Recap %>% dplyr::select(Gene_id,HOG,paste0("Node_label_",UP_DOWN))
  
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
  labs(y= "Number of genes shared with M. truncatula", x = "M. truncatula Nod factor response") +
  geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
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
  labs(y= "Number of genes shared with L. japonicus", x = "L. japonicus Nod factor response") +
  geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
  theme_classic(base_size = 18)

pdf(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()

svg(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
print(p)
dev.off()

}


######## Release Nod factor shared intra specific ########
UP_or_DOWN <- c("Up","Down")

for(UP_DOWN in UP_or_DOWN){
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZones",UP_DOWN,".tsv"))
  
  # Medtru Mimpud 
  Multi_sp_traits <- Recap_HOG %>% dplyr::filter(FIIp==UP_DOWN | FIId==UP_DOWN | Mimpud_Release==UP_DOWN | Mimpud_Persist==UP_DOWN) %>%
    dplyr::select(HOG, starts_with("Node_label"), FIId, Mimpud_Release, Mimpud_Persist)
  names(Multi_sp_traits) <- c("HOG","Node_label","FIId","Mimpud_Release","Mimpud_Persist")
  
  Node="N1"
  N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)
  
  temp <- N1_data_HOG
  temp[temp==UP_DOWN] <- 1
  temp[temp==""] <- 0
  dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
  rownames(dat) <- N1_data_HOG$HOG
  
  upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                     order.by = "freq",keep.order = TRUE, 
                     mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
  dev.off()
  
  pdf(file=paste0("./results/upsetR_",Node,"_Node_Mimpud_traits_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
  print(upsetPlot)
  grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
  dev.off()
  
  
  # Medtru Mimpud Release
  Multi_sp_traits <- Recap_HOG %>% dplyr::filter(FIIp==UP_DOWN | FIId==UP_DOWN | Mimpud_Release==UP_DOWN | Mimpud_Persist==UP_DOWN) %>%
    dplyr::select(HOG, starts_with("Node_label"), FIId, FIIp, Mimpud_Release, Mimpud_Persist)
  names(Multi_sp_traits) <- c("HOG","Node_label","FIId","FIIp","Mimpud_Release","Mimpud_Persist")
  
  Node="N1"
  N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)
  
  temp <- N1_data_HOG
  temp[temp==UP_DOWN] <- 1
  temp[temp==""] <- 0
  dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
  rownames(dat) <- N1_data_HOG$HOG
  
  upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                     order.by = "freq",keep.order = TRUE, 
                     mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
  dev.off()
  
  pdf(file=paste0("./results/upsetR_",Node,"_Node_Medru_Mimpud_intracel_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
  print(upsetPlot)
  grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
  dev.off()
  
  
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  tree <- as(treeObj, "phylo4")
  
  Mimpud_nodes <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
  Medtru_nodes <- c(names(ancestors(tree, "Medtru")),"Medtru")
  Lotjap_nodes <- c(names(ancestors(tree, "Lotjap")),"Lotjap")
  
  
  # 2 sp
  MpMt <- intersect(Medtru_nodes,Mimpud_nodes)
  data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MpMt)
  data_HOG[data_HOG=="Up"]<- 1
  data_HOG[data_HOG==""]  <- 0
  temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(SumMt = sum(c(FIId,FIIp)),
           SumMp = sum(c(Mimpud_Release,Mimpud_Persist)))
  temp$SumMt[temp$SumMt==2] <- 1
  temp$SumMp[temp$SumMp==2] <- 1
  temp <- temp %>% mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(Sum = sum(c(SumMt,SumMp))) %>% dplyr::select(-SumMt,-SumMp)
  temp$Node_label <- data_HOG$Node_label
  temp <- dplyr::filter(temp, Sum>1)
  
  summary_2sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
  summary_2sp$Group <- "2species"
  
  summary <- summary_2sp
  
  
  summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
  mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
             "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
  LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                     "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                        Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                        mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                                "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
  LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)
  
  p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
    labs(y= "", x = "Intracellular response") +
    geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  
  
  pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_Release_",UP_DOWN,".pdf"), width=12, height=8)
  print(p)
  dev.off()
  
  svg(file = paste0(out_path,"../results/Barplot_shared_Nodes_Release_",UP_DOWN,".svg"), width=12, height=8)
  print(p)
  dev.off()
  
  
  # Medtru Mimpud Nod Factor
  Multi_sp_traits <- Recap_HOG %>% dplyr::filter(Mimpud_NodFactor ==UP_DOWN | Lotjap_NodFactor ==UP_DOWN | Medtru_NodFactor ==UP_DOWN) %>%
    dplyr::select(HOG, starts_with("Node_label"),Mimpud_NodFactor, Lotjap_NodFactor, Medtru_NodFactor)
  names(Multi_sp_traits) <- c("HOG","Node_label","Mimpud_NodFactor","Lotjap_NodFactor","Medtru_NodFactor")
  
  Node="N1"
  N1_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label==Node)
  
  temp <- N1_data_HOG
  temp[temp==UP_DOWN] <- 1
  temp[temp==""] <- 0
  dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
  rownames(dat) <- N1_data_HOG$HOG
  
  upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                     order.by = "freq",keep.order = TRUE, 
                     mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
  dev.off()
  
  pdf(file=paste0("./results/upsetR_",Node,"_Node_Medru_Mimpud_NodFactor_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
  print(upsetPlot)
  grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
  dev.off()
  
  
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  tree <- as(treeObj, "phylo4")
  
  Mimpud_nodes <- c(names(ancestors(tree, "Mimpud")),"Mimpud")
  Medtru_nodes <- c(names(ancestors(tree, "Medtru")),"Medtru")
  Lotjap_nodes <- c(names(ancestors(tree, "Lotjap")),"Lotjap")
  
  # 3 sp
  MpMtLj <- intersect(intersect(Mimpud_nodes,Medtru_nodes),Lotjap_nodes)
  MtLj <- intersect(Medtru_nodes,Lotjap_nodes)
  data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MpMtLj)
  
  theespnodes_data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MtLj)
  
  temp <- theespnodes_data_HOG
  temp[temp==UP_DOWN] <- 1
  temp[temp==""] <- 0
  # temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
  dat <- as.data.frame(dplyr::select(temp,-HOG, -Node_label)%>% mutate_if(is.character, as.numeric))
  rownames(dat) <- theespnodes_data_HOG$HOG
  
  upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                     order.by = "freq",keep.order = TRUE, 
                     mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
  dev.off()
  
  pdf(file=paste0("./results/upsetR_Multi_Nodes_Medru_Mimpud_Lotjap_NodFactor_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
  print(upsetPlot)
  grid.text(paste0("Overlaping NF ",UP_DOWN," HOG among species."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
  dev.off()
  
  
  temp <- temp  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(Sum = sum(c(Mimpud_NodFactor,Lotjap_NodFactor,Medtru_NodFactor)))
  temp$Node_label <- theespnodes_data_HOG$Node_label
  # temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
  temp <- dplyr::filter(temp, Sum>1)
  summary_3sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
  summary_3sp$Group <- "3species"
  
  summary <- summary_3sp
  
  
  summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
  mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
             "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
  LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                     "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                        Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                        mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                                "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
  LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)
  
  p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
    labs(y= "", x = "Nod factor response") +
    geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  
  
  
  pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
  print(p)
  dev.off()
  
  
  data_HOG[data_HOG=="Up"]<- 1
  data_HOG[data_HOG==""]  <- 0
  temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label) %>% mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(Sum = sum(c(Mimpud_NodFactor,Lotjap_NodFactor,Medtru_NodFactor)))
  temp$Node_label <- data_HOG$Node_label
  temp <- dplyr::filter(temp, Sum>1)
  
  summary_3sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
  
  # 2 sp
  MtLj <- intersect(Medtru_nodes,Lotjap_nodes)
  data_HOG <- Multi_sp_traits %>% dplyr::filter(Node_label %in% MtLj)
  data_HOG[data_HOG=="Up"]<- 1
  data_HOG[data_HOG==""]  <- 0
  temp <- data_HOG  %>% dplyr::select(-HOG, -Node_label,-Mimpud_NodFactor) %>% mutate_if(is.character, as.numeric) %>%
    rowwise() %>%
    mutate(Sum = sum(c(Lotjap_NodFactor,Medtru_NodFactor)))
  temp$Node_label <- data_HOG$Node_label
  temp <- dplyr::filter(temp, Sum>1)
  
  # temp$Node_label[temp$Node_label %in% setdiff(MtLj,MpMtLj)] <- "N7"
  summary_2sp <- temp %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Sum))
  
  
  summary_3sp$Group <- "3species"
  summary_2sp$Group <- "2species"
  
  summary <- rbind(summary_3sp,summary_2sp)
  
  summary$Node_label <- factor(summary$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
  mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
             "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")
  LEGENDS <- data.frame(Node_label=c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Papilionoideae clade",
                                     "Dalbergioid+Hologalegina clade","Rosales clade","Hologalegina clade","Dalbergioid clade", "Species specific"),
                        Node_number=c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"),
                        mycol=c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N7"="darkseagreen2","N12" = "darkseagreen3",
                                "N15"="coral1","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00"))
  LEGENDS <- dplyr::filter(LEGENDS,Node_number %in% summary$Node_label)
  
  p <- ggplot(data=summary, aes(x=Group, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = LEGENDS$mycol,labels = LEGENDS$Node_label) + 
    labs(y= "", x = "Nod factor response") +
    geom_text(data = summary, aes(y = Overlap, label = Overlap),size=6,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  
  
  
  pdf(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
  print(p)
  dev.off()
  
  svg(file = paste0(out_path,"../results/Barplot_shared_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
  print(p)
  dev.off()
  
  
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
    labs(y= "Number of genes shared with M. truncatula", x = "M. truncatula Nod factor response") +
    geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
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
    labs(y= "Number of genes shared with L. japonicus", x = "L. japonicus Nod factor response") +
    geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".pdf"), width=12, height=8)
  print(p)
  dev.off()
  
  svg(file = paste0(out_path,"../results/Barplot_Lotjap_multi_Nodes_NF_",UP_DOWN,".svg"), width=12, height=8)
  print(p)
  dev.off()
  
}



######## Release intra specific ########

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
HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
UP_or_DOWN <- c("Up")
UP_DOWN <- c("Up")
for(UP_DOWN in UP_or_DOWN){
  # Mimpud
  Mimpud_cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))
  Mimpud_NODES <- c("N1","N2","Mimpud")
  
  SignalP <- fread("../ncr/output.gff3",h=F)
  SignalP <- SignalP %>% filter(str_detect(V1, "^Mimpud_")) %>% dplyr::select(V1,V3)
  names(SignalP) <- c("Gene_id","SignalP")
  
  # Mimpud Nod
  Mimpud_AllUp <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(Mimpud_AllUp) <- "Gene_id"
  Mimpud_AllUp <- left_join(Mimpud_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  Mimpud_AllUp <- data.frame(Gene_id=unique(Mimpud_AllUp$Gene_id[!is.na(Mimpud_AllUp$HOG)]),Mimpud="Up")
  Mimpud_AllUp <- left_join(Mimpud_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  
  HOG_id_NODE <- Nodes_Up %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN))
  names(HOG_id_NODE) <- c("HOG","Node_label")
  HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% Mimpud_NODES)
  Mimpud_NODE_genes <- left_join(Mimpud_AllUp,HOG_id_NODE,by="HOG")
  
  # Specific HOG
  HOG_with_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_")) %>% dplyr::select(HOG)
  HOG_with_Mimpud <- unique(HOG_with_Mimpud$HOG)
  
  HOG_without_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Mimpud <- unique(HOG_without_Mimpud$HOG)
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, HOG_without_Mimpud))
  Mimpud_spe_NODE_genes <- Mimpud_NODE_genes %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  Mimpud_spe_NODE_genes$Node_label <- "Species_specific"
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, Mimpud_spe_NODE_genes$HOG))
  Mimpud_spe <- Mimpud_NODE_genes %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  
  Mimpud_NODE_genes <- rbind(Mimpud_spe_NODE_genes, Mimpud_spe)

  TRAITS = c("Release","Persist")
  for(trait in TRAITS){
    trait_df <- fread(input = paste0("./trait_analysis/",trait,".txt"),header=TRUE)
    names(trait_df)[1] <- "Gene_id"
    trait_df <- trait_df %>% dplyr::filter(Up_Down==paste0(UP_DOWN)) %>% dplyr::select(Gene_id)
    trait_df$trait <- "Up";names(trait_df) <- c("Gene_id",trait)
    Mimpud_NODE_genes <- dplyr::left_join(Mimpud_NODE_genes,trait_df,by="Gene_id")
    Mimpud_NODE_genes <- Mimpud_NODE_genes %>% dplyr::distinct()
  }

  Mimpud_NODE_genes <- left_join(Mimpud_NODE_genes,SignalP,by="Gene_id")

#### Release ####
  # PLS-DA
  # Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_infos.txt"), sep="\t", h=TRUE, row.names=1)
  Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
  Mimpud_prot <- dplyr::select(Mimpud_prot,-Acronym)
  Mimpud_prot$SignalP[Mimpud_prot$SignalP=="signal_peptide"] <- 1
  Mimpud_prot$SignalP[Mimpud_prot$SignalP==""] <- 0; Mimpud_prot$SignalP <- as.numeric(Mimpud_prot$SignalP)
  
  X<-data.matrix(t(Mimpud_prot), rownames.force = NA)
  colors<-rep(c("black","blue1"))
  colors2<-rep(c("black","blue1"), c(10,28))
  Y = data.frame(Gene_id=row.names(Mimpud_prot))
  TRAIT_df <- Mimpud_NODE_genes %>% dplyr::filter(Release=="Up") %>% dplyr::select(Gene_id,Release,Node_label)
  TRAIT_df <- TRAIT_df %>% dplyr::filter(Node_label=="Species_specific") %>% dplyr::select(Gene_id,Release)
  Y <- dplyr::left_join(Y,TRAIT_df)
  Y[is.na(Y)] <- "noRelease"; Y$Release[Y$Release=="Up"] <- "Release"
  Y = as.factor(Y$Release) #les 10 premièeres colonnes sont des sp LFA, les 28 suivantes sont des NLFA
  
  # PLS DA
  res_plsda <- plsda(t(X), Y=Y, ncomp=4, near.zero.var=FALSE) 
  # plotIndiv(res_plsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
  plotVar(res_plsda, comp=c(1,2),cex=3)
  plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Release species specific") # graphe de contribution des IPR à la séparation des LFA et nLFA 
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Loadings_Mimpud_Release_",UP_DOWN,"_aminoacids.pdf"), width=12, height=8)
  plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Release species specific")
  dev.off()
  
  Biplot_data_signalP <- Mimpud_prot %>% dplyr::filter(SignalP=="1") %>% dplyr::select(Pep_length,P_prop,H_prop)
  Biplot_data_signalP$Gene_id <- row.names(Biplot_data_signalP)
  Biplot_data_signalP <- dplyr::left_join(Biplot_data_signalP,TRAIT_df)
  Biplot_data_signalP$Release[Biplot_data_signalP$Release=="Up"] <- "Release"
  Biplot_data_signalP$Release[is.na(Biplot_data_signalP$Release)] <- "noRelease"

 # Signal P Proline proportion and size
  scatterPlot <- ggplot(Biplot_data_signalP,aes(P_prop, Pep_length, color=Release)) + 
    geom_point() + expand_limits(x=0)+
    scale_color_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
    geom_vline(xintercept = 10)+
    geom_hline(yintercept = 200)
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Biplot_data_signalP, aes(P_prop, fill=Release)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")

  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=Release)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")

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
          axis.ticks = element_blank()
    )
  pdffile = paste0("./results/PLS-DA/Mimpud_Release_signalP_Proline_length.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=15,width=15)
  # par(mar=c(5,5,5,5))
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=200&Biplot_data_signalP$P_prop>=10] <- "Long_HighP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<200&Biplot_data_signalP$P_prop>=10] <- "Short_HighP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<200&Biplot_data_signalP$P_prop<10] <- "Short_LowP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=200&Biplot_data_signalP$P_prop<10] <- "Long_LowP"
  
  
  HOG_all_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(Release=="noRelease") %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "noRelease"
  trait_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(Release=="Release") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_df$Trait <- "Release"
  HOG_all_df <- rbind(HOG_all_df,trait_df)

  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  
  Nodes <- unique(Full_summary$Class)
  i="Release"
  for(Node in Nodes){
      CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
      if(length(CloneNode)==0){CloneNode<-0}
      CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noRelease" & Full_summary$Class==Node]
      M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                    Clones_sum$Sum[Clones_sum$Trait=="noRelease"]-CtaiNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
      Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
      Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
  }
  
  Full_summary$Total[Full_summary$Trait=="noRelease" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noRelease"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Trait[Full_summary$Trait=="noRelease"] <- "Ctai"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Ctai","Release"))
 
  Full_summary_Up <- Full_summary
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP")) + 
    labs(y= "Propotion of proteins", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Mimpud_Release_signalP_Proline_length.pdf"), width=18, height=8)
  print(p)
  dev.off()
  
#### Without signal P Proline proportion and size ####
  Biplot_data_nosignalP <- Mimpud_prot %>% dplyr::filter(SignalP=="0") %>% dplyr::select(Pep_length,P_prop,H_prop)
  Biplot_data_nosignalP$Gene_id <- row.names(Biplot_data_nosignalP)
  Biplot_data_nosignalP <- dplyr::left_join(Biplot_data_nosignalP,TRAIT_df)
  Biplot_data_nosignalP$Release[Biplot_data_nosignalP$Release=="Up"] <- "Release"
  Biplot_data_nosignalP$Release[is.na(Biplot_data_nosignalP$Release)] <- "noRelease"
  
  scatterPlot <- ggplot(Biplot_data_nosignalP,aes(P_prop, Pep_length, color=Release)) + 
    geom_point() + expand_limits(x=0)+
    scale_color_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
    geom_vline(xintercept = 10)+
    geom_hline(yintercept = 200)
  
  scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Biplot_data_nosignalP, aes(P_prop, fill=Release)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=Release)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
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
          axis.ticks = element_blank()
    )
  pdffile = paste0("./results/PLS-DA/Mimpud_Release_withoutsignalP_Proline_length.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=15,width=15)
  # par(mar=c(5,5,5,5))
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=200&Biplot_data_nosignalP$P_prop>=10] <- "Long_HighP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<200&Biplot_data_nosignalP$P_prop>=10] <- "Short_HighP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<200&Biplot_data_nosignalP$P_prop<10] <- "Short_LowP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=200&Biplot_data_nosignalP$P_prop<10] <- "Long_LowP"
  
  
  HOG_all_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(Release=="noRelease") %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "noRelease"
  trait_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(Release=="Release") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_df$Trait <- "Release"
  HOG_all_df <- rbind(HOG_all_df,trait_df)
  
  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    }
  
  Nodes <- unique(Full_summary$Class)
  i="Release"
  for(Node in Nodes){
    CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
    if(length(CloneNode)==0){CloneNode<-0}
    CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noRelease" & Full_summary$Class==Node]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Trait=="noRelease"]-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
    Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
  }
  
  Full_summary$Total[Full_summary$Trait=="noRelease" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noRelease"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Trait[Full_summary$Trait=="noRelease"] <- "Ctai"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Ctai","Release"))
  
  Full_summary_Up <- Full_summary
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP")) + 
    labs(y= "Propotion of proteins", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Mimpud_Release_withoutsignalP_Proline_length.pdf"), width=18, height=8)
  print(p)
  dev.off()
  
#### Persist ####
  X<-data.matrix(t(Mimpud_prot), rownames.force = NA)
  colors<-rep(c("black","blue1"))
  colors2<-rep(c("black","blue1"), c(10,28))
  Y = data.frame(Gene_id=row.names(Mimpud_prot))
  TRAIT_df <- Mimpud_NODE_genes %>% dplyr::filter(Persist=="Up") %>% dplyr::select(Gene_id,Persist,Node_label)
  TRAIT_df <- TRAIT_df %>% dplyr::filter(Node_label=="Species_specific") %>% dplyr::select(Gene_id,Persist)
  Y <- dplyr::left_join(Y,TRAIT_df)
  Y[is.na(Y)] <- "noPersist"; Y$Persist[Y$Persist=="Up"] <- "Persist"
  Y = as.factor(Y$Persist) #les 10 premièeres colonnes sont des sp LFA, les 28 suivantes sont des NLFA
  
  # PLS DA
  res_plsda <- plsda(t(X), Y=Y, ncomp=4, near.zero.var=FALSE) 
  # plotIndiv(res_plsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
  plotVar(res_plsda, comp=c(1,2),cex=3)
  plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Persist species specific") # graphe de contribution des IPR à la séparation des LFA et nLFA 
  
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Loadings_Mimpud_Persist_",UP_DOWN,"_aminoacids.pdf"), width=12, height=8)
  plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Mimpud Persist species specific")
  dev.off()
  
  Biplot_data_signalP <- Mimpud_prot %>% dplyr::filter(SignalP=="1") %>% dplyr::select(Pep_length,P_prop,H_prop)
  Biplot_data_signalP$Gene_id <- row.names(Biplot_data_signalP)
  Biplot_data_signalP <- dplyr::left_join(Biplot_data_signalP,TRAIT_df)
  Biplot_data_signalP$Persist[Biplot_data_signalP$Persist=="Up"] <- "Persist"
  Biplot_data_signalP$Persist[is.na(Biplot_data_signalP$Persist)] <- "noPersist"
  
  # Signal P Proline proportion and size
  scatterPlot <- ggplot(Biplot_data_signalP,aes(P_prop, Pep_length, color=Persist)) + 
    geom_point() + expand_limits(x=0)+
    scale_color_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
    geom_vline(xintercept = 10)+
    geom_hline(yintercept = 200)
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Biplot_data_signalP, aes(P_prop, fill=Persist)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=Persist)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
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
          axis.ticks = element_blank()
    )
  pdffile = paste0("./results/PLS-DA/Mimpud_Persist_signalP_Proline_length.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=15,width=15)
  # par(mar=c(5,5,5,5))
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=200&Biplot_data_signalP$P_prop>=10] <- "Long_HighP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<200&Biplot_data_signalP$P_prop>=10] <- "Short_HighP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<200&Biplot_data_signalP$P_prop<10] <- "Short_LowP"
  Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=200&Biplot_data_signalP$P_prop<10] <- "Long_LowP"
  
  
  HOG_all_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(Persist=="noPersist") %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "noPersist"
  trait_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(Persist=="Persist") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_df$Trait <- "Persist"
  HOG_all_df <- rbind(HOG_all_df,trait_df)
  
  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  
  Nodes <- unique(Full_summary$Class)
  i="Persist"
  for(Node in Nodes){
    CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
    if(length(CloneNode)==0){CloneNode<-0}
    CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noPersist" & Full_summary$Class==Node]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Trait=="noPersist"]-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
    Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
  }
  
  Full_summary$Total[Full_summary$Trait=="noPersist" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noPersist"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Trait[Full_summary$Trait=="noPersist"] <- "Ctai"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Ctai","Persist"))
  
  Full_summary_Up <- Full_summary
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP")) + 
    labs(y= "Propotion of proteins", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Mimpud_Persist_signalP_Proline_length.pdf"), width=18, height=8)
  print(p)
  dev.off()

#### Without signal P Proline proportion and size ####
  Biplot_data_nosignalP <- Mimpud_prot %>% dplyr::filter(SignalP=="0") %>% dplyr::select(Pep_length,P_prop,H_prop)
  Biplot_data_nosignalP$Gene_id <- row.names(Biplot_data_nosignalP)
  Biplot_data_nosignalP <- dplyr::left_join(Biplot_data_nosignalP,TRAIT_df)
  Biplot_data_nosignalP$Persist[Biplot_data_nosignalP$Persist=="Up"] <- "Persist"
  Biplot_data_nosignalP$Persist[is.na(Biplot_data_nosignalP$Persist)] <- "noPersist"
  
  scatterPlot <- ggplot(Biplot_data_nosignalP,aes(P_prop, Pep_length, color=Persist)) + 
    geom_point() + expand_limits(x=0)+
    scale_color_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
    geom_vline(xintercept = 10)+
    geom_hline(yintercept = 200)
  
  scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Biplot_data_nosignalP, aes(P_prop, fill=Persist)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=Persist)) + 
    geom_density(alpha=.5) + 
    scale_fill_manual(values = c('#999999','#E69F00')) + 
    theme(legend.position = "none")
  
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
          axis.ticks = element_blank()
    )
  pdffile = paste0("./results/PLS-DA/Mimpud_Persist_withoutsignalP_Proline_length.pdf")
  #adjust the size and margin of the pdf file
  pdf(pdffile,height=15,width=15)
  # par(mar=c(5,5,5,5))
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=200&Biplot_data_nosignalP$P_prop>=10] <- "Long_HighP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<200&Biplot_data_nosignalP$P_prop>=10] <- "Short_HighP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<200&Biplot_data_nosignalP$P_prop<10] <- "Short_LowP"
  Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=200&Biplot_data_nosignalP$P_prop<10] <- "Long_LowP"
  
  
  HOG_all_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(Persist=="noPersist") %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "noPersist"
  trait_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(Persist=="Persist") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_df$Trait <- "Persist"
  HOG_all_df <- rbind(HOG_all_df,trait_df)
  
  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  
  Nodes <- unique(Full_summary$Class)
  i="Persist"
  for(Node in Nodes){
    CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
    if(length(CloneNode)==0){CloneNode<-0}
    CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noPersist" & Full_summary$Class==Node]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Trait=="noPersist"]-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
    Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
  }
  
  Full_summary$Total[Full_summary$Trait=="noPersist" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noPersist"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Trait[Full_summary$Trait=="noPersist"] <- "Ctai"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Ctai","Persist"))
  
  Full_summary_Up <- Full_summary
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("Long_HighP","Long_LowP","Short_LowP","Short_HighP")) + 
    labs(y= "Propotion of proteins", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Mimpud_Persist_withoutsignalP_Proline_length.pdf"), width=18, height=8)
  print(p)
  dev.off()
  
  
  # sparse PLS DA
  # list.keepX = c(5, 5, 11) # garder 100 IPR les + explicatifs pour la séparation des LFA et nLFA (100 pour chaque composante)
  # res_splsda <- splsda(t(X), Y=Y, ncomp=20, near.zero.var=FALSE, keepX = list.keepX) ##reduit nb de var qui peuvent construire les composantes (100 chacune cf keepX.)
  # plotIndiv(res_splsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
  # plotVar(res_splsda, comp=c(2:3),cex=3)
  # barplot((res_splsda$prop_expl_var$X)*100)
  
  # selected_protcarac_1 <- selectVar(res_splsda, comp=1) # protcarac sélectionnés pour la première composante
  # selected_protcarac_2 <- selectVar(res_splsda, comp=2) # protcarac sélectionnés pour la deuxième composante
  # selected_protcarac_3 <- selectVar(res_splsda, comp=3) # protcarac sélectionnés pour la troisième composante
  # 
  # plotLoadings(res_splsda, contrib = 'max', method = 'mean') # graphe de contribution des IPR à la séparation des LFA et nLFA 
  # plotIndiv(res_splsda, style="3d")
  
  
  # HOG_all_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Mimpud=="Up") %>% dplyr::summarise(Overlap=length(Gene_id))
  # HOG_all_df$Trait <- "Mimpud"
  HOG_all_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Mimpud=="Up" & SignalP=="signal_peptide") %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "Mimpud"
  # trait_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Release=="Up") %>% dplyr::summarise(Overlap=length(Gene_id))
  # trait_df$Trait <- "Release"
  trait_sigP_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Release=="Up" & SignalP=="signal_peptide") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_sigP_df$Trait <- "Release_signalP"
  # HOG_all_df <- rbind(HOG_all_df,trait_df)
  HOG_all_df <- rbind(HOG_all_df,trait_sigP_df)
  
  # trait_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Persist=="Up") %>% dplyr::summarise(Overlap=length(Gene_id))
  # trait_df$Trait <- "Persist"
  trait_sigP_df <- Mimpud_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Persist=="Up" & SignalP=="signal_peptide") %>% dplyr::summarise(Overlap=length(Gene_id))
  trait_sigP_df$Trait <- "Persist_signalP"
  # HOG_all_df <- rbind(HOG_all_df,trait_df)
  HOG_all_df <- rbind(HOG_all_df,trait_sigP_df)
  
  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  Mimpud_NODES <- unique(Full_summary$Node_label)
  for(Node in Mimpud_NODES){
    for(i in unique(Full_summary$Trait)){
      CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Node_label ==Node]
      ControlNode <- Full_summary$Overlap[Full_summary$Trait=="Mimpud" & Full_summary$Node_label ==Node]
      M <- as.table(rbind(c(CloneNode, ControlNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                       Clones_sum$Sum[Clones_sum$Trait=="Mimpud"]-ControlNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- test$p.value
      Full_summary$FishTest[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- fishertest$estimate
      Full_summary$Total[Full_summary$Node_label =="N1" & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
    }
  }
  
  Full_summary$Total[Full_summary$Trait=="Mimpud" & Full_summary$Node_label =="N1"] <- Clones_sum$Sum[Clones_sum$Trait=="Mimpud"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than Mimpud"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than Mimpud"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than Mimpud"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than Mimpud"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than Mimpud"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than Mimpud"
  
  Mimpud_Full_summary <- Full_summary
  
# Medtru
Medtru_cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, "^Medtru_"))
Medtru_NODES <- c("N1","N2","N7","N12","N21","Medtru")

FIId <- data.frame(Gene_id=fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_FIId_LogFC1.5.txt"),h=FALSE)$V1,FIId="Up")

SignalP <- fread("../ncr/output.gff3",h=F)
SignalP <- SignalP %>% filter(str_detect(V1, "^Medtru_")) %>% dplyr::select(V1,V3)
names(SignalP) <- c("Gene_id","SignalP")

# Medtru Nod
Medtru_AllUp <- fread(input = paste0("./rnaseq/Medtru_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
names(Medtru_AllUp) <- "Gene_id"
Medtru_AllUp <- left_join(Medtru_AllUp,Medtru_cor_Genes_HOG,by="Gene_id")
Medtru_AllUp <- data.frame(Gene_id=unique(Medtru_AllUp$Gene_id[!is.na(Medtru_AllUp$HOG)]),Medtru="Up")
Medtru_AllUp <- left_join(Medtru_AllUp,Medtru_cor_Genes_HOG,by="Gene_id")

HOG_id_NODE <- Nodes_Up %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN))
names(HOG_id_NODE) <- c("HOG","Node_label")
HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% Medtru_NODES)
Medtru_NODE_genes <- left_join(Medtru_AllUp,HOG_id_NODE,by="HOG")
Medtru_NODE_genes <- left_join(Medtru_NODE_genes,FIId,by="Gene_id")
Medtru_NODE_genes <- left_join(Medtru_NODE_genes,SignalP,by="Gene_id")

# Specific HOG
HOG_with_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_")) %>% dplyr::select(HOG)
HOG_with_Medtru <- unique(HOG_with_Medtru$HOG)

HOG_without_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_", negate = TRUE)) %>% dplyr::select(HOG)
HOG_without_Medtru <- unique(HOG_without_Medtru$HOG)

Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, HOG_without_Medtru))
Medtru_spe_NODE_genes <- Medtru_NODE_genes %>% dplyr::filter(HOG %in% Medtru_spe$HOG)
Medtru_spe_NODE_genes$Node_label <- "Species_specific"

Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, Medtru_spe_NODE_genes$HOG))
Medtru_spe <- Medtru_NODE_genes %>% dplyr::filter(HOG %in% Medtru_spe$HOG)

Medtru_NODE_genes <- rbind(Medtru_spe_NODE_genes, Medtru_spe)

# PLS-DA
Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_infos.txt"), sep="\t", h=TRUE, row.names=1)
Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
# Medtru_prot <- Medtru_prot %>% dplyr::filter(SignalP=="signal_peptide")
Medtru_prot <- dplyr::select(Medtru_prot,-Acronym)
Medtru_prot$SignalP[Medtru_prot$SignalP=="signal_peptide"] <- 1
Medtru_prot$SignalP[Medtru_prot$SignalP==""] <- 0; Medtru_prot$SignalP <- as.numeric(Medtru_prot$SignalP)

X<-data.matrix(t(Medtru_prot), rownames.force = NA)
colors<-rep(c("black","blue1"))
colors2<-rep(c("black","blue1"), c(10,28))
Y = data.frame(Gene_id=row.names(Medtru_prot))
TRAIT_df <- Medtru_NODE_genes %>% dplyr::filter(FIId=="Up") %>% dplyr::select(Gene_id,FIId,Node_label)
TRAIT_df <- TRAIT_df %>% dplyr::filter(Node_label=="Species_specific") %>% dplyr::select(Gene_id,FIId)
Y <- dplyr::left_join(Y,TRAIT_df)
Y[is.na(Y)] <- "noFIId"; Y$FIId[Y$FIId=="Up"] <- "FIId"
Y = as.factor(Y$FIId) #les 10 premièeres colonnes sont des sp LFA, les 28 suivantes sont des NLFA

# PLS DA
res_plsda <- plsda(t(X), Y=Y, ncomp=4, near.zero.var=FALSE) 
# plotIndiv(res_plsda,ind.names=FALSE, col = colors, comp = c(1,2), cex=3, X.label='PC1', Y.label='PC2',abline.line=TRUE, ellipse=TRUE,star=TRUE, ellipse.level=0.80)
plotVar(res_plsda, comp=c(1,2),cex=3)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Medtru FIId species specific") # graphe de contribution des IPR à la séparation des LFA et nLFA 


pdf(file = paste0(out_path,"../results/PLS-DA/Loadings_Medtru_FIId_",UP_DOWN,"_aminoacids.pdf"), width=12, height=8)
plotLoadings(res_plsda, contrib = 'max', method = 'mean', title="Medtru FIId species specific")
dev.off()
#### FIId ####
Biplot_data_signalP <- Medtru_prot %>% dplyr::filter(SignalP=="1") %>% dplyr::select(Pep_length,C_prop,S_prop)
Biplot_data_signalP$Gene_id <- row.names(Biplot_data_signalP)
Biplot_data_signalP <- dplyr::left_join(Biplot_data_signalP,TRAIT_df)
Biplot_data_signalP$FIId[Biplot_data_signalP$FIId=="Up"] <- "FIId"
Biplot_data_signalP$FIId[is.na(Biplot_data_signalP$FIId)] <- "noFIId"

# Signal P Cysteine proportion and size
scatterPlot <- ggplot(Biplot_data_signalP,aes(C_prop, Pep_length, color=FIId)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 8)+
  geom_hline(yintercept = 75)

# scatterPlot
# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(Biplot_data_signalP, aes(C_prop, fill=FIId)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position = "none")

# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=FIId)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position = "none")

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
        axis.ticks = element_blank()
  )
pdffile = paste0("./results/PLS-DA/Medtru_FIId_signalP_Cysteine_length.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=75&Biplot_data_signalP$C_prop>=6] <- "Long_HighC"
Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<75&Biplot_data_signalP$C_prop>=6] <- "Short_HighC"
Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length<75&Biplot_data_signalP$C_prop<6] <- "Short_LowC"
Biplot_data_signalP$Class[Biplot_data_signalP$Pep_length>=75&Biplot_data_signalP$C_prop<6] <- "Long_LowC"


HOG_all_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(FIId=="noFIId") %>% dplyr::summarise(Overlap=length(Gene_id))
HOG_all_df$Trait <- "noFIId"
trait_df <- Biplot_data_signalP %>% group_by(Class) %>% dplyr::filter(FIId=="FIId") %>% dplyr::summarise(Overlap=length(Gene_id))
trait_df$Trait <- "FIId"
HOG_all_df <- rbind(HOG_all_df,trait_df)

Full_summary <- HOG_all_df

Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
for(i in unique(Full_summary$Trait)){
  Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
  
}

Nodes <- unique(Full_summary$Class)
i="FIId"
for(Node in Nodes){
  CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
  if(length(CloneNode)==0){CloneNode<-0}
  CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noFIId" & Full_summary$Class==Node]
  M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                Clones_sum$Sum[Clones_sum$Trait=="noFIId"]-CtaiNode)))
  test <- chisq.test(M)
  fishertest <- fisher.test(M)
  Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
  Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
  Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
  Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
}

Full_summary$Total[Full_summary$Trait=="noFIId" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noFIId"]

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"

Full_summary$Trait[Full_summary$Trait=="noFIId"] <- "Medtru"
Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighC","Long_LowC","Short_LowC","Short_HighC"))
Full_summary$Trait <- factor(Full_summary$Trait,
                             levels = c("Medtru","FIId"))

Full_summary_Up <- Full_summary
mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = mycol,labels = c("Long_HighC","Long_LowC","Short_LowC","Short_HighC")) + 
  labs(y= "Propotion of proteins", x = "Traits") +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 18)

pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Medtru_FIId_signalP_Cysteine_length.pdf"), width=18, height=8)
print(p)
dev.off()

#### Without signal P Cysteine proportion and size ####
Biplot_data_nosignalP <- Medtru_prot %>% dplyr::filter(SignalP=="0") %>% dplyr::select(Pep_length,P_prop,C_prop)
Biplot_data_nosignalP$Gene_id <- row.names(Biplot_data_nosignalP)
Biplot_data_nosignalP <- dplyr::left_join(Biplot_data_nosignalP,TRAIT_df)
Biplot_data_nosignalP$FIId[Biplot_data_nosignalP$FIId=="Up"] <- "FIId"
Biplot_data_nosignalP$FIId[is.na(Biplot_data_nosignalP$FIId)] <- "noFIId"

scatterPlot <- ggplot(Biplot_data_nosignalP,aes(C_prop, Pep_length, color=FIId)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_hline(yintercept = 75)

scatterPlot
# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(Biplot_data_nosignalP, aes(C_prop, fill=FIId)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position = "none")

# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(Biplot_data_signalP, aes(y=Pep_length, fill=FIId)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#E69F00','#999999')) + 
  theme(legend.position = "none")

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
        axis.ticks = element_blank()
  )
pdffile = paste0("./results/PLS-DA/Medtru_FIId_withoutsignalP_Cysteine_length.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=75&Biplot_data_nosignalP$C_prop>=6] <- "Long_HighC"
Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<75&Biplot_data_nosignalP$C_prop>=6] <- "Short_HighC"
Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length<75&Biplot_data_nosignalP$C_prop<6] <- "Short_LowC"
Biplot_data_nosignalP$Class[Biplot_data_nosignalP$Pep_length>=75&Biplot_data_nosignalP$C_prop<6] <- "Long_LowC"


HOG_all_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(FIId=="noFIId") %>% dplyr::summarise(Overlap=length(Gene_id))
HOG_all_df$Trait <- "noFIId"
trait_df <- Biplot_data_nosignalP %>% group_by(Class) %>% dplyr::filter(FIId=="FIId") %>% dplyr::summarise(Overlap=length(Gene_id))
trait_df$Trait <- "FIId"
HOG_all_df <- rbind(HOG_all_df,trait_df)

Full_summary <- HOG_all_df

Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
for(i in unique(Full_summary$Trait)){
  Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
  
}

Nodes <- unique(Full_summary$Class)
i="FIId"
for(Node in Nodes){
  CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Class==Node]
  if(length(CloneNode)==0){CloneNode<-0}
  CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="noFIId" & Full_summary$Class==Node]
  M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                Clones_sum$Sum[Clones_sum$Trait=="noFIId"]-CtaiNode)))
  test <- chisq.test(M)
  fishertest <- fisher.test(M)
  Full_summary$Chi2[Full_summary$Class==Node & Full_summary$Trait==i] <- test$p.value
  Full_summary$FishTest[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$p.value
  Full_summary$oddsratio[Full_summary$Class==Node & Full_summary$Trait==i] <- fishertest$estimate
  Full_summary$Total[Full_summary$Class==unique(Full_summary$Class)[1] & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
}

Full_summary$Total[Full_summary$Trait=="noFIId" & Full_summary$Class==unique(Full_summary$Class)[1]] <- Clones_sum$Sum[Clones_sum$Trait=="noFIId"]

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"

Full_summary$Trait[Full_summary$Trait=="noFIId"] <- "Medtru"
Full_summary$Class <- factor(Full_summary$Class,levels = c("Long_HighC","Long_LowC","Short_LowC","Short_HighC"))
Full_summary$Trait <- factor(Full_summary$Trait,
                             levels = c("Medtru","FIId"))

Full_summary_Up <- Full_summary
mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Class)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = mycol,labels = c("Long_HighC","Long_LowC","Short_LowC","Short_HighC")) + 
  labs(y= "Propotion of proteins", x = "Traits") +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
  theme_classic(base_size = 18)

pdf(file = paste0(out_path,"../results/PLS-DA/Barplot_Medtru_FIId_withoutsignalP_Cysteine_length.pdf"), width=18, height=8)
print(p)
dev.off()



HOG_all_df <- Medtru_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(Medtru=="Up" & SignalP=="signal_peptide") %>% dplyr::summarise(Overlap=length(Gene_id))
HOG_all_df$Trait <- "Medtru"
# trait_df <- Medtru_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(FIId=="Up") %>% dplyr::summarise(Overlap=length(Gene_id))
# trait_df$Trait <- "FIId"
trait_sigP_df <- Medtru_NODE_genes %>% group_by(Node_label) %>% dplyr::filter(FIId=="Up" & SignalP=="signal_peptide") %>% dplyr::summarise(Overlap=length(Gene_id))
trait_sigP_df$Trait <- "FIId_signalP"
# HOG_all_df <- rbind(HOG_all_df,trait_df)
HOG_all_df <- rbind(HOG_all_df,trait_sigP_df)

Full_summary <- HOG_all_df

Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
for(i in unique(Full_summary$Trait)){
  Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
  
}

Medtru_NODES <- unique(Full_summary$Node_label)
for(Node in Medtru_NODES){
  for(i in unique(Full_summary$Trait)){
    CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Node_label ==Node]
    ControlNode <- Full_summary$Overlap[Full_summary$Trait=="Medtru" & Full_summary$Node_label ==Node]
    M <- as.table(rbind(c(CloneNode, ControlNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Trait=="Medtru"]-ControlNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- test$p.value
    Full_summary$FishTest[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$Node_label ==Node & Full_summary$Trait==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$Node_label =="N1" & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
  }
}

Full_summary$Total[Full_summary$Trait=="Medtru" & Full_summary$Node_label =="N1"] <- Clones_sum$Sum[Clones_sum$Trait=="Medtru"]

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than Medtru"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than Medtru"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than Medtru"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than Medtru"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than Medtru"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than Medtru"

# End Medtru
# fwrite(Full_summary,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Medtru_specific_Medtru_HOG_Zones_",UP_DOWN,".txt"),sep="\t")
Full_summary_Up <- rbind(Full_summary,Mimpud_Full_summary)

Full_summary_Up$Node_label[Full_summary_Up$Node_label=="Medtru"] <- "Species specific shared"
Full_summary_Up$Node_label[Full_summary_Up$Node_label=="Mimpud"] <- "Species specific shared"
Full_summary_Up$Node_label <- factor(Full_summary_Up$Node_label,levels = c("N1","N2","N7","N12","N21","Species specific shared","Species_specific"))
Full_summary_Up$Trait <- factor(Full_summary_Up$Trait,
                             levels = c("Medtru","FIId_signalP","Mimpud","Release_signalP","Persist_signalP"))

mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N7"="darkseagreen2",
           "N12" = "darkseagreen3","N21" = "darkseagreen4", "Species specific shared" = "#E1AF00","Species_specific" = "grey70")

p <-ggplot(data=Full_summary_Up, aes(x=Trait, y=Overlap, fill=Node_label)) +
        geom_bar(stat="identity", position="fill",  width = 0.65) +
        scale_fill_manual(values = mycol,labels = c("NFN clade", "Fabales clade", "Papilionoideae clade",
                                                    "Dalbergioid+Hologalegina clade","Hologalegina clade", "Species specific shared","Species_specific")) + 
        #scale_fill_manual(values = c("#3B9AB2","#78B7C5","darkseagreen2","darkseagreen3", "darkseagreen4","#E1AF00")) +
        labs(y= "Number of genes shared", x = "Traits") +
        geom_text(data = Full_summary_Up, aes(y = 1.05, label = Total)) +
        geom_text(data = Full_summary_Up, aes(y = Overlap, label = Signif),size=8,position = position_fill(vjust = 0.45)) + theme_classic(base_size = 16)

pdf(file = paste0(out_path,"../trait_results/Barplot_Proportion_shared_with_Release_HOG_signalP_",UP_DOWN,".pdf"), width=12, height=8)
print(p)
dev.off()


}



