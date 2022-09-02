# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","reshape2",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ggpubr","wesanderson","phangorn"), require, character.only = TRUE)

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

######## cinetics ########
UP_or_DOWN <- c("Up","Down")
species <- c("Aeseve", "Arahyp", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")
treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

tables <- paste0(list.files(path="./rnaseq_timepoints/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables);tables <- gsub("./rnaseq_timepoints/","",tables)

for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZones",UP_DOWN,".tsv"))
  # Load DEGs Mimpud WT and fonctional traits

  for(SP in species){
  subtables <- grep(paste0(SP,".*",UP_DOWN),tables, value=TRUE)
  Summary_count <- data.frame()
  
  DEGs <- fread(paste0("./rnaseq/",list.files(path="./rnaseq/",pattern=paste0(SP,".*",UP_DOWN,".txt$"),full.names=FALSE)),h=FALSE)
  names(DEGs) <- "Gene_id"
  temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
  names(temp) <- c(SP,"HOG")
  temp[is.na(temp)] <- ""
  cor_Genes_HOG <- dplyr::left_join(temp,Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN)))
  
  names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
  NODES <- c(names(ancestors(tree, SP)),SP)
  cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
  cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
  cor_Genes_HOG$species <- SP
  cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==SP] <- "Species specific"
  cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
  Summary_count <- rbind(Summary_count,cor_Genes_HOG)

  for(i in subtables){
    DEGs <- fread(paste0("./rnaseq_timepoints/",i,".txt"),h=FALSE)
    names(DEGs) <- "Gene_id"
    temp <- as.data.frame(dplyr::left_join(DEGs,HOG))
    names(temp) <- c(i,"HOG")
    temp[is.na(temp)] <- ""
    cor_Genes_HOG <- dplyr::left_join(temp,Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN)))

    names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
    NODES <- c(names(ancestors(tree, SP)),SP)
    cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
    cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
    cor_Genes_HOG$species <- i
    cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==SP] <- "Species specific"
    cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
    Summary_count <- rbind(Summary_count,cor_Genes_HOG)
  }

  HOGs_tree <- keep.tip(treeObj,species)
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
    labs(y= "Number of genes shared", x = "Cinetics") +
    geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
    theme_classic(base_size = 18)

  pdf(file = paste0(out_path,"../results/Barplot_Proportion_shared_with_multi_Nodes_",SP,"_cinetics_",UP_DOWN,".pdf"), width=20, height=8)
  print(p)
  dev.off()

  }

}


