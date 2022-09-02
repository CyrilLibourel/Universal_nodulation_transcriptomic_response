# Script to count HOGs on species tree
library("data.table");library("ape");library("seqinr");library("adephylo");library("tidyverse");library("tidyr");library("VennDiagram")
library("UpSetR");library("gaston");library("bio3d");library("phytools");library("ggtree");library("phylobase");library("reshape2")
library("ggpubr");library("wesanderson");library("phangorn")
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

######## Multi species ########

UP_or_DOWN <- c("Up","Down")
species <- c("Aeseve", "Arahyp", "Datglo", "Hiprha", "Lotjap", "Lupalb", "Medtru","Mimpud","Parand")
treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
tree <- as(treeObj, "phylo4")

for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_",UP_DOWN,".tsv"))
  # Load DEGs Mimpud WT and fonctional traits
  Summary_count <- data.frame()
  for(i in species){
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
    # cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, paste0("^",i,"_")))
    # cor_Genes_HOG <- dplyr::left_join(cor_Genes_HOG,Recap_HOG,by="HOG")
    # cor_Genes_HOG <- cor_Genes_HOG  %>% dplyr::select(Gene_id,HOG,paste0("Node_label_",UP_DOWN))
    # names(cor_Genes_HOG) <- c("Gene_id","HOG","Node_label")
    # NODES <- c(names(ancestors(tree, i)),i)
    # cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)
    # cor_Genes_HOG <- cor_Genes_HOG %>% group_by(Node_label) %>% dplyr::summarise(Overlap=length(Gene_id))
    # cor_Genes_HOG$species <- i
    # cor_Genes_HOG$Node_label[cor_Genes_HOG$Node_label==i] <- "Species specific"
    # cor_Genes_HOG$Total[cor_Genes_HOG$Node_label=="N1"] <- sum(cor_Genes_HOG$Overlap)
    # Summary_count <- rbind(Summary_count,cor_Genes_HOG)
  }
  
  HOGs_tree <- keep.tip(treeObj,species)
  Summary_count$species <- factor(Summary_count$species,levels = HOGs_tree$tip.label)
  Summary_count$Node_label <- factor(Summary_count$Node_label,levels = c("N1","N2","N3","N7","N12","N15","N21","N22","Species specific"))
  mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N3"="brown1", "N15"="coral1","N7"="darkseagreen2",
             "N12" = "darkseagreen3","N21" = "darkseagreen4","N22"="aquamarine3", "Species specific" = "#E1AF00")

  p <- ggplot(data=Summary_count, aes(x=species, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Rosales clade","Papilionoideae clade",
                                                "Dalbergioid+Hologalegina clade","Hologalegina clade","Dalbergioid clade", "Species specific")) + 
    labs(y= "Proportion of genes shared with the different nodes", x = "Traits") +
    geom_text(data = Summary_count, aes(y = 1.05, label = Total)) +
    theme_classic(base_size = 18)
    
    pdf(file = paste0(out_path,"../results/Barplot_Proportion_shared_with_multi_Nodes_all_species_",UP_DOWN,".pdf"), width=12, height=8)
    print(p)
    dev.off()
    
    p <- ggplot(data=Summary_count, aes(x=species, y=Overlap, fill=Node_label)) +
      geom_bar(stat="identity", position="stack",  width = 0.65) +
      scale_fill_manual(values = mycol,labels = c("NFN clade", "Fabales clade", "Rosales+Cucurbitales clade","Rosales clade","Papilionoideae clade",
                                                  "Dalbergioid+Hologalegina clade","Hologalegina clade","Dalbergioid clade", "Species specific")) + 
      labs(y= "Number of genes shared with the different nodes", x = "Traits") +
      geom_text(data = Summary_count, aes(y =  max(Total*1.05,na.rm=TRUE), label = Total)) +
      theme_classic(base_size = 18)
    
    pdf(file = paste0(out_path,"../results/Barplot_count_shared_with_multi_Nodes_all_species_",UP_DOWN,".pdf"), width=12, height=8)
    print(p)
    dev.off()
}


