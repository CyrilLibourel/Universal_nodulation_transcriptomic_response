######################
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","phangorn"), require, character.only = TRUE)

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")

UPDOWN <- c("Up","Down")

  data = fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking.tsv"))
  data[data == "Up"] <- 1;data[data == "Down"] <- 1;data[data == "ns"] <- 0
  fwrite(data,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries.txt"),sep="\t",quote=FALSE)
  data = fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries.txt"))
  All_data <- data %>% dplyr::select(HOG,Acronym,Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                                     Ranking_Down,Ranking_Up,Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
  # allNod HOG
  to_select <- c(grep("_Nod_",names(data), value=TRUE))
  SPECIES <- data.frame(class=to_select)
  SPECIES <- tidyr::separate(SPECIES,class,into=c("Species","Nod","FDR","logFC","UpDown"),sep = "_")
  SPECIES <- unique(SPECIES$Species)

  for(i in SPECIES){
    species_to_select <- grep(i,to_select, value=TRUE)
    data_select_yes <- data %>% dplyr::select(HOG,all_of(species_to_select))
    data_select <- data_select_yes %>% dplyr::select(HOG,ends_with("_Up"))
    data_select$sum <- data_select_yes %>% dplyr::select(ends_with("_Up")) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
    data_select <- data_select %>% dplyr::filter(sum>0) %>% dplyr::select(HOG)
    data_select <- data_select %>% mutate(HOG2 = 1)
    names(data_select) <- c("HOG",paste0(i,"_Nod_Up"))
    All_data <- left_join(All_data,data_select)
    
    data_select <- data_select_yes %>% dplyr::select(HOG,ends_with("_Down"))
    data_select$sum <- data_select_yes %>% dplyr::select(ends_with("_Down")) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
    data_select <- data_select %>% dplyr::filter(sum>0) %>% dplyr::select(HOG)
    data_select <- data_select %>% mutate(HOG2 = 1)
    names(data_select) <- c("HOG",paste0(i,"_Nod_Down"))
    All_data <- left_join(All_data,data_select)
  }
  
  All_data[is.na(All_data)] <- ""
  fwrite(All_data,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries_Nodulation_data.txt"),sep="\t",quote=FALSE)
  

HOG <- All_data
UPDOWN <- c("Up","Down")

  All_data <-  as.data.frame(fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries_Nodulation_data.txt"),sep="\t",quote=FALSE))
  All_data[is.na(All_data)] <- 0
  for(updown in UPDOWN){
    HOG <- All_data %>% dplyr::select(HOG,ends_with(paste0("Nod_",updown)))
    HOG$sum <- HOG %>% dplyr::select(ends_with(paste0("_",updown))) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
    HOG <- HOG %>% dplyr::filter(sum>0) %>% dplyr::select(-sum)
    names(HOG) <- c("HOG",gsub(paste0("_Nod_",updown),"",names(HOG)[-1]))
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
    HOGs_tree <- keep.tip(treeObj,names(HOG)[-1])
    HOG <- HOG %>% dplyr::select(HOG,all_of(HOGs_tree$tip.label))
    HOG <- HOG %>% tidyr::unite("PASTE",all_of(names(HOG)[-1]), sep="_", remove=FALSE)
    Node_HOG <- data.frame()
    pb <- txtProgressBar(min = 0, max = length(unique(HOG$PASTE)), style = 3)

    for(i in unique(HOG$PASTE)){
      subset_HOG <- dplyr::filter(HOG,PASTE==i)
      x = as.integer(as.character(subset_HOG[1,] %>% dplyr::select(all_of(HOGs_tree$tip.label))))+1
      names(x) <- HOGs_tree$tip.label
      # reconstruct node states via MPR
      results <- asr_max_parsimony(HOGs_tree, x, 2,weight_by_scenarios=TRUE)
      node_states <- data.frame(State=as.integer(as.character(max.col(results$ancestral_likelihoods-0.01, ties.method="first")))-1,
                                Node_label=HOGs_tree$node.label)
      
      To_delete <- c()
      for(NODE in unique(node_states$Node[node_states$State==1])){
        Nodes_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$node.label[-1]
        Species_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$tip.label
        To_delete <- unique(c(To_delete,Nodes_to_delete,Species_to_delete))
      }

      node_states <- rbind(node_states,data.frame(State=x-1,Node_label=HOGs_tree$tip.label))
      node_states <- node_states %>% filter(State==1) %>% filter(!Node_label %in% To_delete) %>% dplyr::select(Node_label)
      
      node_states$PASTE <- i
      subset_HOG <- dplyr::left_join(subset_HOG,node_states, by="PASTE")
      Node_HOG <- rbind(Node_HOG,subset_HOG)
      setTxtProgressBar(pb, match(i,unique(HOG$PASTE)))
    }
    Node_HOG <- dplyr::arrange(Node_HOG,HOG)
    Node_HOG <- dplyr::select(Node_HOG,-PASTE)
    
    fwrite(Node_HOG,paste0("./results/Recap_HOG_by_node_weightParsimony_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
    summary_Node_HOG <- Node_HOG  %>% group_by(Node_label) %>% summarise(Count = length(HOG))
    #summary_Node_HOG
    fwrite(summary_Node_HOG,paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
    
    # Plot phylogenetic tree
    summary_Node_HOG_corrected <- fread(paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"))
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    HOGs_tree <- keep.tip(treeObj,c("Aeseve","Arahyp","Datglo","Hiprha","Lotjap","Medtru","Mimpud","Parand","Lupalb"))
    
    Nodes_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$node.label),Node_label=HOGs_tree$node.label)
    Nodes_to_plot <- left_join(Nodes_to_plot,summary_Node_HOG_corrected)
    Nodes_to_plot$cex <- log(Nodes_to_plot$Count)/4
    Nodes_to_plot <- dplyr:::filter(Nodes_to_plot,Count!="NA")
    
    Species_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$tip.label),Node_label=HOGs_tree$tip.label)
    Species_to_plot <- left_join(Species_to_plot,summary_Node_HOG_corrected)
    
    edges <- as.data.frame(HOGs_tree$edge); names(edges) <- c("NodeID","ID")
    edges$edgeID <- c(1:nrow(edges));edges <- dplyr::select(edges,-NodeID)
    Edges_to_plot <- left_join(Species_to_plot,edges)
    Edges_to_plot$cex <- log(Edges_to_plot$Count)/4
    Edges_to_plot <- dplyr:::filter(Edges_to_plot,Count!="NA")
    
    pdf(file = paste0("./results/Count_HOG_by_node_weightParsimony_",updown,".pdf"),width=20, height=20)
    {plot.phylo(HOGs_tree,cex=2.5)
      for(i in Nodes_to_plot$ID){
        nodelabels(Nodes_to_plot$Count[Nodes_to_plot$ID==i], i, frame = "c", bg = "coral1", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
      }
      for(i in Edges_to_plot$edgeID){
        edgelabels(Edges_to_plot$Count[Edges_to_plot$edgeID==i], i, frame = "c", bg = "darkgoldenrod1", font = 2, cex=Edges_to_plot$cex[Edges_to_plot$edgeID==i])
      }}
    dev.off()
    
    # Full tree
    summary_Node_HOG_corrected <- fread(paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"))
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
    
    Nodes_to_plot <- data.frame(ID=nodeid(treeObj, treeObj$node.label),Node_label=treeObj$node.label)
    Nodes_to_plot <- left_join(Nodes_to_plot,summary_Node_HOG_corrected)
    Nodes_to_plot$cex <- log(Nodes_to_plot$Count)/4
    Nodes_to_plot <- dplyr:::filter(Nodes_to_plot,Count!="NA")
    
    Species_to_plot <- data.frame(ID=nodeid(treeObj, treeObj$tip.label),Node_label=treeObj$tip.label)
    Species_to_plot <- left_join(Species_to_plot,summary_Node_HOG_corrected)
    
    edges <- as.data.frame(treeObj$edge); names(edges) <- c("NodeID","ID")
    edges$edgeID <- c(1:nrow(edges));edges <- dplyr::select(edges,-NodeID)
    Edges_to_plot <- left_join(Species_to_plot,edges)
    Edges_to_plot$cex <- log(Edges_to_plot$Count)/4
    Edges_to_plot <- dplyr:::filter(Edges_to_plot,Count!="NA")
    
    pdf(file = paste0("./results/Count_HOG_by_node_weightParsimony_",updown,"_fullTree.pdf"),width=20, height=20)
    {plot(treeObj,cex=2.5)
      for(i in Nodes_to_plot$ID){
        nodelabels(Nodes_to_plot$Count[Nodes_to_plot$ID==i], i, frame = "c", bg = "coral1", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
      }
      for(i in Edges_to_plot$edgeID){
        edgelabels(Edges_to_plot$Count[Edges_to_plot$edgeID==i], i, frame = "c", bg = "darkgoldenrod1", font = 2, cex=Edges_to_plot$cex[Edges_to_plot$edgeID==i])
      }}
    dev.off()
    
    
  }


treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
#treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
Nodes_to_plot <- data.frame(ID=nodeid(treeObj, treeObj$node.label),Node_label=treeObj$node.label)
pdf(file = paste0("./results/Tree_with_node_number.pdf"),width=20, height=20)
{plot(treeObj,cex=2.5)
  for(i in Nodes_to_plot$ID){
    nodelabels(Nodes_to_plot$Node_label[Nodes_to_plot$ID==i], i, frame = "c", bg = "white", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
  }
  }
dev.off()


# Add Nodes info to HOG
Resume_HOGs <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking.tsv"))

Nodes_Up <- fread(paste0("./results/Recap_HOG_by_node_weightParsimony_Up.txt"))
Nodes_Up <- dplyr::select(Nodes_Up,HOG,Node_label);names(Nodes_Up) <- c("HOG","Node_label_Up")
Resume_HOGs_Up <- dplyr::select(Resume_HOGs,HOG,Acronym, Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                                Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
Resume_HOGs_Up <- left_join(Resume_HOGs_Up,Nodes_Up)
Resume_HOGs_Up2 <- dplyr::select(Resume_HOGs,HOG,ends_with("_Up"))
Resume_HOGs_Up <- left_join(Resume_HOGs_Up,Resume_HOGs_Up2)
fwrite(Resume_HOGs_Up,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_Up.tsv"),sep="\t",quote=FALSE)


Nodes_Down <- fread(paste0("./results/Recap_HOG_by_node_weightParsimony_Down.txt"))
Nodes_Down <- dplyr::select(Nodes_Down,HOG,Node_label);names(Nodes_Down) <- c("HOG","Node_label_Down")
Resume_HOGs_Down <- dplyr::select(Resume_HOGs,HOG,Acronym, Category,Description,Mean_seq_number_HOG,Mean_seq_number_NodspecieswithDEGs,
                                  Percentage_of_Nodspecies_in_HOG,Percentage_of_nonNodspecies_in_HOG)
Resume_HOGs_Down <- left_join(Resume_HOGs_Down,Nodes_Down)
Resume_HOGs_Down2 <- dplyr::select(Resume_HOGs,HOG,ends_with("_Down"))
Resume_HOGs_Down <- left_join(Resume_HOGs_Down,Resume_HOGs_Down2)
fwrite(Resume_HOGs_Down,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_Down.tsv"),sep="\t",quote=FALSE)

######## Possible nodes for each HOG based on presence/absence ##############
HOG <- as.data.frame(fread(paste0("../Results_Apr14/Phylogenetic_Hierarchical_Orthogroups/N0_corrected.tsv"),h=T))
HOG_binary <- HOG %>% dplyr::select(all_of(SPECIES))
HOG_binary[HOG_binary!=""] <- 1; HOG_binary[HOG_binary==""] <- 0
HOG_binary <- HOG_binary %>% mutate_if(is.character, as.numeric)
HOG_binary$HOG <- HOG$HOG

HOG <- HOG_binary

  HOG$sum <- HOG %>% dplyr::select(-HOG) %>% dplyr::mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
  HOG <- HOG %>% dplyr::filter(sum>0) %>% dplyr::select(-sum)
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
  HOGs_tree <- keep.tip(treeObj,names(HOG %>% dplyr::select(all_of(SPECIES))))
  HOG <- HOG %>% dplyr::select(HOG,all_of(HOGs_tree$tip.label))
  HOG <- HOG %>% tidyr::unite("PASTE",all_of(names(HOG)[-1]), sep="_", remove=FALSE)
  Node_HOG <- data.frame()
  pb <- txtProgressBar(min = 0, max = length(unique(HOG$PASTE)), style = 3)
  
  for(i in unique(HOG$PASTE)){
    subset_HOG <- dplyr::filter(HOG,PASTE==i)
    x = as.integer(as.character(subset_HOG[1,] %>% dplyr::select(all_of(HOGs_tree$tip.label))))+1
    names(x) <- HOGs_tree$tip.label
    # reconstruct node states via MPR
    results <- asr_max_parsimony(HOGs_tree, x, 2,weight_by_scenarios=TRUE)
    node_states <- data.frame(State=as.integer(as.character(max.col(results$ancestral_likelihoods-0.01, ties.method="first")))-1,
                              Node_label=HOGs_tree$node.label)
    
    To_delete <- c()
    for(NODE in unique(node_states$Node[node_states$State==1])){
      Nodes_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$node.label[-1]
      Species_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$tip.label
      To_delete <- unique(c(To_delete,Nodes_to_delete,Species_to_delete))
    }
    
    node_states <- rbind(node_states,data.frame(State=x-1,Node_label=HOGs_tree$tip.label))
    node_states <- node_states %>% filter(State==1) %>% filter(!Node_label %in% To_delete) %>% dplyr::select(Node_label)
    
    node_states$PASTE <- i
    subset_HOG <- dplyr::left_join(subset_HOG,node_states, by="PASTE")
    Node_HOG <- rbind(Node_HOG,subset_HOG)
    setTxtProgressBar(pb, match(i,unique(HOG$PASTE)))
  }
  Node_HOG <- dplyr::arrange(Node_HOG,HOG)
  Node_HOG <- dplyr::select(Node_HOG,-PASTE)
  
  fwrite(Node_HOG,paste0("./results/Recap_HOG_node_weightParsimony.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  summary_Node_HOG <- Node_HOG  %>% group_by(Node_label) %>% summarise(Count = length(HOG))
  #summary_Node_HOG
  fwrite(summary_Node_HOG,paste0("./results/Summary_count_HOG_node_weightParsimony.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  
  # Plot phylogenetic tree
  summary_Node_HOG_corrected <- fread(paste0("./results/Summary_count_HOG_node_weightParsimony.txt"))
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
  HOGs_tree <- keep.tip(treeObj,names(HOG %>% dplyr::select(all_of(SPECIES))))
  Nodes_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$node.label),Node_label=HOGs_tree$node.label)
  Nodes_to_plot <- left_join(Nodes_to_plot,summary_Node_HOG_corrected)
  Nodes_to_plot$cex <- log(Nodes_to_plot$Count)/4
  Nodes_to_plot <- dplyr:::filter(Nodes_to_plot,Count!="NA")
  
  Species_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$tip.label),Node_label=HOGs_tree$tip.label)
  Species_to_plot <- left_join(Species_to_plot,summary_Node_HOG_corrected)
  
  edges <- as.data.frame(HOGs_tree$edge); names(edges) <- c("NodeID","ID")
  edges$edgeID <- c(1:nrow(edges));edges <- dplyr::select(edges,-NodeID)
  Edges_to_plot <- left_join(Species_to_plot,edges)
  Edges_to_plot$cex <- log(Edges_to_plot$Count)/4
  Edges_to_plot <- dplyr:::filter(Edges_to_plot,Count!="NA")
  
  pdf(file = paste0("./results/Count_HOG_node_weightParsimony.pdf"),width=20, height=20)
  {plot.phylo(HOGs_tree,cex=2.5)
    for(i in Nodes_to_plot$ID){
      nodelabels(Nodes_to_plot$Count[Nodes_to_plot$ID==i], i, frame = "c", bg = "coral1", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
    }
    for(i in Edges_to_plot$edgeID){
      edgelabels(Edges_to_plot$Count[Edges_to_plot$edgeID==i], i, frame = "c", bg = "darkgoldenrod1", font = 2, cex=Edges_to_plot$cex[Edges_to_plot$edgeID==i])
    }}
  dev.off()
  


#### random gene_id based on possible genes in nodes ####
UPDOWN <- c("Up","Down")

All_data <-  as.data.frame(fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries_Nodulation_data.txt"),sep="\t",quote=FALSE))
All_data[is.na(All_data)] <- 0

HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

All_HOGs <- data.table(HOG=unique(HOG$HOG))
UPDOWN <- c("Down","Up")

for(updown in UPDOWN){
  Recap_HOG <- fread(paste0("./results/Recap_HOG_by_node_weightParsimony_",updown,".txt"))
  Recap_HOG <- Recap_HOG %>% dplyr:::select(HOG,Node_label)
  subtables <- grep(updown,tables, value=TRUE)
  all_summary_Node_HOG <- data.frame()
  species_node_summary <- data.frame(Node_label=unique(Recap_HOG$Node_label))
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    names(DEGs) <- "Gene_id"
    temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
    temp <- dplyr::left_join(temp,Recap_HOG,by="HOG")
    temp[is.na(temp)] <- ""
    temp <- temp %>% dplyr::filter(HOG!="")
    sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    tree <- as(treeObj, "phylo4")
    NODES <- c(names(ancestors(tree, sp)),sp)
    temp <- dplyr::filter(temp,Node_label %in% NODES)
    summary_cor_Genes_HOG <- temp %>% group_by(Node_label) %>% summarise(Count = length(Gene_id))
    names(summary_cor_Genes_HOG) <- c("Node_label",sp)
    
    species_node_summary <- dplyr::left_join(species_node_summary,summary_cor_Genes_HOG,by="Node_label")
  }
  
  species_node_summary[is.na(species_node_summary)] <- 0
  fwrite(species_node_summary,paste0("./results/Table_count_genes_by_node_species_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  
  
  recap_Node_HOG_corrected <- fread(paste0("./results/Recap_HOG_node_weightParsimony.txt"))
  recap_Node_HOG_corrected <- recap_Node_HOG_corrected %>% dplyr:::select(HOG,Node_label)
  
  subtables <- grep(updown,tables, value=TRUE)
  all_summary_Node_HOG <- data.frame()
  species_node_summary <- data.frame(Node_label=unique(Recap_HOG$Node_label))
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    names(DEGs) <- "Gene_id"
    temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
    temp <- dplyr::left_join(temp,recap_Node_HOG_corrected,by="HOG")
    temp[is.na(temp)] <- ""
    temp <- temp %>% dplyr::filter(HOG!="")
    sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    tree <- as(treeObj, "phylo4")
    NODES <- c(names(ancestors(tree, sp)),sp)
    temp <- dplyr::filter(temp,Node_label %in% NODES)
    summary_cor_Genes_HOG <- temp %>% group_by(Node_label) %>% summarise(Count = length(Gene_id))
    names(summary_cor_Genes_HOG) <- c("Node_label",sp)
    
    species_node_summary <- dplyr::left_join(species_node_summary,summary_cor_Genes_HOG,by="Node_label")
  }

  species_node_summary[is.na(species_node_summary)] <- 0
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  recap_Node_HOG_corrected <- fread(paste0("./results/Recap_HOG_node_weightParsimony.txt"))
  
for(rep in c(1:1000)){
  All_ranking <- HOG %>% dplyr::select(HOG)

for(i in subtables){
  # DEGs <- fread(paste0(i,".txt"),h=FALSE)
  # names(DEGs) <- "Gene_id"
  # temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
  # temp[is.na(temp)] <- ""
  # temp <- temp %>% dplyr::filter(HOG!="")
  sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
  NODES <- c(names(ancestors(tree, sp)),sp)
  cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
  cor_Genes_HOG <- dplyr:::left_join(cor_Genes_HOG,recap_Node_HOG_corrected,by="HOG")
  cor_Genes_HOG <- dplyr::filter(cor_Genes_HOG,Node_label %in% NODES)

  temp_species_node_summary <- species_node_summary %>% dplyr::select(Node_label,sp) %>% dplyr::filter_at(2,all_vars(. > 0))
  DEGs <- data.frame()
  
  for(sampling in temp_species_node_summary$Node_label){
    size <- as.numeric(temp_species_node_summary %>% dplyr::filter_at(1,all_vars(. %in% sampling)) %>% dplyr::select(sp))
    temp_DEGs <- data.frame(Gene_id=sample(cor_Genes_HOG$Gene_id[cor_Genes_HOG$Node_label==sampling],size))
    DEGs <- rbind(DEGs,temp_DEGs)
  }

  temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
  names(temp) <- c(sp,"HOG")
  temp[is.na(temp)] <- ""
  temp[temp[,1] != "",1] <- 1
  temp <- temp[!duplicated(temp[,c('HOG')]),]
  All_ranking <- dplyr::left_join(All_ranking,temp,by="HOG")
  All_ranking[is.na(All_ranking)] <- ""
  # message(sp," Done")
}
  
  All_ranking <- All_ranking[!duplicated(All_ranking[,c('HOG')]),]
  fwrite(All_ranking,paste0("./results_sampling_by_Node/Sampling_test_Nodulation_data.txt"),sep="\t",quote=FALSE)
  All_ranking <- fread("./results_sampling_by_Node/Sampling_test_Nodulation_data.txt",h=T)
  All_ranking[is.na(All_ranking)] <- 0

  All_ranking$sum <- All_ranking %>% dplyr::select(-HOG) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
  All_ranking <- All_ranking %>% dplyr::filter(sum>0) %>% dplyr::select(-sum)

    #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
    HOGs_tree <- keep.tip(treeObj,names(All_ranking)[-1])
    All_ranking <- All_ranking %>% dplyr::select(HOG,all_of(HOGs_tree$tip.label))
    All_ranking <- All_ranking %>% tidyr::unite("PASTE",all_of(names(All_ranking)[-1]), sep="_", remove=FALSE)
    Node_HOG <- data.frame()
    # message("Node parsimony construction")
    # pb <- txtProgressBar(min = 0, max = length(unique(All_ranking$PASTE)), style = 3)

    for(i in unique(All_ranking$PASTE)){
      subset_HOG <- dplyr::filter(All_ranking,PASTE==i)
      x = as.integer(as.character(subset_HOG[1,] %>% dplyr::select(all_of(HOGs_tree$tip.label))))+1
      names(x) <- HOGs_tree$tip.label
      # reconstruct node states via MPR
      results <- asr_max_parsimony(HOGs_tree, x, 2,weight_by_scenarios=TRUE)
      node_states <- data.frame(State=as.integer(as.character(max.col(results$ancestral_likelihoods-0.01, ties.method="first")))-1,
                                Node_label=HOGs_tree$node.label)
      
      To_delete <- c()
      for(NODE in unique(node_states$Node[node_states$State==1])){
        Nodes_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$node.label[-1]
        Species_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$tip.label
        To_delete <- unique(c(To_delete,Nodes_to_delete,Species_to_delete))
      }
      
      node_states <- rbind(node_states,data.frame(State=x-1,Node_label=HOGs_tree$tip.label))
      node_states <- node_states %>% filter(State==1) %>% filter(!Node_label %in% To_delete) %>% dplyr::select(Node_label)
      
      node_states$PASTE <- i
      subset_HOG <- dplyr::left_join(subset_HOG,node_states, by="PASTE")
      Node_HOG <- rbind(Node_HOG,subset_HOG)
      # setTxtProgressBar(pb, match(i,unique(All_ranking$PASTE)))
    }
    message("")
    Node_HOG <- dplyr::arrange(Node_HOG,HOG)
    Node_HOG <- dplyr::select(Node_HOG,HOG,Node_label)

    fwrite(Node_HOG,paste0("./results_sampling_by_Node/Sampling/Recap_HOG_sampling_",rep,"_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
    summary_Node_HOG <- Node_HOG  %>% group_by(Node_label) %>% summarise(Count = length(HOG))
    summary_Node_HOG$Sampling <- rep
    #summary_Node_HOG
    all_summary_Node_HOG <- rbind(all_summary_Node_HOG,summary_Node_HOG)
message("Sampling ",rep," Done")
message("")
}

fwrite(all_summary_Node_HOG,paste0("./results_sampling_by_Node/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")

}


for(updown in UPDOWN){
  Real_summary_Node_HOG <- fread(paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"))
  all_summary_Node_HOG <- fread(paste0("./results_sampling_by_Node/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"),h=T)
  summary_Node_HOG <- all_summary_Node_HOG  %>% group_by(Node_label) %>% summarise(Q001=round(quantile(Count,0.001,na.rm=TRUE),0),
                                                                                   Q01=round(quantile(Count,0.01,na.rm=TRUE),0),
                                                                                   Q05=round(quantile(Count,probs=0.05,na.rm=TRUE),0),
                                                                                   Q95=round(quantile(Count,probs=0.95,na.rm=TRUE),0),
                                                                                   Q99=round(quantile(Count,probs=0.99,na.rm=TRUE),0),
                                                                                   Q999=round(quantile(Count,probs=0.999,na.rm=TRUE),0))

  summary_mean_Node_HOG <- all_summary_Node_HOG  %>% group_by(Node_label) %>% summarise(MeanCount = round(mean(Count,na.rm=TRUE),0),
                                                                               MedianCount = round(median(Count,na.rm=TRUE),0))
  summary_Node_HOG <- dplyr::left_join(summary_Node_HOG,summary_mean_Node_HOG)
  Real_summary_Node_HOG <- dplyr::left_join(Real_summary_Node_HOG,summary_Node_HOG)
  Real_summary_Node_HOG$Enrich <- Real_summary_Node_HOG$Count/Real_summary_Node_HOG$MeanCount
  fwrite(Real_summary_Node_HOG,paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,"_quantile_sampling_by_nodes.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
  HOGs_tree <- keep.tip(treeObj,c("Aeseve","Arahyp","Datglo","Hiprha","Lotjap","Medtru","Mimpud","Parand","Lupalb"))
  
  Nodes_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$node.label),Node_label=HOGs_tree$node.label)
  Nodes_to_plot <- left_join(Nodes_to_plot,summary_Node_HOG)
  Nodes_to_plot$cex <- log(Nodes_to_plot$MeanCount)/4
  Nodes_to_plot <- dplyr:::filter(Nodes_to_plot,MeanCount!="NA")
  
  Species_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$tip.label),Node_label=HOGs_tree$tip.label)
  Species_to_plot <- left_join(Species_to_plot,summary_Node_HOG)
  
  edges <- as.data.frame(HOGs_tree$edge); names(edges) <- c("NodeID","ID")
  edges$edgeID <- c(1:nrow(edges));edges <- dplyr::select(edges,-NodeID)
  Edges_to_plot <- left_join(Species_to_plot,edges)
  Edges_to_plot$cex <- log(Edges_to_plot$MeanCount)/4
  Edges_to_plot <- dplyr:::filter(Edges_to_plot,MeanCount!="NA")
  
  pdf(file = paste0("./results_sampling_by_Node/Count_sampling_HOG_by_node_weightParsimony_",updown,".pdf"),width=20, height=20)
  {plot.phylo(HOGs_tree,cex=2.5)
    for(i in Nodes_to_plot$ID){
      nodelabels(Nodes_to_plot$MeanCount[Nodes_to_plot$ID==i], i, frame = "c", bg = "coral1", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
    }
    for(i in Edges_to_plot$edgeID){
      edgelabels(Edges_to_plot$MeanCount[Edges_to_plot$edgeID==i], i, frame = "c", bg = "darkgoldenrod1", font = 2, cex=Edges_to_plot$cex[Edges_to_plot$edgeID==i])
    }}
  dev.off()
  
}

#### random gene_id ####
HOG <- All_data
UPDOWN <- c("Up","Down")

Recap_HOG <- fread(paste0("./results/Recap_HOG_node_weightParsimony.txt"))
Recap_HOG <- Recap_HOG %>% dplyr:::select(HOG,Node_label)

All_data <-  as.data.frame(fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries_Nodulation_data.txt"),sep="\t",quote=FALSE))
All_data[is.na(All_data)] <- 0

HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

All_HOGs <- data.table(HOG=unique(HOG$HOG))
UPDOWN <- c("Down","Up")

for(updown in UPDOWN){
  subtables <- grep(updown,tables, value=TRUE)
  all_summary_Node_HOG <- data.frame()
  
  for(i in subtables){
    DEGs <- fread(paste0(i,".txt"),h=FALSE)
    names(DEGs) <- "Gene_id"
    temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
    temp <- dplyr::left_join(temp,Recap_HOG,by="HOG")
    temp[is.na(temp)] <- ""
    temp <- temp %>% dplyr::filter(HOG!="")
    sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    treeObj <- as(treeObj, "phylo4")
    NODES <- c(names(ancestors(treeObj, sp)),sp)
    cor_Genes_HOG <- dplyr::filter(temp,Node_label %in% NODES)
    summary_cor_Genes_HOG <- cor_Genes_HOG  %>% group_by(Node_label) %>% summarise(Count = length(Gene_id))
    
    cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
    DEGs <- data.frame(Gene_id=sample(cor_Genes_HOG$Gene_id,nrow(temp)))
    
    temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
    names(temp) <- c(sp,"HOG")
    temp[is.na(temp)] <- ""
    temp[temp[,1] != "",1] <- 1
    temp <- temp[!duplicated(temp[,c('HOG')]),]
    All_ranking <- dplyr::left_join(All_ranking,temp,by="HOG")
    All_ranking[is.na(All_ranking)] <- ""
    message(sp," Done")
  } 
  
  
  for(rep in c(1:1000)){
    All_ranking <- HOG %>% dplyr::select(HOG)
    
    for(i in subtables){
      DEGs <- fread(paste0(i,".txt"),h=FALSE)
      names(DEGs) <- "Gene_id"
      temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
      temp[is.na(temp)] <- ""
      temp <- temp %>% dplyr::filter(HOG!="")
      sp <- unlist(strsplit(gsub("./rnaseq/","",i),split="_"))[1]
      cor_Genes_HOG <- HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
      DEGs <- data.frame(Gene_id=sample(cor_Genes_HOG$Gene_id,nrow(temp)))
      
      temp <- as.data.frame(dplyr::left_join(DEGs,HOG,by="Gene_id"))
      names(temp) <- c(sp,"HOG")
      temp[is.na(temp)] <- ""
      temp[temp[,1] != "",1] <- 1
      temp <- temp[!duplicated(temp[,c('HOG')]),]
      All_ranking <- dplyr::left_join(All_ranking,temp,by="HOG")
      All_ranking[is.na(All_ranking)] <- ""
      message(sp," Done")
    }
    
    All_ranking <- All_ranking[!duplicated(All_ranking[,c('HOG')]),]
    fwrite(All_ranking,paste0("./results_sampling/Sampling_test_Nodulation_data.txt"),sep="\t",quote=FALSE)
    All_ranking <- fread("./results_sampling/Sampling_test_Nodulation_data.txt",h=T)
    All_ranking[is.na(All_ranking)] <- 0
    
    All_ranking$sum <- All_ranking %>% dplyr::select(-HOG) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
    All_ranking <- All_ranking %>% dplyr::filter(sum>0) %>% dplyr::select(-sum)
    treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
    #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
    HOGs_tree <- keep.tip(treeObj,names(All_ranking)[-1])
    All_ranking <- All_ranking %>% dplyr::select(HOG,all_of(HOGs_tree$tip.label))
    All_ranking <- All_ranking %>% tidyr::unite("PASTE",all_of(names(All_ranking)[-1]), sep="_", remove=FALSE)
    Node_HOG <- data.frame()
    message("Node parsimony construction")
    pb <- txtProgressBar(min = 0, max = length(unique(All_ranking$PASTE)), style = 3)
    
    for(i in unique(All_ranking$PASTE)){
      subset_HOG <- dplyr::filter(All_ranking,PASTE==i)
      x = as.integer(as.character(subset_HOG[1,] %>% dplyr::select(all_of(HOGs_tree$tip.label))))+1
      names(x) <- HOGs_tree$tip.label
      # reconstruct node states via MPR
      results <- asr_max_parsimony(HOGs_tree, x, 2,weight_by_scenarios=TRUE)
      node_states <- data.frame(State=as.integer(as.character(max.col(results$ancestral_likelihoods-0.01, ties.method="first")))-1,
                                Node_label=HOGs_tree$node.label)
      
      To_delete <- c()
      for(NODE in unique(node_states$Node[node_states$State==1])){
        Nodes_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$node.label[-1]
        Species_to_delete <- extract.clade(HOGs_tree,nodeid(HOGs_tree, HOGs_tree$node.label[HOGs_tree$node.label==NODE]))$tip.label
        To_delete <- unique(c(To_delete,Nodes_to_delete,Species_to_delete))
      }
      
      node_states <- rbind(node_states,data.frame(State=x-1,Node_label=HOGs_tree$tip.label))
      node_states <- node_states %>% filter(State==1) %>% filter(!Node_label %in% To_delete) %>% dplyr::select(Node_label)
      
      node_states$PASTE <- i
      subset_HOG <- dplyr::left_join(subset_HOG,node_states, by="PASTE")
      Node_HOG <- rbind(Node_HOG,subset_HOG)
      setTxtProgressBar(pb, match(i,unique(All_ranking$PASTE)))
    }
    message("")
    Node_HOG <- dplyr::arrange(Node_HOG,HOG)
    Node_HOG <- dplyr::select(Node_HOG,HOG,Node_label)
    
    fwrite(Node_HOG,paste0("./results_sampling/Recap_HOG_sampling_",rep,"_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
    summary_Node_HOG <- Node_HOG  %>% group_by(Node_label) %>% summarise(Count = length(HOG))
    summary_Node_HOG$Sampling <- rep
    #summary_Node_HOG
    all_summary_Node_HOG <- rbind(all_summary_Node_HOG,summary_Node_HOG)
    message("Sampling ",rep," Done")
    message("")
  }
  
  fwrite(all_summary_Node_HOG,paste0("./results_sampling/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  
}


for(updown in UPDOWN){
  Real_summary_Node_HOG <- fread(paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"))
  all_summary_Node_HOG <- fread(paste0("./results_sampling/Summary_count_HOG_by_node_weightParsimony_",updown,".txt"),h=T)
  summary_Node_HOG <- all_summary_Node_HOG  %>% group_by(Node_label) %>% summarise(Q001=round(quantile(Count,0.001,na.rm=TRUE),0),
                                                                                   Q01=round(quantile(Count,0.01,na.rm=TRUE),0),
                                                                                   Q05=round(quantile(Count,probs=0.05,na.rm=TRUE),0),
                                                                                   Q95=round(quantile(Count,probs=0.95,na.rm=TRUE),0),
                                                                                   Q99=round(quantile(Count,probs=0.99,na.rm=TRUE),0),
                                                                                   Q999=round(quantile(Count,probs=0.999,na.rm=TRUE),0))
  
  summary_mean_Node_HOG <- all_summary_Node_HOG  %>% group_by(Node_label) %>% summarise(MeanCount = round(mean(Count,na.rm=TRUE),0),
                                                                                        MedianCount = round(median(Count,na.rm=TRUE),0))
  summary_Node_HOG <- dplyr::left_join(summary_Node_HOG,summary_mean_Node_HOG)
  Real_summary_Node_HOG <- dplyr::left_join(Real_summary_Node_HOG,summary_Node_HOG)
  Real_summary_Node_HOG$Enrich <- Real_summary_Node_HOG$Count/Real_summary_Node_HOG$MeanCount
  fwrite(Real_summary_Node_HOG,paste0("./results/Summary_count_HOG_by_node_weightParsimony_",updown,"_quantile_sampling.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
  treeObj <- read.tree(paste0("../Results_Apr14/Species_Tree/SpeciesTree_rooted_node_labels_Lotjap.txt"))
  #treeObj$node.label <- paste0("N",1:length(treeObj$node.label))
  HOGs_tree <- keep.tip(treeObj,c("Aeseve","Arahyp","Datglo","Hiprha","Lotjap","Medtru","Mimpud","Parand","Lupalb"))
  
  Nodes_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$node.label),Node_label=HOGs_tree$node.label)
  Nodes_to_plot <- left_join(Nodes_to_plot,summary_Node_HOG)
  Nodes_to_plot$cex <- log(Nodes_to_plot$MeanCount)/4
  Nodes_to_plot <- dplyr:::filter(Nodes_to_plot,MeanCount!="NA")
  
  Species_to_plot <- data.frame(ID=nodeid(HOGs_tree, HOGs_tree$tip.label),Node_label=HOGs_tree$tip.label)
  Species_to_plot <- left_join(Species_to_plot,summary_Node_HOG)
  
  edges <- as.data.frame(HOGs_tree$edge); names(edges) <- c("NodeID","ID")
  edges$edgeID <- c(1:nrow(edges));edges <- dplyr::select(edges,-NodeID)
  Edges_to_plot <- left_join(Species_to_plot,edges)
  Edges_to_plot$cex <- log(Edges_to_plot$MeanCount)/4
  Edges_to_plot <- dplyr:::filter(Edges_to_plot,MeanCount!="NA")
  
  pdf(file = paste0("./results_sampling/Count_sampling_HOG_by_node_weightParsimony_",updown,".pdf"),width=20, height=20)
  {plot.phylo(HOGs_tree,cex=2.5)
    for(i in Nodes_to_plot$ID){
      nodelabels(Nodes_to_plot$MeanCount[Nodes_to_plot$ID==i], i, frame = "c", bg = "coral1", font = 2, cex=Nodes_to_plot$cex[Nodes_to_plot$ID==i])
    }
    for(i in Edges_to_plot$edgeID){
      edgelabels(Edges_to_plot$MeanCount[Edges_to_plot$edgeID==i], i, frame = "c", bg = "darkgoldenrod1", font = 2, cex=Edges_to_plot$cex[Edges_to_plot$edgeID==i])
    }}
  dev.off()
  
}

#### Add random HOG count ####
All_data <-  as.data.frame(fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_binaries_Nodulation_data.txt"),sep="\t",quote=FALSE))
All_data[is.na(All_data)] <- 0

HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

tables <- paste0(list.files(path="./rnaseq/",pattern=".txt$",full.names=TRUE))
tables <- gsub(".txt","",tables)

All_HOGs <- data.table(HOG=unique(HOG$HOG))
UPDOWN <- c("Down","Up")

for(updown in UPDOWN){
  HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_",updown,".tsv"))
  Recap_HOG <- HOG %>% dplyr::select(HOG,paste0("Node_label_",updown))
  names(Recap_HOG) <- c("HOG","Node_label")
  Recap_HOG <- Recap_HOG %>% dplyr::filter(Node_label!="")
  Recap_HOG <- Recap_HOG %>% tidyr::unite("PASTE",all_of(names(Recap_HOG)), sep="_", remove=FALSE)
  
  names(HOG)[names(HOG) == paste0("Node_label_",updown)] <- "Node_label"
  HOG <- HOG %>% tidyr::unite("PASTE",HOG,Node_label, sep="_", remove=FALSE)

  # fileNames <- paste0(list.files(path="./results_sampling_by_Node/Sampling/",pattern=paste0(updown,".txt$"),full.names=TRUE))

  for (i in c(1:1000)){
    sampling <- fread(paste0("./results_sampling_by_Node/Sampling/Recap_HOG_sampling_",i,"_",updown,".txt"))
    sampling <- sampling %>% tidyr::unite("PASTE",all_of(names(sampling)), sep="_", remove=TRUE)
    sampling$sample <- 1; names(sampling) <- c("PASTE",i)
    Recap_HOG <- dplyr::left_join(Recap_HOG,sampling,by="PASTE")
  }

  Recap_HOG[is.na(Recap_HOG)] <- 0
  Recap_HOG$sum <- Recap_HOG %>% dplyr::select(c(4:length(Recap_HOG))) %>% mutate(sum = rowSums(., na.rm = TRUE)) %>% dplyr::select(sum)
  Recap_HOG <- dplyr::select(Recap_HOG,PASTE,sum)

  HOG <- dplyr::left_join(HOG,Recap_HOG,by="PASTE")
  HOG <- dplyr::select(HOG,-PASTE)
  names(HOG)[names(HOG) == "sum"] <- "Sampling_percentage"
  HOG$Sampling_percentage = HOG$Sampling_percentage/1000*100
  names(HOG)[names(HOG) == "Node_label"] <- paste0("Node_label_",updown)

  fwrite(HOG,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",updown,".tsv"),
                quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
}


