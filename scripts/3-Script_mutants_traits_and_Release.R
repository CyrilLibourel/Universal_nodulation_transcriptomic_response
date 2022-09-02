# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","reshape2","ComplexHeatmap",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ggpubr","wesanderson","DECIPHER","Biostrings"), require, character.only = TRUE)

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


######## Mimpud detailed traits ########
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_"))

NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")

for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",UP_DOWN,".tsv"))
  
  # Load DEGs Mimpud WT and fonctional traits
  CtaiWT_AllUp <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(CtaiWT_AllUp) <- "Gene_id"
  CtaiWT_AllUp <- left_join(CtaiWT_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_AllUp <- unique(CtaiWT_AllUp$Gene_id[!is.na(CtaiWT_AllUp$HOG)])
  
  CtaiWT_Other <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",setdiff(UP_or_DOWN,UP_DOWN),".txt"),header=FALSE)
  names(CtaiWT_Other) <- "Gene_id"
  CtaiWT_Other <- left_join(CtaiWT_Other,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_Other <- unique(CtaiWT_Other$Gene_id[!is.na(CtaiWT_Other$HOG)])
  
  HOG_id_NODE <- Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN),starts_with("Mimpud_"))
  names(HOG_id_NODE) <- c("HOG","Node_label","Mimpud")
  HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% NODES)
  Mimpud_NODE_genes <- left_join(HOG_id_NODE,Mimpud_cor_Genes_HOG,by="HOG")
  
  # Overlap CtaiWT_UP_DOWN and NODE
  Mimpud_NODE_genes_Up <- Mimpud_NODE_genes %>% dplyr::filter(Gene_id %in% CtaiWT_AllUp)
  Mimpud_NODE_genes_Up$Ctai <- 1

  TRAITS = c("Release","Persist","NFix")
  for(trait in TRAITS){
    Mimpud_trait <- Mimpud_NODE_genes_Up
    trait_df <- fread(input = paste0("./trait_analysis/",trait,".txt"),header=TRUE)
    names(trait_df)[1] <- "Gene_id"
    trait_df$Cond <- gsub("\\[","",trait_df$Cond)
    trait_df$Cond <- gsub("\\]","",trait_df$Cond)
    trait_df$Cond <- gsub("-","_vs_",trait_df$Cond)
    
    Conditions <- unique(trait_df$Cond)

    all_cond_df <- data.frame()
    for(cond in Conditions){
      cond_df <- trait_df %>% dplyr::filter(Cond==cond) %>% dplyr::select(Gene_id)
      cond_df$Trait <- 1;names(cond_df) <- c("Gene_id",cond)
      Mimpud_trait <- dplyr::left_join(Mimpud_trait,cond_df,by="Gene_id")
      Mimpud_trait <- Mimpud_trait %>% dplyr::distinct()
    }
    Mimpud_trait[is.na(Mimpud_trait)] <- 0
    
    HOG_all_df <- Mimpud_trait %>% distinct(HOG, .keep_all = TRUE) %>% dplyr::select(HOG,Ctai)
    for(cond in Conditions){
      trait_df <- eval(parse(text = paste0("Mimpud_trait %>% group_by(HOG) %>% dplyr::summarise(Overlap=max(",cond,"))")))
      # trait_df <- eval(parse(text = paste0("Mimpud_trait %>% group_by(HOG) %>% dplyr::filter(",trait,"==1) %>% dplyr::summarise(Overlap=max(",trait,"))")))
      names(trait_df) <- c("HOG",paste0(cond))
      
      HOG_all_df <- dplyr::left_join(HOG_all_df,trait_df,by="HOG")
    }
    HOG_all_df[HOG_all_df==1] <- UP_DOWN
    HOG_all_df[HOG_all_df==0] <- ""
    
    HOG_all_df <- left_join(Recap_HOG,HOG_all_df,by="HOG")
    HOG_all_df[is.na(HOG_all_df)] <- ""

    Nodes <- c("N1","N2","Mimpud")
    for(Node in Nodes){
      Node_df <- eval(parse(text = paste0("HOG_all_df %>% dplyr::filter(Node_label_",UP_DOWN,"==Node)")))
      Node_df[Node_df==UP_DOWN] <- 1
      Node_df[Node_df==""] <- 0
      Node_df <- Node_df %>% tidyr::unite("PASTE",all_of(paste0(Conditions)), sep="_", remove=FALSE)
      Node_df <- Node_df %>% dplyr::arrange(desc(PASTE))
      
      #####
      dat <- as.matrix(dplyr::select(Node_df,all_of(paste0(Conditions)))%>% mutate_if(is.character, as.numeric)) 
      rownames(dat) <- Node_df$HOG
      
      if(UP_DOWN=="Up"){my_palette <- c(colorRampPalette(colors=c("cornsilk2","firebrick"))(n=2))
      }else{my_palette <- c(colorRampPalette(colors=c("cornsilk2","royalblue"))(n=2))}
      
      ht1 <- Heatmap(dat, name = "HOG", cluster_rows = FALSE,
                     row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
                     cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
                     row_names_side = "left",row_order=rownames(dat), row_title = NULL,
                     row_split = factor(Node_df$PASTE, levels = unique(Node_df$PASTE)),
                     show_row_names = FALSE, col = my_palette , width = unit(25, "mm"), height = unit(14, "cm"),
                     # row_km = length(unique(Node_df$PASTE)),
                     row_title_rot = 0 ,heatmap_legend_param = list(
                       at = c(0, 1),
                       labels = c("n.s.", UP_DOWN),
                       title = "HOG deregulation",
                       legend_height = unit(4, "cm")
                     ))
      
      ht1 <- draw(ht1)
      dev.off()
      
      pdffile = paste0("./results/Heatmap_",Node,"_Node_Mimpud_",cond,"_conditions",UP_DOWN,".pdf")
      #adjust the size and margin of the pdf file
      pdf(pdffile,height=8,width=8)
      # par(mar=c(5,5,5,5))
      draw(ht1)
      dev.off()
      
      
      Node_df <- eval(parse(text = paste0("HOG_all_df %>% dplyr::filter(Node_label_",UP_DOWN,"==Node)")))
      Node_df[Node_df==UP_DOWN] <- 1
      Node_df[Node_df==""] <- 0
      dat <- as.data.frame(dplyr::select(Node_df,all_of(paste0(Conditions)))%>% mutate_if(is.character, as.numeric))
      rownames(dat) <- Node_df$HOG
      
      upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                         order.by = "freq",keep.order = TRUE, 
                         mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))
      
      pdf(file=paste0("./results/upsetR_",Node,"_Node_Mimpud_",trait,"_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
      print(upsetPlot)
      grid.text(paste0("Overlaping ",UP_DOWN," HOG among ",trait," conditions in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
      dev.off()
      
      
    }
    
    HOG_all_df <- Mimpud_NODE_genes_Up %>% group_by(Node_label) %>% dplyr::filter(Ctai==1) %>% dplyr::summarise(Overlap=length(Gene_id))
    HOG_all_df$Trait <- "Ctai"
    for(cond in Conditions){
      cond_df <- eval(parse(text = paste0("Mimpud_trait %>% group_by(Node_label) %>% dplyr::filter(",cond,"==1) %>% dplyr::summarise(Overlap=length(Gene_id))")))
      cond_df$Trait <- cond
      HOG_all_df <- rbind(HOG_all_df,cond_df)
    }
    
    Full_summary <- HOG_all_df
    
    Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
    for(i in unique(Full_summary$Trait)){
      Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
      
    }
    
    Nodes <- c("N1","N2","Mimpud")
    for(Node in Nodes){
      for(i in Conditions){
        CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Node_label==Node]
        if(length(CloneNode)==0){CloneNode<-0}
        CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="Ctai" & Full_summary$Node_label==Node]
        M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                      Clones_sum$Sum[Clones_sum$Trait=="Ctai"]-CtaiNode)))
        test <- chisq.test(M)
        fishertest <- fisher.test(M)
        Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Trait==i] <- test$p.value
        Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$p.value
        Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$estimate
        Full_summary$Total[Full_summary$Node_label=="Mimpud" & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
      }
    }
    
    Full_summary$Total[Full_summary$Trait=="Ctai" & Full_summary$Node_label=="Mimpud"] <- Clones_sum$Sum[Clones_sum$Trait=="Ctai"]
    
    Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
    Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
    Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
    Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
    Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
    Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
    Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than C.tai"
    Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than C.tai"
    Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than C.tai"
    Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than C.tai"
    Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than C.tai"
    Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than C.tai"
    
    Full_summary$Node_label[Full_summary$Node_label=="Mimpud"] <- "M. pudica specific"
    Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","M. pudica specific"))
    ordered_summary <- Full_summary %>% dplyr::select(Node_label,Trait,Overlap,Total) %>% dplyr::filter(Node_label=="M. pudica specific")
    ordered_summary$Prop <- ordered_summary$Overlap/ordered_summary$Total
    ordered_summary <- ordered_summary %>% arrange(-Prop)
    Full_summary$Trait <- factor(Full_summary$Trait,
                                 levels = unique(ordered_summary$Trait))
    # fwrite(Full_summary,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_Traits_",UP_DOWN,".txt"),sep="\t")
    
    Full_summary_Up <- Full_summary
    mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
    p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Node_label)) +
      geom_bar(stat="identity", position="fill",  width = 0.65) +
      scale_fill_manual(values = mycol,labels = c("NFN node", "Fabales node", "M. pudica")) + 
      labs(y= "Propotion of genes shared with C.tai", x = "Traits") +
      geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
      geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    pdf(file = paste0(out_path,"../trait_results/Barplot_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_",trait,"_",UP_DOWN,".pdf"), width=15, height=8)
    print(p)
    dev.off()
  }
}
 



######## Mimpud traits ########
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))

NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")

for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",UP_DOWN,".tsv"))
  
  # Load DEGs Mimpud WT and fonctional traits
  CtaiWT_AllUp <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(CtaiWT_AllUp) <- "Gene_id"
  CtaiWT_AllUp <- left_join(CtaiWT_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_AllUp <- unique(CtaiWT_AllUp$Gene_id[!is.na(CtaiWT_AllUp$HOG)])
  
  CtaiWT_Other <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",setdiff(UP_or_DOWN,UP_DOWN),".txt"),header=FALSE)
  names(CtaiWT_Other) <- "Gene_id"
  CtaiWT_Other <- left_join(CtaiWT_Other,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_Other <- unique(CtaiWT_Other$Gene_id[!is.na(CtaiWT_Other$HOG)])
  
  HOG_id_NODE <- Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN),starts_with("Mimpud_"))
  names(HOG_id_NODE) <- c("HOG","Node_label","Mimpud")
  HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% NODES)
  Mimpud_NODE_genes <- left_join(HOG_id_NODE,Mimpud_cor_Genes_HOG,by="HOG")
  
  # Overlap CtaiWT_UP_DOWN and NODE
  Mimpud_NODE_genes_Up <- Mimpud_NODE_genes %>% dplyr::filter(Gene_id %in% CtaiWT_AllUp)
  Mimpud_NODE_genes_Up$Ctai <- 1
  
  TRAITS = c("NF","Organogenesis","Release","Persist","NFix")
  for(trait in TRAITS){
    trait_df <- fread(input = paste0("./trait_analysis/",trait,".txt"),header=TRUE)
    names(trait_df)[1] <- "Gene_id"
    trait_df <- trait_df %>% dplyr::filter(Up_Down==paste0(UP_DOWN)) %>% dplyr::select(Gene_id)
    trait_df$trait <- "1";names(trait_df) <- c("Gene_id",trait)
    Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,trait_df,by="Gene_id")
    Mimpud_NODE_genes_Up <- Mimpud_NODE_genes_Up %>% dplyr::distinct()
  }
  Mimpud_NODE_genes_Up[is.na(Mimpud_NODE_genes_Up)] <- 0
  
  

  HOG_all_df <- Mimpud_NODE_genes_Up %>% group_by(Node_label) %>% dplyr::filter(Ctai==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "Ctai"
  for(trait in TRAITS){
    trait_df <- eval(parse(text = paste0("Mimpud_NODE_genes_Up %>% group_by(Node_label) %>% dplyr::filter(",trait,"==1) %>% dplyr::summarise(Overlap=length(Gene_id))")))
    trait_df$Trait <- trait
    HOG_all_df <- rbind(HOG_all_df,trait_df)
  }
  
  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  
  Nodes <- c("N1","N2","Mimpud")
  for(Node in Nodes){
    for(i in TRAITS){
      CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Node_label==Node]
      CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="Ctai" & Full_summary$Node_label==Node]
      M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                    Clones_sum$Sum[Clones_sum$Trait=="Ctai"]-CtaiNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Trait==i] <- test$p.value
      Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$estimate
      Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i]
    }
  }
  
  Full_summary$Total[Full_summary$Trait=="Ctai" & Full_summary$Node_label=="N1"] <- Clones_sum$Sum[Clones_sum$Trait=="Ctai"]
  
  Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
  Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
  Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than C.tai"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than C.tai"

  Full_summary$Node_label[Full_summary$Node_label=="Mimpud"] <- "M. pudica specific"
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","M. pudica specific"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Ctai","NF","Organogenesis","Release","Persist","NFix"))
  fwrite(Full_summary,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_Traits_",UP_DOWN,".txt"),sep="\t")
  
  Full_summary_Up <- Full_summary
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  p <- ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("NFN node", "Fabales node", "M. pudica")) + 
    labs(y= "Propotion of genes shared with C.tai", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 18)
  
  pdf(file = paste0(out_path,"../trait_results/Barplot_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_Traits_",UP_DOWN,".pdf"), width=18, height=8)
  print(p)
  dev.off()
  
}



######## Mimpud traits by HOG ########
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))

NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")

for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",UP_DOWN,".tsv"))
  
  # Load DEGs Mimpud WT and fonctional traits
  CtaiWT_AllUp <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(CtaiWT_AllUp) <- "Gene_id"
  CtaiWT_AllUp <- left_join(CtaiWT_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_AllUp <- unique(CtaiWT_AllUp$Gene_id[!is.na(CtaiWT_AllUp$HOG)])
  
  CtaiWT_Other <- fread(input = paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_",setdiff(UP_or_DOWN,UP_DOWN),".txt"),header=FALSE)
  names(CtaiWT_Other) <- "Gene_id"
  CtaiWT_Other <- left_join(CtaiWT_Other,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_Other <- unique(CtaiWT_Other$Gene_id[!is.na(CtaiWT_Other$HOG)])
  
  HOG_id_NODE <- Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN),starts_with("Mimpud_"))
  names(HOG_id_NODE) <- c("HOG","Node_label","Mimpud")
  HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% NODES)
  Mimpud_NODE_genes <- left_join(HOG_id_NODE,Mimpud_cor_Genes_HOG,by="HOG")
  
  # Overlap CtaiWT_UP_DOWN and NODE
  Mimpud_NODE_genes_Up <- Mimpud_NODE_genes %>% dplyr::filter(Gene_id %in% CtaiWT_AllUp)
  Mimpud_NODE_genes_Up$Ctai <- 1
  
  TRAITS = c("NF_inhib_def","NF","Organogenesis","Release","Persist","NFix")
  for(trait in TRAITS){
    trait_df <- fread(input = paste0("./trait_analysis/",trait,".txt"),header=TRUE)
    names(trait_df)[1] <- "Gene_id"
    trait_df <- trait_df %>% dplyr::filter(Up_Down==paste0(UP_DOWN)) %>% dplyr::select(Gene_id)
    trait_df$trait <- "1";names(trait_df) <- c("Gene_id",trait)
    Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,trait_df,by="Gene_id")
    Mimpud_NODE_genes_Up <- Mimpud_NODE_genes_Up %>% dplyr::distinct()
  }
  Mimpud_NODE_genes_Up[is.na(Mimpud_NODE_genes_Up)] <- 0
  
  HOG_all_df <- Mimpud_NODE_genes_Up %>% distinct(HOG, .keep_all = TRUE) %>% dplyr::select(HOG,Ctai)
  for(trait in TRAITS){
    trait_df <- eval(parse(text = paste0("Mimpud_NODE_genes_Up %>% group_by(HOG) %>% dplyr::summarise(Overlap=max(",trait,"))")))
    # trait_df <- eval(parse(text = paste0("Mimpud_NODE_genes_Up %>% group_by(HOG) %>% dplyr::filter(",trait,"==1) %>% dplyr::summarise(Overlap=max(",trait,"))")))
    names(trait_df) <- c("HOG",paste0("Mimpud_",trait))

    HOG_all_df <- dplyr::left_join(HOG_all_df,trait_df,by="HOG")
  }
  HOG_all_df[HOG_all_df==1] <- UP_DOWN
  HOG_all_df[HOG_all_df==0] <- ""
  
  HOG_all_df <- left_join(Recap_HOG,HOG_all_df,by="HOG")
  HOG_all_df[is.na(HOG_all_df)] <- ""
  
  fwrite(HOG_all_df,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits",UP_DOWN,".tsv"), quote=FALSE,sep="\t")
  
  Nodes <- c("N1","N2","Mimpud")
  for(Node in Nodes){
    Node_df <- eval(parse(text = paste0("HOG_all_df %>% dplyr::filter(Node_label_",UP_DOWN,"==Node)")))
    Node_df[Node_df==UP_DOWN] <- 1
    Node_df[Node_df==""] <- 0
    Node_df <- Node_df %>% tidyr::unite("PASTE",all_of(paste0("Mimpud_",TRAITS)), sep="_", remove=FALSE)
    Node_df <- Node_df %>% dplyr::arrange(desc(PASTE))
    
    #####
    dat <- as.matrix(dplyr::select(Node_df,all_of(paste0("Mimpud_",TRAITS)))%>% mutate_if(is.character, as.numeric)) 
    rownames(dat) <- Node_df$HOG
    
    if(UP_DOWN=="Up"){my_palette <- c(colorRampPalette(colors=c("cornsilk2","firebrick"))(n=2))
    }else{my_palette <- c(colorRampPalette(colors=c("cornsilk2","royalblue"))(n=2))}
    
    ht1 <- Heatmap(dat, name = "HOG", cluster_rows = FALSE,
                   row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
                   cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
                   row_names_side = "left",row_order=rownames(dat), row_title = NULL,
                   row_split = factor(Node_df$PASTE, levels = unique(Node_df$PASTE)),
                   show_row_names = FALSE, col = my_palette , width = unit(25, "mm"), height = unit(14, "cm"),
                   # row_km = length(unique(Node_df$PASTE)),
                   row_title_rot = 0 ,heatmap_legend_param = list(
                     at = c(0, 1),
                     labels = c("n.s.", UP_DOWN),
                     title = "HOG deregulated",
                     legend_height = unit(4, "cm")
                   ))
    
    ht1 <- draw(ht1)
    dev.off()
    
    pdffile = paste0("./results/Heatmap_",Node,"_Node_Mimpud_traits_",UP_DOWN,".pdf")
    #adjust the size and margin of the pdf file
    pdf(pdffile,height=8,width=8)
    # par(mar=c(5,5,5,5))
    draw(ht1)
    dev.off()
    
    # rcl_list <- melt(row_order(ht1))
    # names(rcl_list) <- c("ID","Cluster")
    # rcl_list$Cluster <- paste0("Cluster",rcl_list$Cluster)
    # rcl_list$HOG <-  Node_df$HOG[rcl_list$ID]
    # 
    # Node_df_recap <- left_join(Node_df,rcl_list,by="HOG")
    # Node_df_recap <- dplyr::select(Node_df_recap,-ID)
    # 
    # # fwrite(Node_df_recap,paste0("./results/Summary_cluster_heatmap_NFN_Node_mimpud_traits_",UP_DOWN,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE,sep="\t")
    # 
    # for(i in unique(Node_df_recap$Cluster)){
    #   Mimpud_NODE_genes_cluster <- dplyr::filter(Node_df_recap,Cluster==i)
    #   dat <- as.matrix(dplyr::select(Mimpud_NODE_genes_cluster,-HOG,-Node_label,-Mimpud,-Gene_id,-Ctai,-Ranking,-Acronym,-Category,
    #                                  -Description,-Mean_seq_number_HOG,-Mean_seq_number_NodspecieswithDEGs,-Cluster,-Sampling_percentage))
    #   rownames(dat) <- Mimpud_NODE_genes_cluster$Gene_id
    #   
    #   ht_cluster <- Heatmap(dat, name = "HOG", cluster_rows = TRUE, row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 9),
    #                         cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE, column_names_rot = 90, column_names_side = "top",
    #                         row_names_side = "left", show_row_names = FALSE, col = my_palette , width = unit(25, "mm"), height = unit(15, "cm"),
    #                         row_km = 1,row_title_rot = 0,heatmap_legend_param = list(
    #                           at = c(0, 1),
    #                           labels = c("n.s.", UP_DOWN),
    #                           title = "Gene deregulation",
    #                           legend_height = unit(4, "cm")
    #                         )) +
    #     rowAnnotation(Gene_names = anno_mark(at=c(which(Mimpud_NODE_genes_cluster$Acronym != "")),labels = Mimpud_NODE_genes_cluster$Acronym[Mimpud_NODE_genes_cluster$Acronym != ""]))
    #   
    #   pdffile = paste0("./results/Heatmap_NFN_Node_mimpud_traits_",UP_DOWN,"_",i,".pdf")
    #   #adjust the size and margin of the pdf file
    #   pdf(pdffile,height=20,width=15)
    #   # par(mar=c(5,5,5,5))
    #   draw(ht_cluster)
    #   dev.off()
    #   
    # }
    # 
    
    Node_df <- eval(parse(text = paste0("HOG_all_df %>% dplyr::filter(Node_label_",UP_DOWN,"==Node)")))
    Node_df[Node_df==UP_DOWN] <- 1
    Node_df[Node_df==""] <- 0
    dat <- as.data.frame(dplyr::select(Node_df,all_of(paste0("Mimpud_",TRAITS)))%>% mutate_if(is.character, as.numeric))
    rownames(dat) <- Node_df$HOG
    
    upsetPlot <- upset(dat, sets = unique(names(dat)),point.size = 3.2, line.size = 1.1, mb.ratio = c(0.7, 0.3), nintersects = 91,
                       order.by = "freq",keep.order = TRUE, 
                       mainbar.y.label = "Intersections of HOG",matrix.color = "grey20", sets.x.label = "Number of HOG",text.scale = c(1.7, 1.7, 1.2, 1.2, 1.6, 1.2))

    pdf(file=paste0("./results/upsetR_",Node,"_Node_Mimpud_traits_",UP_DOWN,".pdf"),width = 12, height = 8, onefile=FALSE)
    print(upsetPlot)
    grid.text(paste0("Overlaping ",UP_DOWN," HOG among traits in ",Node," node."),x = 0.65, y=0.95, gp=gpar(fontsize=12))
    dev.off()
    
    
    
}
}
