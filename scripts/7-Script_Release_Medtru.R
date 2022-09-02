# Script to count HOGs on species tree
library("data.table");library("ape");library("seqinr");library("adephylo");library("tidyverse");library("tidyr");library("VennDiagram")
library("UpSetR");library("gaston");library("bio3d");library("phytools");library("ggtree");library("phylobase");library("reshape2")
library("ggpubr");library("wesanderson")
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

######## Medicago zones ########
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Medtru_"))

NODES <- c("N1","N2","N7","N12","N21","Medtru")
UP_or_DOWN <- c("Up","Down")

GENOTYPES <- c("FIIp","FIId","FI","IZ","ZIII")

for(UP_DOWN in UP_or_DOWN){
  FIIp <- fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_FIIp_LogFC1.5.txt"),h=FALSE)$V1
  FIId <- fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_FIId_LogFC1.5.txt"),h=FALSE)$V1
  FI <- fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_FI_LogFC1.5.txt"),h=FALSE)$V1
  IZ <- fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_IZ_LogFC1.5.txt"),h=FALSE)$V1
  ZIII <- fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_ZIII_LogFC1.5.txt"),h=FALSE)$V1
  
  FIIp_spe <- setdiff(FIIp,FIId);FIIp_spe <- setdiff(FIIp_spe,FI);FIIp_spe <- setdiff(FIIp_spe,IZ)
  FIIp_spe <- data.frame(Gene_id=setdiff(FIIp_spe,ZIII),FIIp=1)
  FIIp_spe <- data.frame(Gene_id=FIIp,FIIp=1)
  
  FIId_spe <- setdiff(FIId,FIIp);FIId_spe <- setdiff(FIId_spe,FI);FIId_spe <- setdiff(FIId_spe,IZ)
  FIId_spe <- data.frame(Gene_id=setdiff(FIId_spe,ZIII),FIId=1)
  FIId_spe <- data.frame(Gene_id=FIId,FIId=1)
  
  FI_spe <- setdiff(FI,FIIp);FI_spe <- setdiff(FI_spe,FIId);FI_spe <- setdiff(FI_spe,IZ)
  FI_spe <- data.frame(Gene_id=setdiff(FI_spe,ZIII),FI=1)
  FI_spe <- data.frame(Gene_id=FI,FI=1)
  
  IZ_spe <- setdiff(IZ,FIIp);IZ_spe <- setdiff(IZ_spe,FIId);IZ_spe <- setdiff(IZ_spe,FI)
  IZ_spe <- data.frame(Gene_id=setdiff(IZ_spe,ZIII),IZ=1)
  IZ_spe <- data.frame(Gene_id=IZ,IZ=1)  
  ZIII_spe <- setdiff(ZIII,FIIp);ZIII_spe <- setdiff(ZIII_spe,FIId);ZIII_spe <- setdiff(ZIII_spe,FI)
  ZIII_spe <- data.frame(Gene_id=setdiff(ZIII_spe,IZ),ZIII=1)
  ZIII_spe <- data.frame(Gene_id=ZIII,ZIII=1)
  
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF",UP_DOWN,".tsv"))
  
  # Load DEGs Mimpud WT and fonctional traits
  CtaiWT_AllUp <- fread(input = paste0("./rnaseq/Medtru_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(CtaiWT_AllUp) <- "Gene_id"
  CtaiWT_AllUp <- left_join(CtaiWT_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_AllUp <- unique(CtaiWT_AllUp$Gene_id[!is.na(CtaiWT_AllUp$HOG)])
  
  CtaiWT_Other <- fread(input = paste0("./rnaseq/Medtru_Nod_FDR005_logFC1.5_",setdiff(UP_or_DOWN,UP_DOWN),".txt"),header=FALSE)
  names(CtaiWT_Other) <- "Gene_id"
  CtaiWT_Other <- left_join(CtaiWT_Other,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT_Other <- unique(CtaiWT_Other$Gene_id[!is.na(CtaiWT_Other$HOG)])
  
  HOG_id_NODE <- Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN),starts_with("Medtru_Nod_"))
  names(HOG_id_NODE) <- c("HOG","Node_label","Medtru")
  HOG_id_NODE <- HOG_id_NODE %>% dplyr::filter(Node_label %in% NODES)
  Mimpud_NODE_genes <- left_join(HOG_id_NODE,Mimpud_cor_Genes_HOG,by="HOG")
  
  # Overlap CtaiWT_UP_DOWN and NODE
  Mimpud_NODE_genes_Up <- Mimpud_NODE_genes %>% dplyr::filter(Gene_id %in% CtaiWT_AllUp)
  Mimpud_NODE_genes_Up$Medtru <- 1
  
  Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,FIId_spe,by="Gene_id")
  Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,FIIp_spe,by="Gene_id")
  Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,FI_spe,by="Gene_id")
  Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,ZIII_spe,by="Gene_id")
  Mimpud_NODE_genes_Up <- dplyr::left_join(Mimpud_NODE_genes_Up,IZ_spe,by="Gene_id")
  
  Mimpud_NODE_genes_Up[is.na(Mimpud_NODE_genes_Up)] <- 0
  
  HOGinfo <- dplyr::select(Mimpud_NODE_genes_Up,HOG,FI,FIIp,FIId,IZ,ZIII)
  HOGinfo <- HOGinfo %>% group_by(HOG) %>% dplyr::summarise(FI=sum(FI),FIIp=sum(FIIp),FIId=sum(FIId),IZ=sum(IZ),ZIII=sum(ZIII))
  HOGinfo <- column_to_rownames(HOGinfo, var = "HOG")
  HOGinfo[HOGinfo!="0"] <- UP_DOWN;HOGinfo[HOGinfo=="0"] <- ""
  HOGinfo <- rownames_to_column(HOGinfo, var = "HOG")
  
  HOG_all_df <- left_join(Recap_HOG,HOGinfo,by="HOG")
  HOG_all_df[is.na(HOG_all_df)] <- ""
  HOG_all_df <- as_tibble(HOG_all_df)
  HOG_all_df <- HOG_all_df %>% distinct()
  fwrite(HOG_all_df,paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_Mimpud_traits_MtruLjapMpudNF_and_MedtruZones",UP_DOWN,".tsv"), quote=FALSE,sep="\t")
  
  HOG_all_df <- Mimpud_NODE_genes_Up %>% group_by(Node_label) %>% dplyr::filter(Medtru==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  HOG_all_df$Trait <- "Medtru"
  for(trait in GENOTYPES){
    trait_df <- eval(parse(text = paste0("Mimpud_NODE_genes_Up %>% group_by(Node_label) %>% dplyr::filter(Medtru==1 & ",trait,"==1) %>% dplyr::summarise(Overlap=length(Gene_id))")))
    trait_df$Trait <- trait
    HOG_all_df <- rbind(HOG_all_df,trait_df)
  }

  Full_summary <- HOG_all_df
  
  Clones_sum <- data.frame(Trait=unique(Full_summary$Trait),Sum=0)
  for(i in unique(Full_summary$Trait)){
    Clones_sum$Sum[Clones_sum$Trait==i] <- sum(Full_summary$Overlap[Full_summary$Trait==i])
    
  }
  
  for(Node in NODES){
    for(i in GENOTYPES){
      CloneNode <- Full_summary$Overlap[Full_summary$Trait==i & Full_summary$Node_label==Node]
      CtaiNode <- Full_summary$Overlap[Full_summary$Trait=="Medtru" & Full_summary$Node_label==Node]
      M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Trait==i]-CloneNode,
                                                    Clones_sum$Sum[Clones_sum$Trait=="Medtru"]-CtaiNode)))
      test <- chisq.test(M)
      fishertest <- fisher.test(M)
      Full_summary$Chi2[Full_summary$Node_label==Node & Full_summary$Trait==i] <- test$p.value
      Full_summary$FishTest[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label==Node & Full_summary$Trait==i] <- fishertest$estimate
      Full_summary$Total[Full_summary$Node_label=="N1" & Full_summary$Trait==i] <- Clones_sum$Sum[Clones_sum$Trait==i] 
    }
  }
  
  Full_summary$Total[Full_summary$Trait=="Medtru" & Full_summary$Node_label=="N1"] <- Clones_sum$Sum[Clones_sum$Trait=="Medtru"]
  
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

  Full_summary$Node_label[Full_summary$Node_label=="Medtru"] <- "M. truncatula"
  Full_summary$Node_label <- factor(Full_summary$Node_label,levels = c("N1","N2","N7","N12","N21","M. truncatula"))
  Full_summary$Trait <- factor(Full_summary$Trait,
                               levels = c("Medtru","FI","FIId","FIIp","IZ","ZIII"))
  Full_summary_Up <- Full_summary
  
  fwrite(Full_summary_Up,paste0(out_path,"../trait_results/Summary_stats_Proportion_shared_with_Medtru_release_multi_Nodes_Medtru_Zones_Traits_",UP_DOWN,".txt"),sep="\t")
  
  
  mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]
  mycol <- c("N1" = "#3B9AB2", "N2" = "#78B7C5", "N7"="darkseagreen2",
             "N12" = "darkseagreen3","N21" = "darkseagreen4", "M. truncatula" = "#E1AF00")
  
  pdf(file = paste0(out_path,"../trait_results/Barplot_Proportion_shared_with_Medtru_release_multi_Nodes_Medtru_Zones_Traits_",UP_DOWN,".pdf"), width=15, height=10)
  print(ggplot(data=Full_summary, aes(x=Trait, y=Overlap, fill=Node_label)) +
    geom_bar(stat="identity", position="fill",  width = 0.65) +
    scale_fill_manual(values = mycol,labels = c("NFN clade", "Fabales clade", "Papilionoideae clade",
                                                "Dalbergioid+Hologalegina clade","Hologalegina clade", "M. truncatula")) + 
    #scale_fill_manual(values = c("#3B9AB2","#78B7C5","darkseagreen2","darkseagreen3", "darkseagreen4","#E1AF00")) +
    labs(y= "Number of genes shared with M.truncatula", x = "Traits") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) + theme_classic(base_size = 16))
    dev.off()
  
}

