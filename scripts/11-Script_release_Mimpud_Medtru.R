# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","reshape2","ComplexHeatmap","ggVennDiagram","mixOmics",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ggpubr","wesanderson","DECIPHER","Biostrings","gridExtra","ggvenn","ggdensity"), require, character.only = TRUE)

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

######## Comparaison species specific Release vs RNS response ########
UP_or_DOWN <- c("Up")
UP_DOWN="Up"

  Recap_HOG <- fread(paste0("./results/Nod_final_table_Up.txt"))
  Recap_Medtru <- fread(paste0("./results/Medtru_Nod_final_table.txt"))
  Recap_Medtru <- Recap_Medtru %>% dplyr::filter(Node_label_Up=="Medtru")
  Recap_Mimpud <- fread(paste0("./results/Mimpud_Nod_final_table.txt"))
  Recap_Mimpud <- Recap_Mimpud %>% dplyr::filter(Node_label_Up=="Mimpud")
  
  Medtru_Nod <-   fread(paste0("./rnaseq/Medtru_Nod_FDR005_logFC1.5_Up.txt"),h=FALSE)
  names(Medtru_Nod) <- "Gene_id";Medtru_Nod$Medtru <- "Up"
  Medtru_Nod <- Medtru_Nod %>% dplyr::filter(Gene_id %in% Recap_Medtru$Gene_id)
  Mimpud_Nod <- fread(paste0("./rnaseq/Mimpud_Nod_FDR005_logFC1.5_Up.txt"),h=FALSE)
  names(Mimpud_Nod) <- "Gene_id";Mimpud_Nod$Mimpud <- "Up"
  Mimpud_Nod <- Mimpud_Nod %>% dplyr::filter(Gene_id %in% Recap_Mimpud$Gene_id)
  
  #### 1) differences in signalP #### 
  # Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_infos.txt"), sep="\t", h=TRUE, row.names=1)
  Medtru_prot <- read.table(paste0("./results_NCR/Medtru_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
  Medtru_prot <- dplyr::select(Medtru_prot,-Acronym)
  Medtru_prot$Gene_id <- row.names(Medtru_prot)
  Medtru_prot <- Medtru_prot %>% dplyr::filter(Gene_id %in% Recap_Medtru$Gene_id)
  Medtru_prot$SignalP[Medtru_prot$SignalP=="signal_peptide"] <- 1
  Medtru_prot$SignalP[Medtru_prot$SignalP==""] <- 0; Medtru_prot$SignalP <- as.numeric(Medtru_prot$SignalP)
  
  # FIId <- data.frame(Gene_id=fread(paste0("./nod_laser_medtru/All_",UP_DOWN,"_Medtru_FIId_LogFC1.5.txt"),h=FALSE)$V1,FIId="Up")
  FIId <- Recap_Medtru %>% dplyr::select(Gene_id,Medtru_FIId_Up) %>% dplyr::filter(Medtru_FIId_Up=="Up")
  names(FIId) <- c("Gene_id","FIId")
  FIId <- FIId %>% distinct()
  
  Medtru_df <- dplyr::left_join(Medtru_Nod,FIId,by="Gene_id")
  Medtru_df <- dplyr::left_join(Medtru_df,Medtru_prot,by="Gene_id")
  Medtru_df[is.na(Medtru_df)] <- ""
  
  Full_summary <- data.frame(Class=c("Medtru","FIId","Medtru","FIId"),SignalP=c("noSignalP","noSignalP","SignalP","SignalP"),
                             Overlap=c(length(Medtru_df$Gene_id[Medtru_df$Medtru=="Up"&Medtru_df$FIId==""&Medtru_df$SignalP==0]),
                                       length(Medtru_df$Gene_id[Medtru_df$Medtru=="Up"&Medtru_df$FIId=="Up"&Medtru_df$SignalP==0]),
                                       length(Medtru_df$Gene_id[Medtru_df$Medtru=="Up"&Medtru_df$FIId==""&Medtru_df$SignalP==1]),
                                       length(Medtru_df$Gene_id[Medtru_df$Medtru=="Up"&Medtru_df$FIId=="Up"&Medtru_df$SignalP==1])))
  M <- data.frame(noSignalP=c(Full_summary$Overlap[Full_summary$Class=="Medtru"&Full_summary$SignalP=="noSignalP"],
                              Full_summary$Overlap[Full_summary$Class=="FIId"&Full_summary$SignalP=="noSignalP"]),
                  SignalP=c(Full_summary$Overlap[Full_summary$Class=="Medtru"&Full_summary$SignalP=="SignalP"],
                            Full_summary$Overlap[Full_summary$Class=="FIId"&Full_summary$SignalP=="SignalP"]))
  row.names(M) <-c("Medtru","FIId")
  Signal_P_fishertest <- fisher.test(M)
  
  Full_summary$FishTest[Full_summary$Class=="FIId" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$p.value
  Full_summary$oddsratio[Full_summary$Class=="FIId" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$estimate
  
  Full_summary$Total[Full_summary$Class=="FIId" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="FIId"])
  Full_summary$Total[Full_summary$Class=="Medtru" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="Medtru"])
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*";Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**";Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=SignalP)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("cornflowerblue","coral1")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  p
  dev.off()
  
  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Medru_FIId_signalP_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Medru_FIId_signalP_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  
  # Mimpud Release
  Mimpud_prot <- read.table(paste0("./results_NCR/Mimpud_proteins_amino_acid_info.txt"), sep="\t", h=TRUE, row.names=1)
  Mimpud_prot <- dplyr::select(Mimpud_prot,-Acronym)
  Mimpud_prot$Gene_id <- row.names(Mimpud_prot)
  Mimpud_prot <- Mimpud_prot %>% dplyr::filter(Gene_id %in% Recap_Mimpud$Gene_id)
  Mimpud_prot$SignalP[Mimpud_prot$SignalP=="signal_peptide"] <- 1
  Mimpud_prot$SignalP[Mimpud_prot$SignalP==""] <- 0; Mimpud_prot$SignalP <- as.numeric(Mimpud_prot$SignalP)
  
  Release <- Recap_Mimpud %>% dplyr::select(Gene_id,Mimpud_Release_Up) %>% dplyr::filter(Mimpud_Release_Up=="Up")
  names(Release) <- c("Gene_id","Release")
  Release <- Release %>% distinct()
  Persist <- Recap_Mimpud %>% dplyr::select(Gene_id,Mimpud_Persist_Up) %>% dplyr::filter(Mimpud_Persist_Up=="Up")
  names(Persist) <- c("Gene_id","Persist")
  Persist <- Persist %>% distinct()
  
  Mimpud_df <- dplyr::left_join(Mimpud_Nod,Release,by="Gene_id")
  Mimpud_df <- dplyr::left_join(Mimpud_df,Persist,by="Gene_id")
  Mimpud_df <- dplyr::left_join(Mimpud_df,Mimpud_prot,by="Gene_id")
  Mimpud_df[is.na(Mimpud_df)] <- ""
  
  Full_summary <- data.frame(Class=c("Mimpud","Release","Mimpud","Release"),SignalP=c("noSignalP","noSignalP","SignalP","SignalP"),
                             Overlap=c(length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Release==""&Mimpud_df$SignalP==0]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Release=="Up"&Mimpud_df$SignalP==0]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Release==""&Mimpud_df$SignalP==1]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Release=="Up"&Mimpud_df$SignalP==1])))
  M <- data.frame(noSignalP=c(Full_summary$Overlap[Full_summary$Class=="Mimpud"&Full_summary$SignalP=="noSignalP"],
                              Full_summary$Overlap[Full_summary$Class=="Release"&Full_summary$SignalP=="noSignalP"]),
                  SignalP=c(Full_summary$Overlap[Full_summary$Class=="Mimpud"&Full_summary$SignalP=="SignalP"],
                            Full_summary$Overlap[Full_summary$Class=="Release"&Full_summary$SignalP=="SignalP"]))
  row.names(M) <-c("Mimpud","Release")
  Signal_P_fishertest <- fisher.test(M)
  
  Full_summary$FishTest[Full_summary$Class=="Release" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$p.value
  Full_summary$oddsratio[Full_summary$Class=="Release" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$estimate
  
  Full_summary$Total[Full_summary$Class=="Release" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="Release"])
  Full_summary$Total[Full_summary$Class=="Mimpud" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud"])
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=SignalP)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("cornflowerblue","coral1")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  p
  dev.off()
  
  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  # Mimpud Persist
  Full_summary <- data.frame(Class=c("Mimpud","Persist","Mimpud","Persist"),SignalP=c("noSignalP","noSignalP","SignalP","SignalP"),
                             Overlap=c(length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Persist==""&Mimpud_df$SignalP==0]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Persist=="Up"&Mimpud_df$SignalP==0]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Persist==""&Mimpud_df$SignalP==1]),
                                       length(Mimpud_df$Gene_id[Mimpud_df$Mimpud=="Up"&Mimpud_df$Persist=="Up"&Mimpud_df$SignalP==1])))
  M <- data.frame(noSignalP=c(Full_summary$Overlap[Full_summary$Class=="Mimpud"&Full_summary$SignalP=="noSignalP"],
                              Full_summary$Overlap[Full_summary$Class=="Persist"&Full_summary$SignalP=="noSignalP"]),
                  SignalP=c(Full_summary$Overlap[Full_summary$Class=="Mimpud"&Full_summary$SignalP=="SignalP"],
                            Full_summary$Overlap[Full_summary$Class=="Persist"&Full_summary$SignalP=="SignalP"]))
  row.names(M) <-c("Mimpud","Persist")
  Signal_P_fishertest <- fisher.test(M)
  
  Full_summary$FishTest[Full_summary$Class=="Persist" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$p.value
  Full_summary$oddsratio[Full_summary$Class=="Persist" & Full_summary$SignalP=="SignalP"] <- Signal_P_fishertest$estimate
  
  Full_summary$Total[Full_summary$Class=="Persist" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="Persist"])
  Full_summary$Total[Full_summary$Class=="Mimpud" & Full_summary$SignalP=="SignalP"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud"])
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=SignalP)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("cornflowerblue","coral1")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)
  p
  dev.off()
  
  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  #### 2) compare signalP sequences #### 
  # Medtru
  Medtru_df_signalP <- Medtru_df %>% dplyr::filter(SignalP==1)
  categories <- names(Medtru_df_signalP)[c(-1,-2,-3,-25)]
  Stats_summary=data.frame(Category=categories)
  for(i in categories){
    sub_medtru <- as.vector(as.data.frame(Medtru_df_signalP %>% dplyr::filter(Medtru=="Up"&FIId=="") %>% dplyr::select(paste0(i))))
    sub_trait <- as.vector(as.data.frame(Medtru_df_signalP %>% dplyr::filter(Medtru=="Up"&FIId=="Up") %>% dplyr::select(paste0(i))))
    t_test <- t.test(as.numeric(sub_trait[,1]),as.numeric(sub_medtru[,1]))
    Stats_summary$ttest_p.val[Stats_summary$Category==i] <- t_test$p.value
    Stats_summary$Mean_Release[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))
    Stats_summary$Mean_Mimpud[Stats_summary$Category==i] <- mean(as.numeric(sub_medtru[,1]))
    Stats_summary$Ratio_means[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))/mean(as.numeric(sub_medtru[,1]))
  }
  
  Stats_summary <- arrange(Stats_summary, ttest_p.val)
  Stats_summary$ttest_p.val <- p.adjust(Stats_summary$ttest_p.val, method="fdr")
  fwrite(Stats_summary,paste0("./results/IntracellAccomodation/Barplot_Medtru_FIId_Prot_specificity_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  scatterPlot <- ggplot(Medtru_df_signalP,aes(as.numeric(C_prop), as.numeric(Pep_length), color=FIId)) + 
    geom_point() + expand_limits(x=0)+xlab("Cysteine proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) +
    geom_hline(yintercept=100, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(C_prop), as.numeric(Pep_length), fill=FIId),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 25))+
    scale_y_continuous(limits = c(0, 1000))
  

  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Medtru_df_signalP, aes(as.numeric(C_prop), fill=FIId)) + 
    geom_density(alpha=.5)+xlab("") + 
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 25))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Medtru_df_signalP, aes(y=Pep_length, fill=FIId)) + 
    geom_density(alpha=.5)+ylab("") + 
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=100, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Medtru_signalP_Cysteine_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()

  scatterPlot <- ggplot(Medtru_df_signalP,aes(as.numeric(P_prop), as.numeric(Pep_length), color=FIId)) + 
    geom_point() + expand_limits(x=0)+xlab("Proline proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) +
    geom_hline(yintercept=100, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(P_prop), as.numeric(Pep_length), fill=FIId),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 55))+
    scale_y_continuous(limits = c(0, 1000))
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Medtru_df_signalP, aes(as.numeric(P_prop), fill=FIId)) + 
    geom_density(alpha=.5)+xlab("") + 
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 55))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Medtru_df_signalP, aes(y=Pep_length, fill=FIId)) + 
    geom_density(alpha=.5)+ylab("") + 
    geom_hline(yintercept=100, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Medtru_signalP_Proline_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  # Mimpud Release
  Mimpud_df$Pep_length <- as.numeric(Mimpud_df$Pep_length)
  Mimpud_df_signalP <- Mimpud_df %>% dplyr::filter(SignalP==1)
  # Mimpud_df_signalP$PH_prop <- as.numeric(Mimpud_df_signalP$P_prop)+as.numeric(Mimpud_df_signalP$H_prop)
  # Mimpud_df_signalP$PHD_prop <- (as.numeric(Mimpud_df_signalP$P_prop)+as.numeric(Mimpud_df_signalP$H_prop)+as.numeric(Mimpud_df_signalP$D_prop))
  categories <- names(Mimpud_df_signalP)[c(-1,-2,-3,-4,-26)]
  Stats_summary=data.frame(Category=categories)
  for(i in categories){
    sub_Mimpud <- as.vector(as.data.frame(Mimpud_df_signalP %>% dplyr::filter(Mimpud=="Up"&Release=="") %>% dplyr::select(paste0(i))))
    sub_trait <- as.vector(as.data.frame(Mimpud_df_signalP %>% dplyr::filter(Mimpud=="Up"&Release=="Up") %>% dplyr::select(paste0(i))))
    t_test <- t.test(as.numeric(sub_trait[,1]),as.numeric(sub_Mimpud[,1]))
    Stats_summary$ttest_p.val[Stats_summary$Category==i] <- t_test$p.value
    Stats_summary$Mean_Release[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))
    Stats_summary$Mean_Mimpud[Stats_summary$Category==i] <- mean(as.numeric(sub_Mimpud[,1]))
    Stats_summary$Ratio_means[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))/mean(as.numeric(sub_Mimpud[,1]))
  }
  
  Stats_summary <- arrange(Stats_summary, ttest_p.val)
  Stats_summary$ttest_p.val <- p.adjust(Stats_summary$ttest_p.val, method="fdr")
  fwrite(Stats_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Realease_Prot_specificity_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  Mimpud_df_signalP$P_prop <- as.numeric(Mimpud_df_signalP$P_prop)
  Mimpud_df_signalP$H_prop <- as.numeric(Mimpud_df_signalP$H_prop)
  Mimpud_df_signalP$H_prop <- as.numeric(Mimpud_df_signalP$D_prop)
  
  scatterPlot <- ggplot(Mimpud_df_signalP,aes(as.numeric(P_prop), as.numeric(Pep_length), color=Release)) + 
    geom_point() + expand_limits(x=0)+xlab("Proline proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(P_prop), as.numeric(Pep_length), fill=Release),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 55))+
    scale_y_continuous(limits = c(0, 1000))
  
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Mimpud_df_signalP, aes(as.numeric(P_prop), fill=Release)) + 
    geom_density(alpha=.5)+xlab("") +
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 55))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Mimpud_df_signalP, aes(y=Pep_length, fill=Release)) + 
    geom_density(alpha=.5)+ylab("") + 
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Mimpud_Release_Proline_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  scatterPlot <- ggplot(Mimpud_df_signalP,aes(as.numeric(C_prop), as.numeric(Pep_length), color=Release)) + 
    geom_point() + expand_limits(x=0)+xlab("Cysteine proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(C_prop), as.numeric(Pep_length), fill=Release),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 25))+
    scale_y_continuous(limits = c(0, 1000))
  
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Mimpud_df_signalP, aes(as.numeric(C_prop), fill=Release)) + 
    geom_density(alpha=.5)+xlab("") + 
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 25))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Mimpud_df_signalP, aes(y=Pep_length, fill=Release)) + 
    geom_density(alpha=.5)+ylab("") + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Mimpud_Release_Cysteine_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  # Mimpud Persist 
  Mimpud_df$Pep_length <- as.numeric(Mimpud_df$Pep_length)
  Mimpud_df_signalP <- Mimpud_df %>% dplyr::filter(SignalP==1)
  # Mimpud_df_signalP$PH_prop <- as.numeric(Mimpud_df_signalP$P_prop)+as.numeric(Mimpud_df_signalP$H_prop)
  # Mimpud_df_signalP$PHD_prop <- (as.numeric(Mimpud_df_signalP$P_prop)+as.numeric(Mimpud_df_signalP$H_prop)+as.numeric(Mimpud_df_signalP$D_prop))
  categories <- names(Mimpud_df_signalP)[c(-1,-2,-3,-4,-26)]
  Stats_summary=data.frame(Category=categories)
  for(i in categories){
    sub_Mimpud <- as.vector(as.data.frame(Mimpud_df_signalP %>% dplyr::filter(Mimpud=="Up"&Persist=="") %>% dplyr::select(paste0(i))))
    sub_trait <- as.vector(as.data.frame(Mimpud_df_signalP %>% dplyr::filter(Mimpud=="Up"&Persist=="Up") %>% dplyr::select(paste0(i))))
    t_test <- t.test(as.numeric(sub_trait[,1]),as.numeric(sub_Mimpud[,1]))
    Stats_summary$ttest_p.val[Stats_summary$Category==i] <- t_test$p.value
    Stats_summary$Mean_Release[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))
    Stats_summary$Mean_Mimpud[Stats_summary$Category==i] <- mean(as.numeric(sub_Mimpud[,1]))
    Stats_summary$Ratio_means[Stats_summary$Category==i] <- mean(as.numeric(sub_trait[,1]))/mean(as.numeric(sub_Mimpud[,1]))
  }
  
  Stats_summary <- arrange(Stats_summary, ttest_p.val)
  Stats_summary
  fwrite(Stats_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_Prot_specificity_Up_speciesspecific.txt"),quote=FALSE, sep="\t")
  
  Mimpud_df_signalP$P_prop <- as.numeric(Mimpud_df_signalP$P_prop)
  Mimpud_df_signalP$H_prop <- as.numeric(Mimpud_df_signalP$H_prop)
  
  scatterPlot <- ggplot(Mimpud_df_signalP,aes(as.numeric(P_prop), as.numeric(Pep_length), color=Persist)) + 
    geom_point() + expand_limits(x=0)+xlab("Proline proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(P_prop), as.numeric(Pep_length), fill=Persist),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 55))+
    scale_y_continuous(limits = c(0, 1000))
  
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Mimpud_df_signalP, aes(as.numeric(P_prop), fill=Persist)) + 
    geom_density(alpha=.5)+xlab("") + 
    geom_vline(xintercept=10, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 55))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Mimpud_df_signalP, aes(y=Pep_length, fill=Persist)) + 
    geom_density(alpha=.5)+ylab("") + 
    geom_hline(yintercept=100, linetype="dashed", color = "black", size=0.5) +
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    theme(legend.position = "none")+
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Mimpud_Persist_Proline_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  scatterPlot <- ggplot(Mimpud_df_signalP,aes(as.numeric(C_prop), as.numeric(Pep_length), color=Persist)) + 
    geom_point() + expand_limits(x=0)+xlab("Cysteine proportion")+ylab("Protein length") + 
    scale_color_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position=c(0.9,1), legend.justification=c(1,1)) +
    geom_hdr(method="kde",aes(as.numeric(C_prop), as.numeric(Pep_length), fill=Persist),probs=c(0.6), alpha = .4)+
    scale_x_continuous(limits = c(0, 25))+
    scale_y_continuous(limits = c(0, 1000))
  
  
  # scatterPlot
  # Courbe de densité marginale de x (panel du haut)
  xdensity <- ggplot(Mimpud_df_signalP, aes(as.numeric(C_prop), fill=Persist)) + 
    geom_density(alpha=.5)+xlab("") + 
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    geom_vline(xintercept=6, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position = "none")+
    scale_x_continuous(limits = c(0, 25))
  
  # Courbe de densité marginale de y (panel de droite)
  ydensity <- ggplot(Mimpud_df_signalP, aes(y=Pep_length, fill=Persist)) + 
    geom_density(alpha=.5)+ylab("") + 
    scale_fill_manual(values = c("cornflowerblue","coral1")) + 
    geom_hline(yintercept=150, linetype="dashed", color = "black", size=0.5) +
    theme(legend.position = "none")
    scale_y_continuous(limits = c(0, 1000))
  
  blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(plot.background = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.text.x = element_blank(), 
          axis.text.y = element_blank(), axis.ticks = element_blank())
  
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  pdf(file=paste0("./results/IntracellAccomodation/Biplot_Mimpud_Persist_Cysteine_Protlength_Up_speciesspecific.pdf"),width = 12, height = 8, onefile=FALSE)
  grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
               ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
  dev.off()
  
  #### 3) species specific DEOG vs HOG ####
  # Medtru FIId
  # Specific HOG
  HOG_with_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_")) %>% dplyr::select(HOG)
  HOG_with_Medtru <- unique(HOG_with_Medtru$HOG)
  
  HOG_without_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Medtru <- unique(HOG_without_Medtru$HOG)
  
  Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, HOG_without_Medtru))
  Medtru_full_df <- dplyr::left_join(Medtru_df,Recap_Medtru %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Medtru_spe_NODE_genes <- Medtru_full_df %>% dplyr::filter(HOG %in% Medtru_spe$HOG)
  Medtru_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, Medtru_spe_NODE_genes$HOG))
  Medtru_spe <- Medtru_full_df %>% dplyr::filter(HOG %in% Medtru_spe$HOG)
  
  Medtru_full_df <- rbind(Medtru_spe_NODE_genes, Medtru_spe)
  Medtru_full_df$Node_label_Up[Medtru_full_df$Node_label_Up=="Medtru"] <- "Medtru_DEOG"
  Medtru_full_df$Node_label_Up[Medtru_full_df$Node_label_Up=="Species_specific"] <- "Medtru_specific_OG"
  
  Medtru_full_df$CysRich[Medtru_full_df$FIId=="Up" & Medtru_full_df$SignalP==1 & Medtru_full_df$Pep_length<100 & Medtru_full_df$C_prop>6] <- "CysRich"
  Medtru_full_df[is.na(Medtru_full_df)] <- ""
  
  fwrite(Medtru_full_df %>% dplyr::filter(CysRich=="CysRich"),paste0("./results/IntracellAccomodation/Medtru_FIId_signalP_Up_speciesspecific_C6prot100.txt"),sep="\t")
  
  
  Medtru_global <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global$Class <- "Medtru_WSR"
  Medtru_global_signalP <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global_signalP$Class <- "Medtru_WSR_signalP"
  Medtru_FIId_signalP_NCR <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="CysRich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_FIId_signalP_NCR$Class <- "FIId_CysRich"
  Medtru_global_signalP_woFIId <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="",SignalP==1) %>%
  dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global_signalP_woFIId$Class <- "Medtru_WSR_signalP_woFIIdCysRich"
  
  Full_summary <- rbind(Medtru_global,Medtru_global_signalP)
  Full_summary <- rbind(Full_summary,Medtru_FIId_signalP_NCR)
  Full_summary <- rbind(Full_summary,Medtru_global_signalP_woFIId)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Medtru_WSR"])
    Comparaisons <- "Medtru_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Medtru_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
      
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR_signalP_woFIIdCysRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Medtru_WSR_signalP_woFIIdCysRich"])
    Comparaisons <- "FIId_CysRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Medtru_WSR_signalP_woFIIdCysRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR_signalP_woFIIdCysRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"

  Full_summary$Class <- factor(Full_summary$Class,levels = c("Medtru_WSR","Medtru_WSR_signalP","Medtru_WSR_signalP_woFIIdCysRich","FIId_CysRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Medtru_FIId_signalP_Up_speciesspecific_C6prot100.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Medtru_FIId_signalP_Up_speciesspecific_C6prot100.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  
  # Proline rich Medtru
  HOG_with_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_")) %>% dplyr::select(HOG)
  HOG_with_Medtru <- unique(HOG_with_Medtru$HOG)
  
  HOG_without_Medtru <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Medtru_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Medtru <- unique(HOG_without_Medtru$HOG)
  
  Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, HOG_without_Medtru))
  Medtru_full_df <- dplyr::left_join(Medtru_df,Recap_Medtru %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Medtru_spe_NODE_genes <- Medtru_full_df %>% dplyr::filter(HOG %in% Medtru_spe$HOG)
  Medtru_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Medtru_spe <- data.frame(HOG=setdiff(HOG_with_Medtru, Medtru_spe_NODE_genes$HOG))
  Medtru_spe <- Medtru_full_df %>% dplyr::filter(HOG %in% Medtru_spe$HOG)
  
  Medtru_full_df <- rbind(Medtru_spe_NODE_genes, Medtru_spe)
  Medtru_full_df$Node_label_Up[Medtru_full_df$Node_label_Up=="Medtru"] <- "Medtru_DEOG"
  Medtru_full_df$Node_label_Up[Medtru_full_df$Node_label_Up=="Species_specific"] <- "Medtru_specific_OG"
  
  Medtru_full_df$ProRich[Medtru_full_df$FIId=="Up" & Medtru_full_df$SignalP==1 & Medtru_full_df$Pep_length<100 & Medtru_full_df$P_prop>10] <- "ProRich"
  Medtru_full_df[is.na(Medtru_full_df)] <- ""
  
  Medtru_global <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global$Class <- "Medtru_WSR"
  Medtru_global_signalP <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global_signalP$Class <- "Medtru_WSR_signalP"
  Medtru_FIId_signalP_NCR <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(ProRich=="ProRich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Medtru_FIId_signalP_NCR$Class <- "FIId_ProRich"
  Medtru_global_signalP_woFIId <- Medtru_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(ProRich=="",SignalP==1) %>%
    dplyr::summarise(Overlap=length(Gene_id))
  Medtru_global_signalP_woFIId$Class <- "Medtru_WSR_signalP_woFIIdProRich"
  
  Full_summary <- rbind(Medtru_global,Medtru_global_signalP)
  Full_summary <- rbind(Full_summary,Medtru_FIId_signalP_NCR)
  Full_summary <- rbind(Full_summary,Medtru_global_signalP_woFIId)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Medtru_WSR"])
    Comparaisons <- "Medtru_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Medtru_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR_signalP_woFIIdProRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Medtru_WSR_signalP_woFIIdProRich"])
    Comparaisons <- "FIId_ProRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Medtru_WSR_signalP_woFIIdProRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Medtru_WSR_signalP_woFIIdProRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Medtru_WSR","Medtru_WSR_signalP","Medtru_WSR_signalP_woFIIdProRich","FIId_ProRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Medtru_FIId_signalP_Up_speciesspecific_P10prot100.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Medtru_FIId_signalP_Up_speciesspecific_P10prot100.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  # Mimpud Release
  # Specific HOG
  HOG_with_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_")) %>% dplyr::select(HOG)
  HOG_with_Mimpud <- unique(HOG_with_Mimpud$HOG)
  
  HOG_without_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Mimpud <- unique(HOG_without_Mimpud$HOG)
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, HOG_without_Mimpud))
  Mimpud_full_df <- dplyr::left_join(Mimpud_df,Recap_Mimpud %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Mimpud_spe_NODE_genes <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  Mimpud_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, Mimpud_spe_NODE_genes$HOG))
  Mimpud_spe <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  
  Mimpud_full_df <- rbind(Mimpud_spe_NODE_genes, Mimpud_spe)
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Mimpud"] <- "Mimpud_DEOG"
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Species_specific"] <- "Mimpud_specific_OG"
  
  Mimpud_full_df$Prich[Mimpud_full_df$Release=="Up" & Mimpud_full_df$SignalP==1 & Mimpud_full_df$Pep_length<150 & Mimpud_full_df$P_prop>10] <- "Prorich"
  Mimpud_full_df[is.na(Mimpud_full_df)] <- ""
  
  fwrite(Mimpud_full_df %>% dplyr::filter(Prich=="Prorich"),paste0("./results/IntracellAccomodation/Mimpud_Release_signalP_Up_speciesspecific_P10prot150.txt"),sep="\t")
  
  
  Mimpud_global <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global$Class <- "Mimpud_WSR"
  Mimpud_global_signalP <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP$Class <- "Mimpud_WSR_signalP"
  Mimpud_Release_signalP_Prich <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(Prich=="Prorich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_Release_signalP_Prich$Class <- "Release_ProRich"
  Mimpud_global_signalP_woRelease <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(Prich=="",SignalP==1) %>%
    dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP_woRelease$Class <- "Mimpud_WSR_signalP_woReleaseProRich"
  
  Full_summary <- rbind(Mimpud_global,Mimpud_global_signalP)
  Full_summary <- rbind(Full_summary,Mimpud_Release_signalP_Prich)
  Full_summary <- rbind(Full_summary,Mimpud_global_signalP_woRelease)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR"])
    Comparaisons <- "Mimpud_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woReleaseProRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woReleaseProRich"])
    Comparaisons <- "Release_ProRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woReleaseProRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woReleaseProRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Mimpud_WSR","Mimpud_WSR_signalP","Mimpud_WSR_signalP_woReleaseProRich","Release_ProRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific_P10prot150.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific_P10prot150.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  # Mimpud Persist
  # Specific HOG
  HOG_with_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_")) %>% dplyr::select(HOG)
  HOG_with_Mimpud <- unique(HOG_with_Mimpud$HOG)
  
  HOG_without_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Mimpud <- unique(HOG_without_Mimpud$HOG)
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, HOG_without_Mimpud))
  Mimpud_full_df <- dplyr::left_join(Mimpud_df,Recap_Mimpud %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Mimpud_spe_NODE_genes <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  Mimpud_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, Mimpud_spe_NODE_genes$HOG))
  Mimpud_spe <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  
  Mimpud_full_df <- rbind(Mimpud_spe_NODE_genes, Mimpud_spe)
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Mimpud"] <- "Mimpud_DEOG"
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Species_specific"] <- "Mimpud_specific_OG"
  
  Mimpud_full_df$Prich[Mimpud_full_df$Persist=="Up" & Mimpud_full_df$SignalP==1 & Mimpud_full_df$Pep_length<150 & Mimpud_full_df$P_prop>10] <- "Prorich"
  Mimpud_full_df[is.na(Mimpud_full_df)] <- ""
  
  fwrite(Mimpud_full_df %>% dplyr::filter(Prich=="Prorich"),paste0("./results/IntracellAccomodation/Mimpud_Persist_signalP_Up_speciesspecific_P10prot150.txt"),sep="\t")
  
  
  Mimpud_global <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global$Class <- "Mimpud_WSR"
  Mimpud_global_signalP <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP$Class <- "Mimpud_WSR_signalP"
  Mimpud_Persist_signalP_Prich <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(Prich=="Prorich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_Persist_signalP_Prich$Class <- "Persist_ProRich"
  Mimpud_global_signalP_woPersist <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(Prich=="",SignalP==1) %>%
    dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP_woPersist$Class <- "Mimpud_WSR_signalP_woPersistProRich"
  
  Full_summary <- rbind(Mimpud_global,Mimpud_global_signalP)
  Full_summary <- rbind(Full_summary,Mimpud_Persist_signalP_Prich)
  Full_summary <- rbind(Full_summary,Mimpud_global_signalP_woPersist)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR"])
    Comparaisons <- "Mimpud_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woPersistProRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woPersistProRich"])
    Comparaisons <- "Persist_ProRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woPersistProRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woPersistProRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Mimpud_WSR","Mimpud_WSR_signalP","Mimpud_WSR_signalP_woPersistProRich","Persist_ProRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific_P10prot150.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific_P10prot150.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  # Cysteine rich Mimpud
  # Mimpud Release
  # Specific HOG
  HOG_with_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_")) %>% dplyr::select(HOG)
  HOG_with_Mimpud <- unique(HOG_with_Mimpud$HOG)
  
  HOG_without_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Mimpud <- unique(HOG_without_Mimpud$HOG)
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, HOG_without_Mimpud))
  Mimpud_full_df <- dplyr::left_join(Mimpud_df,Recap_Mimpud %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Mimpud_spe_NODE_genes <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  Mimpud_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, Mimpud_spe_NODE_genes$HOG))
  Mimpud_spe <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  
  Mimpud_full_df <- rbind(Mimpud_spe_NODE_genes, Mimpud_spe)
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Mimpud"] <- "Mimpud_DEOG"
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Species_specific"] <- "Mimpud_specific_OG"
  
  Mimpud_full_df$CysRich[Mimpud_full_df$Release=="Up" & Mimpud_full_df$SignalP==1 & Mimpud_full_df$Pep_length<150 & Mimpud_full_df$C_prop>6] <- "CysRich"
  Mimpud_full_df[is.na(Mimpud_full_df)] <- ""
  
  Mimpud_global <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global$Class <- "Mimpud_WSR"
  Mimpud_global_signalP <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP$Class <- "Mimpud_WSR_signalP"
  Mimpud_Release_signalP_Prich <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="CysRich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_Release_signalP_Prich$Class <- "Release_CysRich"
  Mimpud_global_signalP_woRelease <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="",SignalP==1) %>%
    dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP_woRelease$Class <- "Mimpud_WSR_signalP_woReleaseCysRich"
  
  Full_summary <- rbind(Mimpud_global,Mimpud_global_signalP)
  Full_summary <- rbind(Full_summary,Mimpud_Release_signalP_Prich)
  Full_summary <- rbind(Full_summary,Mimpud_global_signalP_woRelease)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR"])
    Comparaisons <- "Mimpud_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woReleaseCysRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woReleaseCysRich"])
    Comparaisons <- "Release_CysRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woReleaseCysRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woReleaseCysRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Mimpud_WSR","Mimpud_WSR_signalP","Mimpud_WSR_signalP_woReleaseCysRich","Release_CysRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  
  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific_C6prot150.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Release_signalP_Up_speciesspecific_C6prot150.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  
  # Mimpud Persist
  # Specific HOG
  HOG_with_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_")) %>% dplyr::select(HOG)
  HOG_with_Mimpud <- unique(HOG_with_Mimpud$HOG)
  
  HOG_without_Mimpud <- HOG %>% dplyr::filter(str_detect(Gene_id, "^Mimpud_", negate = TRUE)) %>% dplyr::select(HOG)
  HOG_without_Mimpud <- unique(HOG_without_Mimpud$HOG)
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, HOG_without_Mimpud))
  Mimpud_full_df <- dplyr::left_join(Mimpud_df,Recap_Mimpud %>% dplyr::select(Gene_id,HOG,Node_label_Up),by="Gene_id")
  
  Mimpud_spe_NODE_genes <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  Mimpud_spe_NODE_genes$Node_label_Up <- "Species_specific"
  
  Mimpud_spe <- data.frame(HOG=setdiff(HOG_with_Mimpud, Mimpud_spe_NODE_genes$HOG))
  Mimpud_spe <- Mimpud_full_df %>% dplyr::filter(HOG %in% Mimpud_spe$HOG)
  
  Mimpud_full_df <- rbind(Mimpud_spe_NODE_genes, Mimpud_spe)
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Mimpud"] <- "Mimpud_DEOG"
  Mimpud_full_df$Node_label_Up[Mimpud_full_df$Node_label_Up=="Species_specific"] <- "Mimpud_specific_OG"
  
  Mimpud_full_df$CysRich[Mimpud_full_df$Persist=="Up" & Mimpud_full_df$SignalP==1 & Mimpud_full_df$Pep_length<150 & Mimpud_full_df$C_prop>6] <- "CysRich"
  Mimpud_full_df[is.na(Mimpud_full_df)] <- ""
  
  Mimpud_global <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global$Class <- "Mimpud_WSR"
  Mimpud_global_signalP <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(SignalP==1) %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP$Class <- "Mimpud_WSR_signalP"
  Mimpud_Persist_signalP_Prich <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="CysRich") %>% dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_Persist_signalP_Prich$Class <- "Persist_CysRich"
  Mimpud_global_signalP_woPersist <- Mimpud_full_df %>% group_by(Node_label_Up) %>% dplyr::filter(CysRich=="",SignalP==1) %>%
    dplyr::summarise(Overlap=length(Gene_id))
  Mimpud_global_signalP_woPersist$Class <- "Mimpud_WSR_signalP_woPersistCysRich"
  
  Full_summary <- rbind(Mimpud_global,Mimpud_global_signalP)
  Full_summary <- rbind(Full_summary,Mimpud_Persist_signalP_Prich)
  Full_summary <- rbind(Full_summary,Mimpud_global_signalP_woPersist)
  
  Nodes <- unique(Full_summary$Node_label_Up)
  for(Node in Nodes){
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR"])
    Comparaisons <- "Mimpud_WSR_signalP"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woPersistCysRich"] <- sum(Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woPersistCysRich"])
    Comparaisons <- "Persist_CysRich"
    for(i in Comparaisons){
      Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- sum(Full_summary$Overlap[Full_summary$Class==i])
      ClassNode <- Full_summary$Overlap[Full_summary$Class==i & Full_summary$Node_label_Up==Node]
      if(length(ClassNode)==0){ClassNode<-0}
      CtrlNode <- Full_summary$Overlap[Full_summary$Class=="Mimpud_WSR_signalP_woPersistCysRich" & Full_summary$Node_label_Up==Node]
      M <- as.table(rbind(c(ClassNode, CtrlNode), c(Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class==i]-ClassNode,
                                                    Full_summary$Total[Full_summary$Node_label_Up==Node & Full_summary$Class=="Mimpud_WSR_signalP_woPersistCysRich"]-CtrlNode)))
      fishertest <- fisher.test(M)
      Full_summary$FishTest[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$p.value
      Full_summary$oddsratio[Full_summary$Node_label_Up==Node & Full_summary$Class==i] <- fishertest$estimate
    }
  }
  
  Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
  Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
  Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More"
  Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More"
  Full_summary$Class <- factor(Full_summary$Class,levels = c("Mimpud_WSR","Mimpud_WSR_signalP","Mimpud_WSR_signalP_woPersistCysRich","Persist_CysRich"))
  
  p <- ggplot(data=Full_summary, aes(x=Class, y=Overlap, fill=Node_label_Up)) +
    geom_bar(stat="identity", position="fill",  width = 0.65)+
    scale_fill_manual(values = c("#efb800ff","#e65d00ff")) +
    labs(y= "", x = "Intracellular response") +
    geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
    geom_text(data = Full_summary, aes(y = Overlap, label = SignifFisher),size=8,position = position_fill(vjust = 0.45)) +
    theme_classic(base_size = 16)

  fwrite(Full_summary,paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific_C6prot150.txt"),quote=FALSE, sep="\t")
  
  pdf(file=paste0("./results/IntracellAccomodation/Barplot_Mimpud_Persist_signalP_Up_speciesspecific_C6prot150.pdf"),width = 12, height = 8, onefile=FALSE)
  print(p)
  dev.off()
  


