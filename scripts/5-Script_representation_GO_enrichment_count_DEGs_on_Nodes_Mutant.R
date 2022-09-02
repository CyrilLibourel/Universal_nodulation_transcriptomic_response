# Script to count HOGs on species tree
lapply(c("data.table","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","UpSetR","gaston","bio3d","phytools","ggtree","phylobase",
         "reshape2","topGO","purrr","RVenn","ggplot2","Rgraphviz","semtree","tidyverse","rlist","ggplot2","hrbrthemes","gplots",
       "pheatmap","RColorBrewer","ComplexHeatmap","circlize","wesanderson"), require, character.only = TRUE)

# path to input and output folder

setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")

out_path <- "//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3/results/"

# create output folder if necessary
if(dir.exists(out_path)==FALSE){dir.create(out_path,recursive = T)}
# load orthogroup file
HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)

# HOG <- as.data.frame(fread(paste0("//194.199.55.66/evo/commun/projects/nodMimosa/Results_Aug04/Phylogenetic_Hierarchical_Orthogroups/N0_corrected.tsv"),h=T))
# GO_HOG <- as.data.frame(fread(paste0("//194.199.55.66/evo/commun/projects/nodMimosa/Results_Aug04/Phylogenetic_Hierarchical_Orthogroups/N0_goTerms.tsv"),h=T))
# GO_HOG_Mimpud <- dplyr::select(GO_HOG,HOG,Mimpud)
# write.table(GO_HOG_Mimpud,paste0("//194.199.55.66/evo/commun/projects/nodMimosa/Results_Aug04/Phylogenetic_Hierarchical_Orthogroups/Mimpud_goTerms.tsv"), quote=FALSE,row.names=FALSE, col.names=FALSE, sep="\t")
# GO terms Mimpud
geneID2GO <- readMappings(file = paste0("//194.199.55.66/evo/commun/projects/nodMimosa/interproscan/GO/GO_aggregate/Table_Gene_id_GO_aggregate_Mimpud.txt"),sep="\t",IDsep=",")
geneUniverse <- names(geneID2GO)

# Correspondance Mimpud - HOG
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))

# 1st => Up/Down, 2nd => nodes to investigate, 3rd => differences among clones
NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")
for(UP_DOWN in UP_or_DOWN){
        if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",UP_DOWN,".tsv"))

# Load DEGs Mimpud WT and clones

CtaiWT_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
names(CtaiWT_AllUp) <- "Gene_id"
CtaiWT_AllUp <- left_join(CtaiWT_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
CtaiWT_AllUp <- unique(CtaiWT_AllUp$Gene_id[!is.na(CtaiWT_AllUp$HOG)])

CtaiWT_Other <- fread(input = paste0("./Mutant_analysis/Mimpud_Nod_FDR005_logFC1.5_",setdiff(UP_or_DOWN,UP_DOWN),".txt"),header=FALSE)
names(CtaiWT_Other) <- "Gene_id"
CtaiWT_Other <- left_join(CtaiWT_Other,Mimpud_cor_Genes_HOG,by="Gene_id")
CtaiWT_Other <- unique(CtaiWT_Other$Gene_id[!is.na(CtaiWT_Other$HOG)])

GMI1000_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolGMI1000All",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(GMI1000_AllUp) <- "Gene_id"
GMI1000_AllUp <- left_join(GMI1000_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
GMI1000_AllUp <- unique(GMI1000_AllUp$Gene_id[!is.na(GMI1000_AllUp$HOG)])

GMI1000pRalta_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_GMI1000pRaltaAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(GMI1000pRalta_AllUp) <- "Gene_id"
GMI1000pRalta_AllUp <- left_join(GMI1000pRalta_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
GMI1000pRalta_AllUp <- unique(GMI1000pRalta_AllUp$Gene_id[!is.na(GMI1000pRalta_AllUp$HOG)])

hrcV_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_hrcVAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(hrcV_AllUp) <- "Gene_id"
hrcV_AllUp <- left_join(hrcV_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
hrcV_AllUp <- unique(hrcV_AllUp$Gene_id[!is.na(hrcV_AllUp$HOG)])

hrpG_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_hrpGAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(hrpG_AllUp) <- "Gene_id"
hrpG_AllUp <- left_join(hrpG_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
hrpG_AllUp <- unique(hrpG_AllUp$Gene_id[!is.na(hrpG_AllUp$HOG)])

hrpGefpR_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_efpRAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(hrpGefpR_AllUp) <- "Gene_id"
hrpGefpR_AllUp <- left_join(hrpGefpR_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
hrpGefpR_AllUp <- unique(hrpGefpR_AllUp$Gene_id[!is.na(hrpGefpR_AllUp)])

nifH_AllUp <- fread(input = paste0("./Mutant_analysis/Mimpud_nifHAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
names(nifH_AllUp) <- "Gene_id"
nifH_AllUp <- left_join(nifH_AllUp,Mimpud_cor_Genes_HOG,by="Gene_id")
nifH_AllUp <- unique(nifH_AllUp$Gene_id[!is.na(nifH_AllUp$HOG)])

Raw_count <- data.frame(Clone=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"),
                        Overlap=c(length(unique(GMI1000_AllUp)),length(unique(GMI1000pRalta_AllUp)),length(unique(hrcV_AllUp)),
                        length(unique(hrpG_AllUp)),length(unique(hrpGefpR_AllUp)),length(unique(nifH_AllUp)),length(unique(CtaiWT_AllUp))))
Raw_count$Up_Down <- UP_DOWN

write.table(Raw_count,paste0(out_path,"../mutant_results/Table_mutant_genecount_different_",UP_DOWN,".txt"),
            quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")

for(NODE in NODES){
        if(dir.exists(paste0(out_path,"GO_enrichment/",NODE))==FALSE){dir.create(paste0(out_path,"GO_enrichment/",NODE),recursive = T)}
        
        # Ancestral Node to check for
        HOG_id_NODE <- Recap_HOG %>% dplyr::select(HOG,paste0("Node_label_",UP_DOWN),starts_with("Mimpud_"))
        names(HOG_id_NODE) <- c("HOG","Node_label","Mimpud")
        HOG_id_NODE <- data.frame(HOG=unique(Recap_HOG$HOG[Recap_HOG$Node_label==NODE & Recap_HOG$Mimpud==UP_DOWN]))
        
        Mimpud_NODE_genes <- left_join(HOG_id_NODE,Mimpud_cor_Genes_HOG,by="HOG")$Gene_id

        # Overlap CtaiWT_UP_DOWN and NODE
        Mimpud_NODE_genes_Up <- intersect(Mimpud_NODE_genes,CtaiWT_AllUp)
        write.table(Mimpud_NODE_genes_Up,paste0(out_path,"../mutant_results/DEGs_common_",NODE,"_Mimpud_CtaiWT_All",UP_DOWN,".txt"), quote=FALSE,row.names=FALSE, col.names=FALSE, sep="\t")
        
        # Overlap CtaiWT reverse UP_DOWN and NODE
        Mimpud_NODE_genes_other <- intersect(Mimpud_NODE_genes,CtaiWT_Other)
        Mimpud_NODE_genes_other <- setdiff(Mimpud_NODE_genes_other,Mimpud_NODE_genes_Up)
        
        all_df <- data.frame(HOG=Mimpud_NODE_genes_other,Ctai_other=1)
        Mimpud_NODE_df <- data.frame(HOG=Mimpud_NODE_genes_Up,Ctai=1)
        all_df <- full_join(all_df,Mimpud_NODE_df,by="HOG")
        GMI1000_df <- data.frame(GMI1000=1,HOG=GMI1000_AllUp)
        all_df <- left_join(all_df,GMI1000_df,by="HOG")
        GMI1000pRalta_df <- data.frame(pRalta=1,HOG=GMI1000pRalta_AllUp)
        all_df <- left_join(all_df,GMI1000pRalta_df,by="HOG")
        hrcV_df <- data.frame(hrcV=1,HOG=hrcV_AllUp)
        all_df <- left_join(all_df,hrcV_df,by="HOG")
        hrpG_df <- data.frame(hrpG=1,HOG=hrpG_AllUp)
        all_df <- left_join(all_df,hrpG_df,by="HOG")
        hrpGefpR_df <- data.frame(hrpGefpR=1,HOG=hrpGefpR_AllUp)
        all_df <- left_join(all_df,hrpGefpR_df,by="HOG")
        nifH_df <- data.frame(nifH=1,HOG=nifH_AllUp)
        all_df <- left_join(all_df,nifH_df,by="HOG")
        all_df[is.na(all_df)] <- 0
        
        GMI1000 <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$GMI1000==1],Clone = "GMI1000")
        GMI1000_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$GMI1000==0],Clone="GMI1000_no")
        GMI1000_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$GMI1000==1],Clone="GMI1000_bad")
        
        pRalta <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$pRalta==1],Clone = "pRalta")
        pRalta_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$pRalta==0],Clone="pRalta_no")
        pRalta_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$pRalta==1],Clone="pRalta_bad")

        hrcV <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrcV==1],Clone = "hrcV")
        hrcV_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrcV==0],Clone="hrcV_no")
        hrcV_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$hrcV==1],Clone="hrcV_bad")
        
        hrpG <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrpG==1],Clone = "hrpG")
        hrpG_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrpG==0],Clone="hrpG_no")
        hrpG_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$hrpG==1],Clone="hrpG_bad")

        hrpGefpR <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrpGefpR==1],Clone = "hrpGefpR")
        hrpGefpR_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$hrpGefpR==0],Clone="hrpGefpR_no")
        hrpGefpR_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$hrpGefpR==1],Clone="hrpGefpR_bad")

        nifH <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$nifH==1],Clone = "nifH")
        nifH_no_Node_id <- data.frame(HOG=all_df$HOG[all_df$Ctai==1 & all_df$nifH==0],Clone="nifH_no")
        nifH_bad <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1 & all_df$nifH==1],Clone="nifH_bad")

        Ctai <- data.frame(HOG=all_df$HOG[all_df$Ctai==1],Clone = "Ctai")
        Bad_Ctai <- data.frame(HOG=all_df$HOG[all_df$Ctai_other==1],Clone = "Ctai_other")
        WT_length <- nrow(Ctai)
        Bad_length <- nrow(Bad_Ctai)
        
        pdf(file = paste0(out_path,"../mutant_results/Barplot_overlap_",NODE,"_Mimpud_mutants_",UP_DOWN,".pdf"), width=10, height=6)
        barplot(c(nrow(GMI1000)*100/WT_length,nrow(pRalta)*100/WT_length,nrow(hrcV)*100/WT_length,nrow(hrpG)*100/WT_length,nrow(hrpGefpR)*100/WT_length,nrow(nifH)*100/WT_length),
                names.arg=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH"),ylab="Percentage of genes recruited by each strain",ylim=c(0,100))
        dev.off()
        
        pdf(file = paste0(out_path,"../mutant_results/Barplot_bad_overlap_",NODE,"_Mimpud_mutants_",UP_DOWN,".pdf"), width=10, height=6)
        barplot(c(nrow(GMI1000_bad)*100/WT_length,nrow(pRalta_bad)*100/Bad_length,nrow(hrcV_bad)*100/Bad_length,nrow(hrpG_bad)*100/Bad_length,nrow(hrpGefpR_bad)*100/Bad_length,nrow(nifH_bad)*100/Bad_length),
                names.arg=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH"),ylab="Percentage of genes recruited in wrong way by each strain",ylim=c(0,100))
        dev.off()
        
        listInput <-list(GMI1000=GMI1000$HOG,pRalta=pRalta$HOG,
                         hrcV=hrcV$HOG,
                         hrpG=hrpG$HOG,
                         hrpGefpR=hrpGefpR$HOG,
                         nifH=nifH$HOG,Ctai=Ctai$HOG)
        
        pdf(file = paste0(out_path,"../mutant_results/UpSet_",NODE,"_Mimpud_mutants_",UP_DOWN,".pdf"), width=14, height=8)
        upset(fromList(listInput), keep.order=TRUE, sets=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"),order.by="degree",nintersects=NA,decreasing = FALSE)
        dev.off()
        
        HOG_all_df <- rbind(GMI1000,pRalta,hrcV,hrpG,hrpGefpR,nifH,Ctai)
        write.table(HOG_all_df,paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_",UP_DOWN,".txt"), quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")
        
        HOG_all_df <- rbind(HOG_all_df,GMI1000_no_Node_id,
                            pRalta_no_Node_id,hrcV_no_Node_id,hrpG_no_Node_id,hrpGefpR_no_Node_id,nifH_no_Node_id,
                            GMI1000_bad,pRalta_bad,hrcV_bad,hrpG_bad,hrpGefpR_bad,nifH_bad)
        
        write.table(HOG_all_df,paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_",UP_DOWN,"_with_bad.txt"), quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")
        
        CLONES <- unique(HOG_all_df$Clone)
        # for(clone in CLONES){
        #         geneList <- factor(as.integer(geneUniverse %in% HOG_all_df$HOG[HOG_all_df$Clone==clone]))
        #         names(geneList) <- geneUniverse
        #         if(length(levels(geneList))>1){
        #                 # Biological process
        #                 my_BP_data <- new("topGOdata", description=paste0(clone), ontology="BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        #                 #my_BP_data
        #                 TOP_NODES=length(my_BP_data@graph@nodes)
        #                 resultFisher <- runTest(my_BP_data, algorithm="classic", statistic="fisher")
        #                 resultFisher.weight01 <- runTest(my_BP_data, algorithm='weight01', statistic="fisher")
        #                 #resultttest.weight01 <- runTest(my_BP_data, algorithm = "weight01", statistic = "t")
        #                 allRes_BP <- GenTable(my_BP_data, classicFisher = resultFisher, weight01Fisher = resultFisher.weight01,
        #                                       orderBy = "weight01Fisher", ranksOf = "weight01Fisher", topNodes = TOP_NODES)
        #                 write.table(allRes_BP, file=paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        # 
        #                 # NodesNumber <- allRes_BP[as.numeric(allRes_BP$weight01Fisher) < 0.01,]
        #                 # pdf(file = paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".pdf"), width = 8, height = 8)
        #                 # showSigOfNodes(my_BP_data, score(resultFisher.weight01), firstSigNodes = nrow(NodesNumber), useInfo = 'all')
        #                 # dev.off()
        # 
        #                 write.table(allRes_BP[as.numeric(allRes_BP$weight01Fisher)<0.001,], file=paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".signif"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        # 
        #                 # Molecular function
        #                 my_MF_data <- new("topGOdata", description=paste0(clone), ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        #                 #my_MF_data
        #                 TOP_NODES=length(my_MF_data@graph@nodes)
        #                 resultFisher <- runTest(my_MF_data, algorithm="classic", statistic="fisher")
        #                 resultFisher.weight01 <- runTest(my_MF_data, algorithm='weight01', statistic="fisher")
        #                 #resultttest.weight01 <- runTest(my_MF_data, algorithm = "weight01", statistic = "t")
        # 
        #                 allRes_MF <- GenTable(my_MF_data, classicFisher = resultFisher, weight01Fisher = resultFisher.weight01,
        #                                       orderBy = "weight01Fisher", ranksOf = "weight01Fisher", topNodes = TOP_NODES)
        #                 write.table(allRes_MF, file=paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_MolFunc_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        # 
        #                 # NodesNumber <- allRes_MF[as.numeric(allRes_MF$weight01Fisher) < 0.01,]
        #                 # pdf(file = paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_MolFunc_",UP_DOWN,".pdf"), width = 8, height = 8)
        #                 # showSigOfNodes(my_MF_data, score(resultFisher.weight01), firstSigNodes = nrow(NodesNumber), useInfo = 'all')
        #                 # dev.off()
        # 
        #                 write.table(allRes_MF[as.numeric(allRes_MF$weight01Fisher)<0.001,], file=paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_MolFunc_",UP_DOWN,".signif"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        # 
        #         }
        #         else{}
        # }
        
}
# End UP_Down
}

#### GO enrichment GMI1000 induced genes but not in Ctai ####


geneID2GO <- readMappings(file = paste0("//194.199.55.66/evo/commun/projects/nodMimosa/interproscan/GO/GO_aggregate/Table_Gene_id_GO_aggregate_Mimpud.txt"),sep="\t",IDsep=",")
geneUniverse <- names(geneID2GO)

# Correspondance Mimpud - HOG
Mimpud_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Mimpud_cor_Genes_HOG <- Mimpud_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))

# 1st => Up/Down, 2nd => nodes to investigate, 3rd => differences among clones
NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")
for(UP_DOWN in UP_or_DOWN){
  if(dir.exists(paste0(out_path,"GO_enrichment"))==FALSE){dir.create(paste0(out_path,"GO_enrichment"),recursive = T)}
  Recap_HOG <- fread(paste0("./results/N0_corrected_nodMimosa_project_with_Medtru_annot_and_nsgenes_ranking_Ancestral_nodes_with_sampling_",UP_DOWN,".tsv"))
  
  # Load DEGs Mimpud WT and clones
  CtaiWT <- fread(input = paste0("./Mutant_analysis/Mimpud_Nod_FDR005_logFC1.5_",UP_DOWN,".txt"),header=FALSE)
  names(CtaiWT) <- "Gene_id"
  CtaiWT <- dplyr::left_join(CtaiWT,Mimpud_cor_Genes_HOG,by="Gene_id")
  CtaiWT <- data.frame(Gene_id=unique(CtaiWT$Gene_id[!is.na(CtaiWT$HOG)]),Clone="Ctai")

  GMI1000 <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolGMI1000All",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(GMI1000) <- "Gene_id"
  GMI1000 <- dplyr::left_join(GMI1000,Mimpud_cor_Genes_HOG,by="Gene_id")
  GMI1000 <- data.frame(Gene_id=unique(GMI1000$Gene_id[!is.na(GMI1000$HOG)]),Clone="GMI1000")
  
  pRalta <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolGMI1000pRaltaAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(pRalta) <- "Gene_id"
  pRalta <- dplyr::left_join(pRalta,Mimpud_cor_Genes_HOG,by="Gene_id")
  pRalta <- data.frame(Gene_id=unique(pRalta$Gene_id[!is.na(pRalta$HOG)]),Clone="pRalta")
  
  hrcV <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolhrcVAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(hrcV) <- "Gene_id"
  hrcV <- dplyr::left_join(hrcV,Mimpud_cor_Genes_HOG,by="Gene_id")
  hrcV <- data.frame(Gene_id=unique(hrcV$Gene_id[!is.na(hrcV$HOG)]),Clone="hrcV")
  
  
  hrpG <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolhrpGAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(hrpG) <- "Gene_id"
  hrpG <- dplyr::left_join(hrpG,Mimpud_cor_Genes_HOG,by="Gene_id")
  hrpG <- data.frame(Gene_id=unique(hrpG$Gene_id[!is.na(hrpG$HOG)]),Clone="hrpG")
  
  hrpGefpR <- fread(input = paste0("./Mutant_analysis/Mimpud_RsolhrpGefpRAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(hrpGefpR) <- "Gene_id"
  hrpGefpR <- dplyr::left_join(hrpGefpR,Mimpud_cor_Genes_HOG,by="Gene_id")
  hrpGefpR <- data.frame(Gene_id=unique(hrpGefpR$Gene_id[!is.na(hrpGefpR$HOG)]),Clone="hrpGefpR")
  
  nifH <- fread(input = paste0("./Mutant_analysis/Mimpud_CtainifHAll",UP_DOWN,"_FDR005_logFC1.5.txt"),header=FALSE)
  names(nifH) <- "Gene_id"
  nifH <- dplyr::left_join(nifH,Mimpud_cor_Genes_HOG,by="Gene_id")
  nifH <- data.frame(Gene_id=unique(nifH$Gene_id[!is.na(nifH$HOG)]),Clone="nifH")
  
  All_data <- rbind(CtaiWT,GMI1000,pRalta,hrcV,hrpG,hrpGefpR,nifH)
  clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
  for(CLONE in clones){
  if(dir.exists(paste0(out_path,"GO_enrichment/"))==FALSE){dir.create(paste0(out_path,"GO_enrichment/"),recursive = T)}
  clone_specific <- unique(setdiff(All_data$Gene_id[All_data$Clone==CLONE],All_data$Gene_id[All_data$Clone=="Ctai"]))
  geneList <- factor(as.integer(geneUniverse %in% clone_specific))
            names(geneList) <- geneUniverse
            if(length(levels(geneList))>1){
                    # Biological process
                    my_BP_data <- new("topGOdata", description=paste0("GMI1000"), ontology="BP", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
                    #my_BP_data
                    TOP_NODES=length(my_BP_data@graph@nodes)
                    resultFisher <- runTest(my_BP_data, algorithm="classic", statistic="fisher")
                    resultFisher.weight01 <- runTest(my_BP_data, algorithm='weight01', statistic="fisher")
                    #resultttest.weight01 <- runTest(my_BP_data, algorithm = "weight01", statistic = "t")
                    allRes_BP <- GenTable(my_BP_data, classicFisher = resultFisher, weight01Fisher = resultFisher.weight01,
                                          orderBy = "weight01Fisher", ranksOf = "weight01Fisher", topNodes = TOP_NODES)
                    write.table(allRes_BP, file=paste0(out_path,"GO_enrichment/",CLONE,"_specific_response_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                    allRes_BP$weight01Fisher[allRes_BP$weight01Fisher=="< 1e-30"] <- "1e-30"
                    write.table(allRes_BP[as.numeric(allRes_BP$weight01Fisher)<0.01,], file=paste0(out_path,"GO_enrichment/",CLONE,"_specific_response_GO_enrichment_BioProcess_",UP_DOWN,".signif"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

                    # Molecular function
                    my_MF_data <- new("topGOdata", description=paste0("GMI1000"), ontology="MF", allGenes=geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
                    #my_MF_data
                    TOP_NODES=length(my_MF_data@graph@nodes)
                    resultFisher <- runTest(my_MF_data, algorithm="classic", statistic="fisher")
                    resultFisher.weight01 <- runTest(my_MF_data, algorithm='weight01', statistic="fisher")
                    #resultttest.weight01 <- runTest(my_MF_data, algorithm = "weight01", statistic = "t")

                    allRes_MF <- GenTable(my_MF_data, classicFisher = resultFisher, weight01Fisher = resultFisher.weight01,
                                          orderBy = "weight01Fisher", ranksOf = "weight01Fisher", topNodes = TOP_NODES)
                    write.table(allRes_MF, file=paste0(out_path,"GO_enrichment/",CLONE,"_specific_response_GO_enrichment_MolFunc_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                    allRes_MF$weight01Fisher[allRes_MF$weight01Fisher=="< 1e-30"] <- "1e-30"
                    write.table(allRes_MF[as.numeric(allRes_MF$weight01Fisher)<0.01,], file=paste0(out_path,"GO_enrichment/",CLONE,"_specific_response_GO_enrichment_MolFunc_",UP_DOWN,".signif"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

            }else{}
}
  # End UP_Down
}


#### Merge Up and Down ####

Overlap_Up <- data.frame()
for(NODE in NODES){
  df <- fread(paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_Up.txt"))
  df <- df %>% group_by(Clone) %>% dplyr::summarise(Overlap=length(HOG))
  df$HOG_Node <- NODE; df$Overlap <- (df$Overlap)
  Overlap_Up <- rbind(Overlap_Up,df)
}


Overlap_Down <- data.frame()
for(NODE in NODES){
  df <- fread(paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_Down.txt"))
  df <- df %>% group_by(Clone) %>% dplyr::summarise(Overlap=length(HOG))
  df$HOG_Node <- NODE; df$Overlap <- (df$Overlap)
  Overlap_Down <- rbind(Overlap_Down,df)
}



Overlap_Up$Up_Down <- "Up"
Overlap_Down$Up_Down <- "Down"
Overlap_Down$Overlap <- Overlap_Down$Overlap*(-1)

Overlap <- rbind(Overlap_Up,Overlap_Down)
Overlap$Mutant <- factor(Overlap$Clone, levels=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"))
Overlap$HOG_Node <- factor(Overlap$HOG_Node, levels=c("N1","N2","Mimpud"))
Overlap$Up_Down <- factor(Overlap$Up_Down, levels=c("Up","Down"))

LABEL_NODE <- c(`N1` = "NFN node",
                `N2` = "Fabales node",
                `Mimpud` = "Mimosa pudica")
mycol<-c("coral2","deepskyblue3")
pdf(file = paste0(out_path,"../mutant_results/Barplot_overlap_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
ggplot(Overlap, aes(x=Mutant, y=Overlap, fill=Up_Down)) + 
  geom_bar(stat="identity", position="identity") + facet_wrap(~ HOG_Node, nrow=1, labeller = as_labeller(LABEL_NODE)) + 
  scale_fill_manual(values=mycol) +
  labs(y= "Genes shared with C.tai in the different nodes (in % of C.tai)", x = "Strains") + theme_classic(base_size = 16)
dev.off()


########## Proportions of genes in each nodes 

Overlap_Up <- data.frame()
for(NODE in NODES){
  df <- fread(paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_Up_with_bad.txt"))
  df <- df %>% group_by(Clone) %>% dplyr::summarise(Overlap=length(HOG))
  df$HOG_Node <- NODE
  # df$Overlap <- (df$Overlap*100)/df$Overlap[df$Clone=="Ctai"]
  Overlap_Up <- rbind(Overlap_Up,df)
}
# Overlap_Up$Overlap <- (Overlap_Up$Overlap*100)/sum(Overlap_Up$Overlap[Overlap_Up$Clone=="Ctai"])


Overlap_Down <- data.frame()
for(NODE in NODES){
  df <- fread(paste0(out_path,"../mutant_results/Table_overlap_mutant_different_HOGs_",NODE,"_Down_with_bad.txt"))
  df <- df %>% group_by(Clone) %>% dplyr::summarise(Overlap=length(HOG))
  df$HOG_Node <- NODE
  # df$Overlap <- (df$Overlap*100)/df$Overlap[df$Clone=="Ctai"]
  Overlap_Down <- rbind(Overlap_Down,df)
}

# Overlap_Down$Overlap <- (Overlap_Down$Overlap*100)/sum(Overlap_Down$Overlap[Overlap_Down$Clone=="Ctai"])

Overlap_Up$Up_Down <- "Up"

Full_summary <- Overlap_Up
Clones <- c("Ctai","GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
Clones_sum <- data.frame(Clones=Clones,Sum=0)
for(i in Clones){
  Clones_sum$Sum[Clones_sum$Clones==i] <- sum(Full_summary$Overlap[Full_summary$Clone==i])
  
}

Nodes <- c("N1","N2","Mimpud")
for(Node in Nodes){
  for(i in Clones[-1]){
    CloneNode <- Full_summary$Overlap[Full_summary$Clone==i & Full_summary$HOG_Node==Node]
    CtaiNode <- Full_summary$Overlap[Full_summary$Clone=="Ctai" & Full_summary$HOG_Node==Node]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Clones==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Clones=="Ctai"]-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- test$p.value
    Full_summary$Total[Full_summary$HOG_Node=="N1" & Full_summary$Clone==i] <- Clones_sum$Sum[Clones_sum$Clones==i]
    Full_summary$FishTest[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$HOG_Node=="N1" & Full_summary$Clone==i] <- Clones_sum$Sum[Clones_sum$Clones==i]
  }
}

Full_summary$Total[Full_summary$Clone=="Ctai" & Full_summary$HOG_Node=="N1"] <- Clones_sum$Sum[Clones_sum$Clones=="Ctai"]

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"

Full_summary$HOG_Node[Full_summary$HOG_Node=="Mimpud"] <- "M. pudica specific"
Full_summary$HOG_Node <- factor(Full_summary$HOG_Node,levels = c("N1","N2","M. pudica specific"))
Full_summary$Clone <- factor(Full_summary$Clone,levels = c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"))
Full_summary_Up <- Full_summary

Overlap_Down$Up_Down <- "Down"

Full_summary <- Overlap_Down
Clones <- c("Ctai","GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
Clones_sum <- data.frame(Clones=Clones,Sum=0)
for(i in Clones){
  Clones_sum$Sum[Clones_sum$Clones==i] <- sum(Full_summary$Overlap[Full_summary$Clone==i])
  
}

Nodes <- c("N1","N2","Mimpud")
for(Node in Nodes){
  for(i in Clones[-1]){
    CloneNode <- Full_summary$Overlap[Full_summary$Clone==i & Full_summary$HOG_Node==Node]
    CtaiNode <- Full_summary$Overlap[Full_summary$Clone=="Ctai" & Full_summary$HOG_Node==Node]
    M <- as.table(rbind(c(CloneNode, CtaiNode), c(Clones_sum$Sum[Clones_sum$Clones==i]-CloneNode,
                                                  Clones_sum$Sum[Clones_sum$Clones=="Ctai"]-CtaiNode)))
    test <- chisq.test(M)
    fishertest <- fisher.test(M)
    Full_summary$Chi2[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- test$p.value
    Full_summary$Total[Full_summary$HOG_Node=="N1" & Full_summary$Clone==i] <- Clones_sum$Sum[Clones_sum$Clones==i]
    Full_summary$FishTest[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- fishertest$p.value
    Full_summary$oddsratio[Full_summary$HOG_Node==Node & Full_summary$Clone==i] <- fishertest$estimate
    Full_summary$Total[Full_summary$HOG_Node=="N1" & Full_summary$Clone==i] <- Clones_sum$Sum[Clones_sum$Clones==i]
  }
}

Full_summary$Total[Full_summary$Clone=="Ctai" & Full_summary$HOG_Node=="N1"] <- Clones_sum$Sum[Clones_sum$Clones=="Ctai"]

Full_summary$Signif[Full_summary$Chi2<0.05] <- "*"
Full_summary$Signif[Full_summary$Chi2<0.01] <- "**"
Full_summary$Signif[Full_summary$Chi2<0.001] <- "***"
Full_summary$SignifFisher[Full_summary$FishTest<0.05] <- "*"
Full_summary$SignifFisher[Full_summary$FishTest<0.01] <- "**"
Full_summary$SignifFisher[Full_summary$FishTest<0.001] <- "***"

Full_summary$HOG_Node[Full_summary$HOG_Node=="Mimpud"] <- "M. pudica specific"
Full_summary$HOG_Node <- factor(Full_summary$HOG_Node,levels = c("N1","N2","M. pudica specific"))
Full_summary$Clone <- factor(Full_summary$Clone,levels = c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"))
Full_summary_Down <- Full_summary

Full_summary <- rbind(Full_summary_Up,Full_summary_Down)
Full_summary <- dplyr::filter(Full_summary,Clone %in% all_of(Clones))

#CLASS <- data.frame(class=Full_summary$Clone)
#CLASS <- tidyr::separate(CLASS,class,into=c("Strain","class"),sep = "_")
#CLASS$class[is.na(CLASS$class)] <- "DEGs"

Full_summary$Clone <- factor(Full_summary$Clone, levels=c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH","Ctai"))
Full_summary$HOG_Node <- factor(Full_summary$HOG_Node, levels=c("N1","N2","M. pudica specific"))
#Full_summary_Down$class <- factor(CLASS$class, levels=c("DEGs","no","bad"))
Full_summary$Up_Down <- factor(Full_summary$Up_Down, levels=c("Up","Down"))
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="*"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="**"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio<1 & Full_summary$SignifFisher=="***"] <- "Less than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="*"] <- "More than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="**"] <- "More than C.tai"
Full_summary$Odds[Full_summary$oddsratio>1 & Full_summary$SignifFisher=="***"] <- "More than C.tai"
Full_summary$HOG_Node <- factor(Full_summary$HOG_Node, levels=c("N1","N2","M. pudica specific"))

LABEL_NODE <- c(`Up` = "Up",
                `Down` = "Down")

Full_summary <- Full_summary %>% tidyr::unite("Clone_Up_Down",c("Clone","Up_Down"), sep="_", remove=FALSE)

Raw_count_Up <- fread(paste0(out_path,"../mutant_results/Table_mutant_genecount_different_Up.txt"))
Raw_count_Down <- fread(paste0(out_path,"../mutant_results/Table_mutant_genecount_different_Down.txt"))
Raw_count <- rbind(Raw_count_Up,Raw_count_Down)
Raw_count <- Raw_count %>% tidyr::unite("Clone_Up_Down",c("Clone","Up_Down"), sep="_", remove=FALSE)
Raw_count <- dplyr::select(Raw_count,-Clone,-Up_Down)

Full_summary <- left_join(Full_summary,Raw_count,by="Clone_Up_Down")

mycol<-wes_palette("Zissou1", 5)[c(1,2,4,5)]


pdf(file = paste0(out_path,"../mutant_results/Barplot_count_shared_with_Ctai_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
ggplot(data=Full_summary, aes(x=Clone, y=Overlap.x, fill=HOG_Node)) +
  geom_bar(stat="identity", position="stack",  width = 0.65) +
  scale_fill_manual(values = mycol) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
  labs(y= "Number of genes shared with C.tai", x = "Strains") +
  geom_text(data = Full_summary, aes(y = max(Total*1.05,na.rm=TRUE), label = Total)) +
  theme_classic(base_size = 16) +
  geom_point(aes(x=Clone, y=Overlap.y, col="red")) +
  geom_text(data = Full_summary, aes(y = Overlap.y, label = Overlap.y))
dev.off()

pdf(file = paste0(out_path,"../mutant_results/Barplot_Proportion_shared_with_Ctai_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
ggplot(data=Full_summary, aes(x=Clone, y=Overlap.x, fill=HOG_Node)) +
  geom_bar(stat="identity", position="fill",  width = 0.65) +
  scale_fill_manual(values = mycol,labels = c("NFN node", "Fabales node", "M. pudica")) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
  labs(y= "Number of genes shared with C.tai", x = "Strains") +
  geom_text(data = Full_summary, aes(y = 1.05, label = Total)) +
  geom_text(data = Full_summary, aes(y = Overlap.x, label = SignifFisher),size=8,position = position_fill(vjust = 0.5))+
  scale_colour_manual(values=c("coral2", "black")) +
  theme_classic(base_size = 16) + labs(fill = "Node") 
dev.off()


# 
# pdf(file = paste0(out_path,"../mutant_results/Barplot_count_uninduced_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
# 
# ggplot(data=Overlap%>%dplyr::filter(class=="no"), aes(x=Mutant, y=Overlap, fill=HOG_Node)) +
#     geom_bar(stat="identity",  width = 0.65) +
#     scale_fill_manual(values = mycol) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
#     labs(y= "Number of uninduced genes compared to C.tai", x = "Strains")  + theme_classic(base_size = 16)
# dev.off()
# 
# pdf(file = paste0(out_path,"../mutant_results/Barplot_count_opposite_to_Ctai_multi_Nodes_Mimpud_mutants_Up_Down.pdf"), width=15, height=10)
# ggplot(data=Overlap%>%dplyr::filter(class=="bad"), aes(x=Mutant, y=Overlap, fill=HOG_Node)) +
#     geom_bar(stat="identity",  width = 0.65) +
#     scale_fill_manual(values = mycol) + facet_wrap(~ Up_Down, labeller = as_labeller(LABEL_NODE), nrow=1) + 
#     labs(y= "Number of genes shared in the opposite way of C.tai", x = "Strains")  + theme_classic(base_size = 16)
# dev.off()
# 


#### Concatenation data BioProcess all clones #### 
# Biological processes

NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up")
clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
for(UP_DOWN in UP_or_DOWN){
        for(NODE in NODES){
                WT_Node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                #WT_Node$weight01Fisher <- -log10(WT_Node$weight01Fisher)
                names(WT_Node) <- c("GO.ID","Term","Ctai")
                for(clone in clones){
                        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
                        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
                        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
                        names(Clone_node) <- c("GO.ID",clone)
                        WT_Node <- dplyr::left_join(WT_Node,Clone_node, by="GO.ID")
                        
                }
                WT_Node[is.na(WT_Node)] <- 1
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,paste0(clones),Ctai)
                write.table(WT_Node, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,"_allGO.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
        }
}




NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Down")
clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
for(UP_DOWN in UP_or_DOWN){
        full_data <- data.frame()
        for(NODE in NODES){
                WT_Node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                #WT_Node$weight01Fisher <- -log10(WT_Node$weight01Fisher)
                names(WT_Node) <- c("GO.ID","Term","Ctai")
                for(clone in clones){
                        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
                        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
                        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
                        names(Clone_node) <- c("GO.ID",clone)
                        WT_Node <- dplyr::left_join(WT_Node,Clone_node, by="GO.ID")
                        
                }
                WT_Node[is.na(WT_Node)] <- 1
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,paste0(clones),Ctai)
                write.table(WT_Node, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,"_allGO.txt.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                
        }
}


#### Heatmap GO enrichment Mimpud ####
# Biological processes

NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up")
clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
for(UP_DOWN in UP_or_DOWN){
        for(NODE in NODES){
                WT_Node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                names(WT_Node) <- c("GO.ID","Term","Ctai")
        for(clone in clones){
        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
        names(Clone_node) <- c("GO.ID",clone)
        WT_Node <- dplyr::left_join(WT_Node,Clone_node, by="GO.ID")

        }
                
                WT_Node[is.na(WT_Node)] <- 1
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,paste0(clones),Ctai)
                write.table(WT_Node, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                
                WT_Node <- WT_Node %>% rowwise() %>% dplyr::mutate(MIN_pval=min(paste0(clones),Ctai))
                WT_Node <- dplyr::filter(WT_Node,MIN_pval<0.01)
                dat <- as.matrix(dplyr::select(WT_Node,paste0(clones),Ctai))
                row.names(dat) <- WT_Node$Term

                col_fun_up = colorRamp2(c(0.01, 0.001, 0.000001), c("white", "chocolate2", "firebrick2"))
                ht_up <- Heatmap(dat, name = "Weighted Fisher p-val", cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE,column_names_rot = 90, column_names_side = "top",
                              row_names_side = "left", col = col_fun_up , width = unit(8, "cm"), height = unit(16, "cm"))
                
                pdf(paste0(out_path,"GO_enrichment/Heatmap_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".pdf"),height=14,width=14)
                draw(ht_up)
                dev.off()
               
        }
}




NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Down")
clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")
for(UP_DOWN in UP_or_DOWN){
        for(NODE in NODES){
                WT_Node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                names(WT_Node) <- c("GO.ID","Term","Ctai")
                for(clone in clones){
                        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
                        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
                        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
                        names(Clone_node) <- c("GO.ID",clone)
                        WT_Node <- dplyr::left_join(WT_Node,Clone_node, by="GO.ID")
                        
                }
                
                WT_Node[is.na(WT_Node)] <- 1
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,paste0(clones),Ctai)
                write.table(WT_Node, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                
                WT_Node <- WT_Node %>% rowwise() %>% dplyr::mutate(MIN_pval=min(paste0(clones),Ctai))
                WT_Node <- dplyr::filter(WT_Node,MIN_pval<0.01)
                dat <- as.matrix(dplyr::select(WT_Node,paste0(clones),Ctai))
                row.names(dat) <- WT_Node$Term
                
                col_fun_down = colorRamp2(c(0.01, 0.001, 0.000001), c("white", "royalblue1", "blue4"))
                ht_down <- Heatmap(dat, name = "Weighted Fisher p-val", cluster_rows = TRUE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE,column_names_rot = 90, column_names_side = "top",
                                 row_names_side = "left", col = col_fun_down , width = unit(8, "cm"), height = unit(16, "cm"))
                
                pdf(paste0(out_path,"GO_enrichment/Heatmap_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".pdf"),height=14,width=14)
                draw(ht_down)
                dev.off()
                
        }
}


NODES <- c("N1","N2","Mimpud")
UP_or_DOWN <- c("Up","Down")
clones <- c("GMI1000","pRalta","hrcV","hrpG","hrpGefpR","nifH")

        for(NODE in NODES){
        UP_DOWN="Up"
                WT_Node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                names(WT_Node) <- c("GO.ID","Term","Ctai_Up")
                for(clone in clones){
                        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
                        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
                        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
                        names(Clone_node) <- c("GO.ID",paste0(clone,"_Up"))
                        WT_Node <- dplyr::left_join(WT_Node,Clone_node, by="GO.ID")
                        
                }
                
                WT_Node[is.na(WT_Node)] <- 1
                WT_Node <- dplyr::select(WT_Node,GO.ID,Term,paste0(clones,"_Up"),Ctai_Up)
                write.table(WT_Node, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                
                #WT_Node <- WT_Node %>% rowwise() %>% dplyr::mutate(MIN_pval=min(paste0(clones),Ctai))
                #WT_Node <- dplyr::filter(WT_Node,MIN_pval<0.01)

        UP_DOWN="Down"
                WT_Node_Down <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_Ctai_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                WT_Node_Down <- dplyr::select(WT_Node_Down,GO.ID,Term,weight01Fisher)
                #WT_Node <- dplyr::filter(WT_Node,weight01Fisher<0.01)
                names(WT_Node_Down) <- c("GO.ID","Term","Ctai_Down")
                for(clone in clones){
                        Clone_node <- fread(paste0(out_path,"GO_enrichment/",NODE,"/",NODE,"_",clone,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"))
                        Clone_node <- dplyr::select(Clone_node,GO.ID,weight01Fisher)
                        #Clone_node <- dplyr::filter(Clone_node,weight01Fisher<0.01)
                        #Clone_node$weight01Fisher <- -log10(Clone_node$weight01Fisher)
                        names(Clone_node) <- c("GO.ID",paste0(clone,"_Down"))
                        WT_Node_Down <- dplyr::left_join(WT_Node_Down,Clone_node, by="GO.ID")
                        
                }
                
                WT_Node_Down[is.na(WT_Node_Down)] <- 1
                WT_Node_Down <- dplyr::select(WT_Node_Down,GO.ID,paste0(clones,"_Down"),Ctai_Down)
                write.table(WT_Node_Down, file=paste0(out_path,"GO_enrichment/Summary_",NODE,"_GO_enrichment_BioProcess_",UP_DOWN,".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
                
                WT_Up_Down <- dplyr::left_join(WT_Node,WT_Node_Down, by="GO.ID")
                WT_Up_Down <- WT_Up_Down %>% rowwise() %>% dplyr::mutate(MIN_pval=min(paste0(clones,"_Up"),Ctai_Up,paste0(clones,"_Down"),Ctai_Down))
                WT_Up_Down <- dplyr::filter(WT_Up_Down,as.numeric(MIN_pval)<0.01)
                
                dat_Up <- as.matrix(dplyr::select(WT_Up_Down,paste0(clones,"_Up"),Ctai_Up))
                row.names(dat_Up) <- WT_Up_Down$Term
                col_fun_up = colorRamp2(c(0.05, 0.001, 0.0001), c("grey90", "chocolate2", "firebrick2"))
                ht_up <- Heatmap(dat_Up, name = "Weighted Fisher p-val Up", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE,column_names_rot = 90, column_names_side = "top",
                                 column_title = "Up regulated genes",column_labels = c(paste0(clones),"C.tai"),row_names_side = "left", col = col_fun_up ,
                                 width = unit(8, "cm"), height = unit(16, "cm"))
                
                dat_Down <- as.matrix(dplyr::select(WT_Up_Down,paste0(clones,"_Down"),Ctai_Down))
                row.names(dat_Down) <- WT_Up_Down$Term
                col_fun_down = colorRamp2(c(0.05, 0.001, 0.0001), c("grey90", "royalblue1", "blue4"))
                ht_down <- Heatmap(dat_Down, name = "Weighted Fisher p-val Down", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, show_row_dend = FALSE,column_names_rot = 90, column_names_side = "top",
                                   column_title = "Down regulated genes",column_labels = c(paste0(clones),"C.tai") ,row_names_side = "left", col = col_fun_down ,
                                   width = unit(8, "cm"), height = unit(16, "cm"))
                
                
                pdf(paste0(out_path,"GO_enrichment/Heatmap_",NODE,"_GO_enrichment_BioProcess_Up_and_Down.pdf"),height=10,width=14)
                draw(ht_up+ht_down)
                dev.off()
                

}



