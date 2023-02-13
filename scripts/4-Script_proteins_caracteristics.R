# Script to count HOGs on species tree
lapply(c("data.table","dplyr","stringr","ape","seqinr","adephylo","tidyverse","tidyr","VennDiagram","DBI","grid","mixOmics","DECIPHER","Biostrings","gridExtra",
         "UpSetR","gaston","bio3d","phytools","phylobase","ggtree","ggplot2","castor","ComplexHeatmap","circlize","dbplyr","wesanderson"), library, character.only = TRUE)


setwd("//194.199.55.66/evo/commun/projects/nodMimosa/Analysis_v3")
sqlite_db_out <- "C:/Users/cyril.libourel/OneDrive/mimosaProject/work/Sql_db.sqlite"
outcon <- dbConnect(RSQLite::SQLite(), dbname=sqlite_db_out)
src_dbi(outcon)

# load orthogroup file
HOG <- dbGetQuery(outcon, paste0("select * from ","N0_corrected_Gene_id_HOGs_correspondance"))

##### General proteins information ####
cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
full_pep <- Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve")

Species <- c("Medtru","Mimpud","Aeseve","Arahyp","Datglo","Lotjap","Parand","Lupalb","Hiprha")
FullSignalP <- fread("../ncr/output.gff3",h=F)

for(sp in Species){
  SignalP <- FullSignalP %>% filter(str_detect(V1, paste0("^",sp,"_"))) %>% dplyr::select(V1,V3,V5)
  names(SignalP) <- c("Gene_id","SignalP","endSP")
  SignalP$Gene_id <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",SignalP$Gene_id)
  
  gene_id <- data.frame(Gene_id=names(full_pep)) %>% dplyr::filter(str_detect(Gene_id, paste0("^",sp,"_")))
  
  pep <- full_pep[names(full_pep) %in% gene_id$Gene_id]
  # pep <- full_pep[names(full_pep) %in% SignalP$Gene_id]
  pep_cleaned <- pep
  pep_df <- data.frame(Gene_id=names(pep_cleaned))
  pb <- txtProgressBar(min = 0, max = length(gene_id$Gene_id), style = 3)
  all_aacomp <- data.frame()
  # remove signalP from sequence
  for(i in gene_id$Gene_id){
    subset <- SignalP %>% filter(Gene_id==i)
    if(nrow(subset)>0){pep_cleaned[names(pep_cleaned) %in% i]<- subseq(pep_cleaned[names(pep_cleaned) %in% i], start = subset$endSP+1)
    }else{}
    pep_df$Pep_length[pep_df$Gene_id==i]=width(pep_cleaned[names(pep_cleaned) %in% i])
    pep_df$A_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("A"), as.prob=T))*100
    pep_df$C_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("C"), as.prob=T))*100
    pep_df$D_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("D"), as.prob=T))*100
    pep_df$E_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("E"), as.prob=T))*100
    pep_df$F_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("F"), as.prob=T))*100
    pep_df$G_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("G"), as.prob=T))*100
    pep_df$H_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("H"), as.prob=T))*100
    pep_df$I_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("I"), as.prob=T))*100
    pep_df$K_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("K"), as.prob=T))*100
    pep_df$L_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("L"), as.prob=T))*100
    pep_df$M_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("M"), as.prob=T))*100
    pep_df$N_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("N"), as.prob=T))*100
    pep_df$P_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("P"), as.prob=T))*100
    pep_df$Q_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("Q"), as.prob=T))*100
    pep_df$R_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("R"), as.prob=T))*100
    pep_df$S_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("S"), as.prob=T))*100
    pep_df$T_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("T"), as.prob=T))*100
    pep_df$V_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("V"), as.prob=T))*100
    pep_df$W_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("W"), as.prob=T))*100
    pep_df$Y_prop[pep_df$Gene_id==i]=as.vector(letterFrequency(pep_cleaned[names(pep_cleaned) %in% i],letters=c("Y"), as.prob=T))*100

    setTxtProgressBar(pb, match(i,gene_id$Gene_id))
  }
  
  
  # pep_df <- left_join(pep_df,all_aacomp,by="Gene_id")
  if(sp=="Aeseve"){
    Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/projects/nodMimosa/ncr/Aeseve_NCRs.txt"),h=T)
    Medtru_annot$Gene_id <- paste0("Aeseve_",Medtru_annot$Gene_id)
    Medtru_annot$NCR <- "NCR"
    cds <- pep_df
    cds <- dplyr::left_join(cds,SignalP)
    cds <- dplyr::left_join(cds,Medtru_annot)
    cds$SignalP[is.na(cds$SignalP)] <- ""
    cds$NCR[is.na(cds$NCR)] <- "noNCR"
    cds <- dplyr::select(cds,-endSP)
  }else{
    Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
    Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)
    cds <- pep_df
    cds <- dplyr::left_join(cds,SignalP)
    cds <- dplyr::left_join(cds,Medtru_annot)
    cds$SignalP[is.na(cds$SignalP)] <- ""
    cds$Acronym[is.na(cds$Acronym)] <- ""
    # cds <- dplyr::left_join(cds,cds_NCR)
    cds <- dplyr::select(cds,-endSP)}
  
  fwrite(cds,paste0("./results_NCR/",sp,"_proteins_amino_acid_info.txt"),quote=FALSE, sep="\t")
  message(sp," Done!")
}


for(sp in Species){
results <-fread(paste0("./results_NCR/",sp,"_proteins_infos.txt"))
cds_NCR <- results %>% dplyr::filter(str_detect(Acronym, "^NCR")) %>% dplyr::select(Gene_id)
cds_NCR$NCR <- "NCR"
results <- dplyr::left_join(results,cds_NCR)
results$NCR[is.na(results$NCR)] <- "noNCR"
cds_200bp <- results %>% dplyr:::filter(Pep_length<200)

scatterPlot <- ggplot(cds_200bp,aes(Cysteine_prop, Pep_length, color=NCR)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 7)+
  geom_hline(yintercept = 65)

scatterPlot

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_200bp, aes(Cysteine_prop, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_200bp, aes(y=Pep_length, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/",sp,"_NCRs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()


}

######## Medicago NCR ########
cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
Medtru_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Medtru_cor_Genes_HOG <- Medtru_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Medtru_"))
cds <- cds[str_detect(names(cds), "^Medtru_")]
pep <- Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve")
# cds <- data.frame(Gene_id=names(cds),Length=width(cds)/3,GC=as.vector(letterFrequency(cds,letters=c("GC"), as.prob=T)),
#                   Cysteine=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
#                                                      letters=c("C"), as.prob=F)),
#                   Cysteine_prop=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
#                                                      letters=c("C"), as.prob=T))*100)


SignalP <- fread("../ncr/output.gff3",h=F)
SignalP <- SignalP %>% filter(str_detect(V1, "^Medtru_")) %>% dplyr::select(V1,V3,V5)
names(SignalP) <- c("Gene_id","SignalP","endSP")

pep <- pep[names(pep) %in% SignalP$Gene_id]
pep_cleaned <- pep
pep_df <- data.frame(Gene_id=names(pep))
pb <- txtProgressBar(min = 0, max = length(SignalP$Gene_id), style = 3)

# remove signalP from sequence
for(i in SignalP$Gene_id){
  subset <- SignalP %>% filter(Gene_id==i)
  pep_cleaned[names(pep_cleaned) %in% i]<-subseq(pep_cleaned[names(pep_cleaned) %in% i], start = subset$endSP+1)
  pep_df$mw[pep_df$Gene_id==i]=mw(as.character(pep_cleaned[names(pep_cleaned) %in% i]), monoisotopic = FALSE, avgScale = "expasy",label = "none", aaShift = NULL)
  pep_df$isoelectric_point[pep_df$Gene_id==i]=pI(pep_cleaned[names(pep_cleaned) %in% i], pKscale = "Bjellqvist")
  pep_df$hydrophobicity[pep_df$Gene_id==i]=hydrophobicity(as.character(pep_cleaned[names(pep_cleaned) %in% i]),scale = "KyteDoolittle")
  pep_df$Pep_length[pep_df$Gene_id==i]=width(pep_cleaned[names(pep_cleaned) %in% i])
  aacomp <- data.frame(Gene_id=t(as.data.frame(aaComp(as.character(pep_cleaned[names(pep_cleaned) %in% i]))))[2,])
  aacomp <- as.data.frame(t(aacomp)); aacomp$Gene_id <- i
  pep_df <- left_join(pep_df,aacomp,by="Gene_id")
  setTxtProgressBar(pb, match(i,SignalP$Gene_id))
}

Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)

cds <- pep_df
cds <- dplyr::left_join(cds,SignalP)
cds <- dplyr::left_join(cds,Medtru_annot)
cds$SignalP[is.na(cds$SignalP)] <- ""
cds$Acronym[is.na(cds$Acronym)] <- ""
cds_NCR <- cds %>% dplyr::filter(str_detect(Acronym, "^NCR")) %>% dplyr::select(Gene_id)
cds_NCR$NCR <- "NCR"
cds <- dplyr::left_join(cds,cds_NCR)
cds$NCR[is.na(cds$NCR)] <- "noNCR"
cds_putativeNCR <- cds %>% dplyr::filter(SignalP=="signal_peptide",Length<100,Cysteine_prop<=12,Cysteine_prop>=6) %>% dplyr::select(Gene_id)
cds_putativeNCR$putativeNCR <- "NCR"
cds <- dplyr::left_join(cds,cds_putativeNCR)
cds$putativeNCR[is.na(cds$putativeNCR)] <- "noNCR"

fwrite(cds,"./results_NCR/Medtru_NCRs_and_putativeNCRs.txt",sep="\t")

cds_200bp <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_200bp <- cds_200bp %>% dplyr::filter(Length<200)

p <- ggplot(cds_200bp, aes(Cysteine, Cysteine_prop, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

p <- ggplot(cds_200bp, aes(Cysteine_prop, Length, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

scatterPlot <- ggplot(cds_200bp,aes(Cysteine_prop, Length, color=NCR)) + 
  geom_point() + expand_limits(x=0)+
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_vline(xintercept = 12)+
  geom_hline(yintercept = 100)

scatterPlot 

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_200bp, aes(Cysteine_prop, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_200bp, aes(y=Length, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/Medtru_NCRs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
                  ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

######## Medicago NodGRP ########
cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
Medtru_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Medtru_cor_Genes_HOG <- Medtru_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Medtru_"))
cds <- cds[str_detect(names(cds), "^Medtru_")]
cds <- data.frame(Gene_id=names(cds),Length=width(cds)/3,GC=as.vector(letterFrequency(cds,letters=c("GC"), as.prob=T)),
                  Glycine=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                     letters=c("G"), as.prob=F)),
                  Glycine_prop=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                          letters=c("G"), as.prob=T))*100)

SignalP <- fread("../ncr/output.gff3",h=F)
SignalP <- SignalP %>% filter(str_detect(V1, "^Medtru_")) %>% dplyr::select(V1,V3)
names(SignalP) <- c("Gene_id","SignalP")

Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/results/annotation/Medtru_annotation_and_acronym.txt"),h=T)
Medtru_annot <- Medtru_annot %>% dplyr::select(Gene_id,Acronym)

cds <- dplyr::left_join(cds,SignalP)
cds <- dplyr::left_join(cds,Medtru_annot)
cds$SignalP[is.na(cds$SignalP)] <- ""
cds$Acronym[is.na(cds$Acronym)] <- ""
cds_NodGRP <- cds %>% dplyr::filter(str_detect(Acronym, "^NodGRP")) %>% dplyr::select(Gene_id)
cds_NodGRP$NodGRP <- "NodGRP"
cds <- dplyr::left_join(cds,cds_NodGRP)
cds$NodGRP[is.na(cds$NodGRP)] <- "noNodGRP"

cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Length<200)

p <- ggplot(cds_1kb, aes(Glycine, Length, color = NodGRP=="NodGRP"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

p <- ggplot(cds_1kb, aes(Glycine_prop, Length, color = NodGRP=="NodGRP"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

scatterPlot <- ggplot(cds_1kb,aes(Glycine_prop, Length, color=NodGRP)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() #+ stat_ellipse()
scatterPlot 

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_1kb, aes(Glycine_prop, fill=NodGRP)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_1kb, aes(y=Length, fill=NodGRP)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/Medtru_putative_NodGRPs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 4), heights=c(1.4, 4))
dev.off()


######## Aeseve NCRs ######## 
cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
Medtru_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Medtru_cor_Genes_HOG <- Medtru_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Aeseve_"))
cds <- cds[str_detect(names(cds), "^Aeseve_")]
cds <- data.frame(Gene_id=names(cds),Length=width(cds)/3,GC=as.vector(letterFrequency(cds,letters=c("GC"), as.prob=T)),
                  Cysteine=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                     letters=c("C"), as.prob=F)),
                  Cysteine_prop=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                          letters=c("C"), as.prob=T))*100)

SignalP <- fread("../ncr/output.gff3",h=F)
SignalP <- SignalP %>% filter(str_detect(V1, "^Aeseve_")) %>% dplyr::select(V1,V3)
SignalP$V1 <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",SignalP$V1)
names(SignalP) <- c("Gene_id","SignalP")

Medtru_annot <- fread(paste0("//194.199.55.66/evo/commun/projects/nodMimosa/ncr/Aeseve_NCRs.txt"),h=T)
Medtru_annot$Gene_id <- paste0("Aeseve_",Medtru_annot$Gene_id)
Medtru_annot$NCR <- "NCR"

cds <- dplyr::left_join(cds,SignalP)
cds <- dplyr::left_join(cds,Medtru_annot)
cds$SignalP[is.na(cds$SignalP)] <- ""
cds$NCR[is.na(cds$NCR)] <- "noNCR"

cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Length<200)

p <- ggplot(cds_1kb, aes(Cysteine, Length, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

p <- ggplot(cds_1kb, aes(Cysteine_prop, Length, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

scatterPlot <- ggplot(cds_1kb,aes(Cysteine_prop, Length, color=NCR)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_vline(xintercept = 12)+
  geom_hline(yintercept = 100)
scatterPlot 

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_1kb, aes(Cysteine_prop, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_1kb, aes(y=Length, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/Aeseve_NCRs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()



######## Mimpud putative NCRs ########
cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
Medtru_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
Medtru_cor_Genes_HOG <- Medtru_cor_Genes_HOG %>% filter(str_detect(Gene_id, "^Mimpud_"))
cds <- cds[str_detect(names(cds), "^Mimpud_")]
cds <- data.frame(Gene_id=names(cds),Length=width(cds)/3,GC=as.vector(letterFrequency(cds,letters=c("GC"), as.prob=T)),
                  Cysteine=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                     letters=c("C"), as.prob=F)),
                  Cysteine_prop=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
                                                          letters=c("C"), as.prob=T))*100)

SignalP <- fread("../ncr/output.gff3",h=F)
SignalP <- SignalP %>% filter(str_detect(V1, "^Mimpud_")) %>% dplyr::select(V1,V3)
SignalP$V1 <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",SignalP$V1)
names(SignalP) <- c("Gene_id","SignalP")

cds <- dplyr::left_join(cds,SignalP)
cds$SignalP[is.na(cds$SignalP)] <- ""

cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Length<200)
cds_putativeNCR <- cds_1kb %>% dplyr::filter(SignalP=="signal_peptide",Length<85,Cysteine_prop<=12,Cysteine_prop>=6) %>% dplyr::select(Gene_id)
cds_putativeNCR$NCR <- "NCR"
cds <- dplyr::left_join(cds,cds_putativeNCR)
cds$NCR[is.na(cds$NCR)] <- "noNCR"
cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Length<200)

p <- ggplot(cds_1kb, aes(Cysteine, Length, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

p <- ggplot(cds_1kb, aes(Cysteine_prop, Length, color = NCR=="NCR"))+
  geom_point()
p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
p + stat_ellipse(type = "norm")

scatterPlot <- ggplot(cds_1kb,aes(Cysteine_prop, Length, color=NCR)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = 6)+
  geom_vline(xintercept = 12)+
  geom_hline(yintercept = 85)
scatterPlot 

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_1kb, aes(Cysteine_prop, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity
# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_1kb, aes(y=Length, fill=NCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity
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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/Mimpud_putative_NCRs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

######## All species putative NCRs ########
species=c("Aeseve","Arahyp","Datglo","Medtru","Mimpud","Parand","Lupalb","Hiprha")
all_NCRs <- c()

Length_th=65
Cysteine_prop_th_up=12
Cysteine_prop_th_down=7


for(sp in species){
cds <-fread(paste0("./results_NCR/",sp,"_proteins_infos.txt"))
# cds <- readDNAStringSet(paste0("../Results_Apr14/nodMimpud_cds_db.fas"))
# names(cds) <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",names(cds))
# Medtru_cor_Genes_HOG <- fread("./results/N0_corrected_Gene_id_HOGs_correspondance.txt",h=T)
# Medtru_cor_Genes_HOG <- Medtru_cor_Genes_HOG %>% filter(str_detect(Gene_id, paste0("^",sp,"_")))
# cds <- cds[str_detect(names(cds), paste0("^",sp,"_"))]
# cds <- data.frame(Gene_id=names(cds),Length=width(cds)/3,GC=as.vector(letterFrequency(cds,letters=c("GC"), as.prob=T)),
#                   Cysteine=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
#                                                      letters=c("C"), as.prob=F)),
#                   Cysteine_prop=as.vector(letterFrequency(Biostrings::translate(cds,no.init.codon=TRUE,if.fuzzy.codon="solve"),
#                                                           letters=c("C"), as.prob=T))*100)
# 
# SignalP <- fread("../ncr/output.gff3",h=F)
# SignalP <- SignalP %>% filter(str_detect(V1, paste0("^",sp,"_"))) %>% dplyr::select(V1,V3)
# SignalP$V1 <- gsub("\\.[0-9]+$|\\.t[0-9]+$", "",SignalP$V1)
# names(SignalP) <- c("Gene_id","SignalP")
# 
# cds <- dplyr::left_join(cds,SignalP)
# cds$SignalP[is.na(cds$SignalP)] <- ""

cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Pep_length<200)
cds_putativeNCR <- cds_1kb %>% dplyr::filter(SignalP=="signal_peptide",Pep_length<Length_th,Cysteine_prop>=Cysteine_prop_th_down) %>% dplyr::select(Gene_id)
cds_putativeNCR$putativeNCR <- "NCR"
all_NCRs <- rbind(all_NCRs,cds_putativeNCR)
cds <- dplyr::left_join(cds,cds_putativeNCR, by="Gene_id")
cds$putativeNCR[is.na(cds$putativeNCR)] <- "noNCR"
cds_1kb <- cds %>% dplyr::filter(SignalP=="signal_peptide")
cds_1kb <- cds_1kb %>% dplyr::filter(Pep_length<200)

# p <- ggplot(cds_1kb, aes(Cysteine, Pep_length, color = putativeNCR=="NCR"))+
#   geom_point()
# p + stat_ellipse()
# Changer le type d'ellipses: 
# Valeurs possibles "t", "norm", "euclid"
# p + stat_ellipse(type = "norm")

# p <- ggplot(cds_1kb, aes(Cysteine_prop, Length, color = NCR=="NCR"))+
#   geom_point()
# p + stat_ellipse()
# # Changer le type d'ellipses: 
# # Valeurs possibles "t", "norm", "euclid"
# p + stat_ellipse(type = "norm")

scatterPlot <- ggplot(cds_1kb,aes(Cysteine_prop, Pep_length, color=putativeNCR)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.9,1), legend.justification=c(1,1)) + geom_density_2d() + #+ stat_ellipse()
  geom_vline(xintercept = Cysteine_prop_th_down)+
  # geom_vline(xintercept = 12)+
  geom_hline(yintercept = Length_th)
scatterPlot 

# Courbe de densité marginale de x (panel du haut)
xdensity <- ggplot(cds_1kb, aes(Cysteine_prop, fill=putativeNCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
xdensity

# Courbe de densité marginale de y (panel de droite)
ydensity <- ggplot(cds_1kb, aes(y=Pep_length, fill=putativeNCR)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity

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
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

pdffile = paste0("./results_NCR/",sp,"_putative_NCRs.pdf")
#adjust the size and margin of the pdf file
pdf(pdffile,height=15,width=15)
# par(mar=c(5,5,5,5))
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, 
             ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
dev.off()

fwrite(cds,paste0("./results_NCR/",sp,"_putative_NCRs.txt"), sep="\t", quote=FALSE)

}



