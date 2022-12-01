#' ---
#' title: "Credible Set Size vs. Posterior Probability"
#' subtitle: "CKDGen - Chr X"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(paste0(projectpath,"scripts/"))

#' # Load data ####
#' ***
ToDoList3 = data.table(pheno = c("eGFR_ALL","eGFR_FEMALE","eGFR_MALE","UA_ALL","UA_MALE"),
                       files = c(paste0(path_CS_eGFR_ALL,"topliste_2022-11-29_credSets.txt"),
                                 paste0(path_CS_eGFR_FEMALE,"topliste_2022-11-29_credSets_female.txt"),
                                 paste0(path_CS_eGFR_MALE,"topliste_2022-11-29_credSets_male.txt"),
                                 paste0(path_CS_UA_ALL,"topliste_2022-11-29_credSets.txt"),
                                 paste0(path_CS_UA_MALE,"topliste_2022-11-29_credSets_male.txt")))

tab7 = foreach(i=1:dim(ToDoList3)[1])%do%{
  #i=3
  myRow = ToDoList3[i,]
  
  tab = fread(myRow$files)
  
  stopifnot(tab$litsnp == tab$provided_name)
  
  myNames2 =  c("snp","CredSet","PostProb","SumProb","cyto","pos","tagger","r2_tagger","tagsnp",
                "effect_allele","other_allele","eaf",
                "beta","SE","logP", "I2","nearestgenes","gene_biotype","nearestgene","Eigen","EigenPC","CADD_scaled","DANN",
                "GWAVA_Region","GWAVA_TSS","GWAVA_Unmatched","regulome_score","regulome_score_numeric","regulome_details",
                "corinfo_gwas2","cisgene","transgene","coremine_genes","hgnc4pathway","entrez4pathway",
                "KEGG","reactome","DOSE","GO","tissues")
  myNames = names(tab)[c(1:4,13,15:18,61,62,64,70:72,74,19:42)]
  stopifnot(sum(myNames == myNames2)==34)
  colsOut<-setdiff(colnames(tab),myNames)
  tab[,get("colsOut"):=NULL]
  setcolorder(tab,myNames)
  names(tab) = myNames2
  tab[,phenotype := myRow$pheno]
  if(myRow$pheno == "eGFR_MALE"){
    tab[cyto == "Xq22.1" & CredSet == "Region1", CredSet := "Region7"]
  }
  if(myRow$pheno == "UA_ALL"){
    tab[CredSet == "Region21_SNP1", CredSet := "Region21_rs202138804"]
    tab[CredSet == "Region21_SNP2", CredSet := "Region21_rs7056552"]
    tab[CredSet == "Region22_SNP1", CredSet := "Region22_rs111884516"]
    tab[CredSet == "Region22_SNP2", CredSet := "Region22_rs4328011"]
  }
  tab
  
}

tab7 = rbindlist(tab7)
table(tab7$phenotype)

#' # Get plot data ####
#' ***
#' I need Posterior Probability, Region, Phenotype, N per Region and Phenotype, CADD score, Missense Mutation.
#' 
#' In addition, I only want one cred set per region (best phenotype per region)
#' 
locusDef = fread("../results/01_Locus_Definitions.txt")
locusDef[,dumID := paste0("Region",region,"_",phenotype)]
locusDef[,dumID]

plotData = copy(tab7)
plotData[,dumID := paste0(CredSet,"_",phenotype)]
unique(plotData$dumID)
plotData[,dumID := gsub("Region9_rs.*","Region9_eGFR_ALL",dumID)]
plotData[,dumID := gsub("Region21_rs.*","Region21_UA_ALL",dumID)]
plotData[,dumID := gsub("Region22_rs.*","Region22_UA_ALL",dumID)]
plotData = plotData[dumID %in% c(locusDef$dumID,"Region7_eGFR_FEMALE"),]

length(unique(plotData$dumID))
length(unique(plotData$CredSet))
plotData[,dumID2 := paste0(CredSet,"_",phenotype)]
length(unique(plotData$dumID2))

names(plotData)
plotData[,missense := F]
plotData[grepl("Missense",nearestgenes),missense := T]

plotData = plotData[,c(1:4,41,22,44)]
plotData
plotData[,rsID := gsub(":.*","",snp)]
plotData[,table(rsID=="chr23",missense)]

plotData[,dumID := paste(phenotype,CredSet,sep="_")]
length(unique(plotData$dumID))
dummy = plotData[,.N, by=dumID]
dummy

matched = match(plotData$dumID,dummy$dumID)
plotData[,N := dummy[matched,N]]

plotData[,table(CADD_scaled>10)]
plotData[,table(CADD_scaled>20)]
plotData[,CADD_type := 0]
plotData[CADD_scaled>10,  CADD_type :=CADD_type+1]
plotData[CADD_scaled>20,  CADD_type :=CADD_type+1]

setorder(plotData,CADD_type)
plotData[,CADD_type := as.factor(CADD_type)]
plotData[,table(CADD_type,missense)]

#' # Plotting ####
#' ***
#' 
#' ## Plot 1 ####
#' facet plot per phenotypes (no limitation on x-axis)

myPlot1 = ggplot(plotData, aes(x=N, y=PostProb,color=CADD_type)) + 
  facet_wrap(~phenotype, scales = "free") +
  geom_point(size = 2) +
  theme_bw(base_size = 10)+
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC"),
                      labels=c("CAD<=10","CAD in (10,20]","CAD>20"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score")
myPlot1

#' ## Plot 2 ####
#' facet plot per phenotypes (no limitation on x-axis, no legend)

myPlot2 = ggplot(plotData, aes(x=N, y=PostProb,color=CADD_type)) + 
  facet_wrap(~phenotype, scales = "free") +
  geom_point(size = 2) +
  theme_bw(base_size = 10)+
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC"),
                      labels=c("CAD<=10","CAD in (10,20]","CAD>20"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability") +
  guides(label="none",color="none")
myPlot2

#' ## Plot 3 ####
#' one plot over all phenotypes (no limitation on x-axis)

myPlot3 = ggplot(plotData[CADD_scaled>=10,], aes(x=N, y=PostProb,color=CADD_type)) + 
  # facet_wrap(~phenotype, scales = "free") +
  geom_point(data=plotData[CADD_scaled<10,], aes(x=N, y=PostProb),col="black",size=2.5,alpha=0.5,shape = 16)+
  geom_point(size=3,alpha=0.75) + 
  theme_bw(base_size = 10)+
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC"),
                      labels=c("CAD<=10","CAD in (10,20]","CAD>20"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score")
myPlot3

# tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize.tiff", 
#      width = 2500 , height = 2000, res=300, compression = 'lzw')
# myPlot3
# dev.off()

#' ## Plot 4 ####
#' one plot over all phenotypes (x-axis limited to 400)

myPlot4 = ggplot(plotData[CADD_scaled>=10 & N<400,], aes(x=N, y=PostProb,color=CADD_type)) + 
  # facet_wrap(~phenotype, scales = "free") +
  geom_point(data=plotData[CADD_scaled<10 & N<400,], aes(x=N, y=PostProb),col="black",size=2.5,alpha=0.5,shape = 16)+
  geom_point(size=3,alpha=0.75) + 
  theme_bw(base_size = 10)+
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC"),
                      labels=c("CAD<=10","CAD in (10,20]","CAD>20"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score")
myPlot4

# tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize_N400.tiff", 
#      width = 2500 , height = 2000, res=300, compression = 'lzw')
# myPlot4
# dev.off()

#' ## Plot 5 ####
#' one plot over all phenotypes (x-axis limited to 400, missense mutations labeled)

myPlot5 = ggplot(plotData[CADD_scaled>=10 & N<400,], aes(x=N, y=PostProb,color=CADD_type)) + 
  # facet_wrap(~phenotype, scales = "free") +
  geom_point(data=plotData[CADD_scaled<10 & N<400,], aes(x=N, y=PostProb),col="black",size=2.5,alpha=0.5,shape = 16)+
  geom_point(size=3,alpha=0.75) + 
  theme_bw(base_size = 10)+
  scale_colour_manual(values=c("#000000","#B2182B","#2166AC"),
                      labels=c("CAD<=10","CAD in (10,20]","CAD>20"))+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  geom_label_repel(data = subset(plotData, missense==T & CADD_scaled> 10 & N<400),
                   aes(x=N, y=PostProb, label = rsID),
                   size=4, direction = 'y',
                   xlim = c(150,Inf), ylim = c(0.5,Inf),box.padding = 0.5, max.overlaps = Inf)+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score")+
  guides(label="none")

myPlot5

#' ## Plot 6 ####
#' one plot over all phenotypes (x-axis limited to 400, missense mutations labeled, shape by phenotype) 
#' 
#' THX Carl for plotting support!

plotData2 = copy(plotData)
plotData2[missense==T & CADD_scaled>10,]
SNPInfo = c("SLC25A5 \np.(Leu111Arg)","SLC25A43 \np.(Pro334Leu)","ARHGAP4 \np.(Glu602=)", "TSPAN6 \np.(Ala108Thr)",
            "MID2 \np.(Ala378Asp)", "SERPINA7 \np.(Leu303Phe)","PLXNB3 \np.(Val1596Glu)")
plotData2[missense==T & CADD_scaled>10,rsID := paste(rsID, "\n",SNPInfo)]
plotData2 = plotData2[N<400,]

myPlot6 = ggplot(
  plotData2, 
  aes(x = N, 
      y = PostProb,
      color = CADD_type,
      size = CADD_type,
      alpha = CADD_type,
      shape = phenotype
  )
) + 
  # facet_wrap(~phenotype, scales = "free") +
  geom_point() +
  # scale_fill_manual() +
  scale_shape_manual(values = c(16, 17, 15, 3),
                     labels = c("eGFR ALL","eGFR FEMALE","eGFR MALE","UA ALL")) +
  scale_size_manual(values = c(2.5, 3, 3)) +
  scale_alpha_manual(values = c(0.5, 0.75, 0.75)) +
  
  # geom_point(data=plotData[CADD_scaled<10 & N<400,], aes(x=N, y=PostProb),col="black",size=2.5,alpha=0.5,shape = 16)+
  # geom_point(size=3,alpha=0.75) + 
  
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  geom_label_repel(data = subset(plotData2, missense==T & CADD_scaled> 10 & N<400),
                   aes(x=N, y=PostProb, label = rsID),
                   size=4, direction = 'y',
                   xlim = c(150,Inf), ylim = c(0.5,1.05),box.padding = 0.5, max.overlaps = Inf, 
                   show.legend = FALSE
  )+
  scale_colour_manual(
    values = c("#000000","#B2182B","#2166AC"),
    labels = c("CAD<=10","CAD in (10,20]","CAD>20"))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score")+
  guides(
    label="none",size = "none", alpha = "none", fill = "none")

myPlot6

tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize_N400_MissenseMut.tiff", 
     width = 2500 , height = 2000, res=300, compression = 'lzw')
myPlot6
dev.off()

#' ## Plot 7 ####
#' one plot over all phenotypes (x-axis not limited, missense mutations labeled, shape by phenotype) 

plotData[missense==T & CADD_scaled>10,rsID := paste(rsID, "\n",SNPInfo)]

plotData[missense==T & CADD_scaled>10,]

myPlot7 = ggplot(
  plotData, 
  aes(x = N, 
      y = PostProb,
      color = CADD_type,
      size = CADD_type,
      alpha = CADD_type,
      shape = phenotype
  )
) + 
  # facet_wrap(~phenotype, scales = "free") +
  geom_point() +
  # scale_fill_manual() +
  scale_shape_manual(values = c(16, 17, 15, 3),
                     labels = c("eGFR ALL","eGFR FEMALE","eGFR MALE","UA ALL")) +
  scale_size_manual(values = c(2.5, 3, 3)) +
  scale_alpha_manual(values = c(0.5, 0.75, 0.75)) +
  
  # geom_point(data=plotData[CADD_scaled<10 & N<400,], aes(x=N, y=PostProb),col="black",size=2.5,alpha=0.5,shape = 16)+
  # geom_point(size=3,alpha=0.75) + 
  
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text = element_text(size = 20))+
  geom_label_repel(data = subset(plotData, missense==T & CADD_scaled> 10 & N<400),
                   aes(x=N, y=PostProb, label = rsID),
                   size=4, direction = 'y',
                   xlim = c(150,Inf), ylim = c(0.5,1.05),box.padding = 0.5, max.overlaps = Inf, 
                   show.legend = FALSE
  )+
  geom_label_repel(data = subset(plotData, missense==T & CADD_scaled> 10 & N>1500 ),
                   aes(x=N, y=PostProb, label = rsID),
                   size=3, direction = 'y',
                   xlim = c(1500,1800), ylim = c(0,0.5),
                   box.padding = 0.5, max.overlaps = Inf, 
                   show.legend = FALSE
  )+
  geom_label_repel(data = subset(plotData, missense==T & CADD_scaled> 10 & N>400 & N<1500 ),
                   aes(x=N, y=PostProb, label = rsID),
                   size=3, direction = 'y',
                   xlim = c(800,1300), ylim = c(0,0.5),
                   box.padding = 0.5, max.overlaps = Inf, 
                   show.legend = FALSE
  )+
  scale_colour_manual(
    values = c("#000000","#B2182B","#2166AC"),
    labels = c("CAD<=10","CAD in (10,20]","CAD>20"))+
  labs(x="Credible Set Size", 
       y = "Posterior Probability",
       color = "CADD Score", shape = "Phenotype and \nsetting")+
  guides(
    label="none",size = "none", alpha = "none", fill = "none")

myPlot7

tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize_MissenseMut.tiff", 
     width = 2500 , height = 2000, res=300, compression = 'lzw')
myPlot7
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

