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
                       files = c("../../../10_metaGWAS/01_eGFR_allEth_sex_combined/08_credSets/gwasresults_V5/synopsis/topliste_tabdelim/topliste_2022-11-04_credSets.txt",
                                 "../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_female_V5/synopsis/topliste_tabdelim/topliste_2022-11-04_credSets_female.txt",
                                 ("../../../10_metaGWAS/01_eGFR_allEth_sex_stratified/08_credSets/gwasresults_male_V4/synopsis/topliste_tabdelim/topliste_2022-07-22_credSets_male.txt"),
                                 ("../../../10_metaGWAS/03_uric_acid_allEth_sex_combined/08_credSets/gwasresults_V4/synopsis/topliste_tabdelim/topliste_2022-07-21_credSets.txt"),
                                 ("../../../10_metaGWAS/03_uric_acid_allEth_sex_stratified/08_credSets/gwasresults_male_V4/synopsis/topliste_tabdelim/topliste_2022-07-22_credSets_male.txt")))

tab7 = foreach(i=1:dim(ToDoList3)[1])%do%{
  #i=1
  myRow = ToDoList3[i,]
  
  tab = fread(myRow$files)
  
  stopifnot(tab$litsnp == tab$provided_name)
  
  myNames2 =  c("snp","CredSet","PostProb","SumProb","cyto","pos","tagger","r2_tagger","tagsnp","effect_allele","other_allele","eaf",
                "beta","SE","logP", "I2","nearestgenes","gene_biotype","nearestgene","Eigen","EigenPC","CADD_scaled","DANN",
                "GWAVA_Region","GWAVA_TSS","GWAVA_Unmatched","regulome_score","regulome_score_numeric","regulome_details",
                "corinfo_gwas2","cisgene","transgene","coremine_genes","hgnc4pathway","entrez4pathway","KEGG","reactome","DOSE",
                "GO","tissues")
  myNames = names(tab)[c(1:4,13,15:18,61,62,64,70:72,74,19:42)]
  stopifnot(sum(myNames == myNames2)==34)
  colsOut<-setdiff(colnames(tab),myNames)
  tab[,get("colsOut"):=NULL]
  setcolorder(tab,myNames)
  names(tab) = myNames2
  tab[,phenotype := myRow$pheno]
  tab
}
tab7 = rbindlist(tab7)
table(tab7$phenotype)

#' # Get plot data ####
#' ***
#' I need Posterior Probability, Region, Phenotype, N per Region and Phenotype, CADD score
plotData = copy(tab7)
plotData = plotData[,c(1:4,41,22)]
plotData

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

tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize.tiff", 
     width = 2500 , height = 2000, res=300, compression = 'lzw')
myPlot3
dev.off()

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

tiff(filename = "../figures/MainFigure4_PostProbByCredSetSize_N400.tiff", 
     width = 2500 , height = 2000, res=300, compression = 'lzw')
myPlot4
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

