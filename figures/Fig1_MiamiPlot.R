#' ---
#' title: "Miami Plot"
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
source("../helperFunctions/miamiPlot.R")

setwd(paste0(projectpath,"figures/"))

#' # Load data ####
#' ***
eGFR_male   = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
eGFR_female = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
eGFR_all    = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
UA_male     = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
UA_female   = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
UA_all      = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")

#' # Check sex-specificity ####
#' ***
erg1 = rbind(eGFR_male,eGFR_female,eGFR_all,UA_male,UA_female,UA_all)
mySNPs_eGFR = unique(erg1[grepl("eGFR",phenotype) & P<=5e-8 & invalid_assoc ==F,rsID])
mySNPs_UA = unique(erg1[grepl("UA",phenotype) & P<=5e-8 & invalid_assoc ==F,rsID])

eGFR_male2 = copy(eGFR_male)[rsID %in% mySNPs_eGFR,]
eGFR_female2 = copy(eGFR_female)[rsID %in% mySNPs_eGFR,]
UA_male2 = copy(UA_male)[rsID %in% mySNPs_UA,]
UA_female2 = copy(UA_female)[rsID %in% mySNPs_UA,]
table(eGFR_female2$rsID == eGFR_male2$rsID)
table(UA_female2$rsID == UA_male2$rsID)

dumTab1 = foreach(i = 1:dim(eGFR_female2)[1])%do%{
  #i=1
  test = interactionTest(mean1 = eGFR_male2[i,beta],
                         se1 = eGFR_male2[i,SE],
                         mean2 = eGFR_female2[i,beta],
                         se2 = eGFR_female2[i,SE])
  
  res = data.table(rsID = eGFR_female2[i,rsID],
                   pos = eGFR_female2[i,position],
                   EA = eGFR_female2[i,effect_allele],
                   OA = eGFR_female2[i,other_allele],
                   NStudies_F = eGFR_female2[i,numberOfStudies],
                   NSamples_F = eGFR_female2[i,N],
                   info_F = eGFR_female2[i,infoScore],
                   EAF_F = eGFR_female2[i,EAF],
                   beta_F = eGFR_female2[i,beta],
                   SE_F = eGFR_female2[i,SE],
                   pval_F = eGFR_female2[i,P],
                   I2_F = eGFR_female2[i,I2],
                   NStudies_M = eGFR_male2[i,numberOfStudies],
                   NSamples_M = eGFR_male2[i,N],
                   info_M = eGFR_male2[i,infoScore],
                   EAF_M = eGFR_male2[i,EAF],
                   beta_M = eGFR_male2[i,beta],
                   SE_M = eGFR_male2[i,SE],
                   pval_M = eGFR_male2[i,P],
                   I2_M = eGFR_male2[i,I2], 
                   IA_diff = test$meandiff,
                   IA_SE = test$meandiff_se, 
                   IA_pval = test$meandiff_p)
  
  res
}
IA_eGFR = rbindlist(dumTab1)
dumTab2 = foreach(i = 1:dim(UA_female2)[1])%do%{
  #i=1
  test = interactionTest(mean1 = UA_male2[i,beta],
                         se1 = UA_male2[i,SE],
                         mean2 = UA_female2[i,beta],
                         se2 = UA_female2[i,SE])
  
  res = data.table(rsID = UA_female2[i,rsID],
                   pos = UA_female2[i,position],
                   EA = UA_female2[i,effect_allele],
                   OA = UA_female2[i,other_allele],
                   NStudies_F = UA_female2[i,numberOfStudies],
                   NSamples_F = UA_female2[i,N],
                   info_F = UA_female2[i,infoScore],
                   EAF_F = UA_female2[i,EAF],
                   beta_F = UA_female2[i,beta],
                   SE_F = UA_female2[i,SE],
                   pval_F = UA_female2[i,P],
                   I2_F = UA_female2[i,I2],
                   NStudies_M = UA_male2[i,numberOfStudies],
                   NSamples_M = UA_male2[i,N],
                   info_M = UA_male2[i,infoScore],
                   EAF_M = UA_male2[i,EAF],
                   beta_M = UA_male2[i,beta],
                   SE_M = UA_male2[i,SE],
                   pval_M = UA_male2[i,P],
                   I2_M = UA_male2[i,I2], 
                   IA_diff = test$meandiff,
                   IA_SE = test$meandiff_se, 
                   IA_pval = test$meandiff_p)
  
  res
}
IA_UA = rbindlist(dumTab2)

IA_eGFR[,phenotype := "eGFR"]
IA_UA[,phenotype := "UA"]
IA_eGFR[,IA_pval_adj := p.adjust(p=IA_pval, method = "fdr")]
IA_UA[,IA_pval_adj := p.adjust(p=IA_pval, method = "fdr")]
IA_eGFR[,IA_pval_adj2 := p.adjust(p=IA_pval, method = "bonferroni")]
IA_UA[,IA_pval_adj2 := p.adjust(p=IA_pval, method = "bonferroni")]

IA = rbind(IA_eGFR,IA_UA)
myFDR<- addHierarchFDR(pvalues = IA[,IA_pval], 
                       categs = IA[,phenotype],quiet = F)
IA[,hierFDR := myFDR$hierarch_fdr5proz]
IA[,table(hierFDR,phenotype)]
IA[,table(hierFDR,IA_pval_adj<=0.05)]
IA[,table(hierFDR,IA_pval_adj2<=0.05)]

IA[hierFDR==T & pval_F <5e-8,sexType := "female-specific"]
IA[hierFDR==T & pval_M <5e-8,sexType := "male-specific"]
IA[hierFDR==F ,sexType := NA]

#' # Filter data ####
#' ***
#' Exclude all invalid associations (info, maf, I2, number of studies) and keep only SNPs with p<=0.05 for plotting
eGFR_male   = eGFR_male[invalid_assoc==F,]
eGFR_female = eGFR_female[invalid_assoc==F,]
eGFR_all    = eGFR_all[invalid_assoc==F,]
UA_male     = UA_male[invalid_assoc==F,]
UA_female   = UA_female[invalid_assoc==F,]
UA_all      = UA_all[invalid_assoc==F,]

eGFR_male   = eGFR_male[P<=0.05,]
eGFR_female = eGFR_female[P<=0.05,]
eGFR_all    = eGFR_all[P<=0.05,]
UA_male     = UA_male[P<=0.05,]
UA_female   = UA_female[P<=0.05,]
UA_all      = UA_all[P<=0.05,]

#' # Merge data ####
#' ***
#' Step 1: generate eGFR and UA object with best p-value per SNP over all three settings
eGFR<-rbind(eGFR_male,eGFR_female,eGFR_all)
setorder(eGFR,"P")
eGFR = eGFR[!duplicated(rsID),]
eGFR[,table(phenotype)]

UA<-rbind(UA_male,UA_female,UA_all)
setorder(UA,"P")
UA = UA[!duplicated(rsID),]
UA[,table(phenotype)]

#' Step 2: merge the two traits 
plotData = rbind(eGFR,UA)

myNames<-c("rsID","chromosome","position","P","logP","phenotype")
colsOut<-setdiff(colnames(plotData),myNames)
plotData[,get("colsOut"):=NULL]
head(plotData)
plotData[grepl("eGFR",phenotype),flag:="top"]
plotData[grepl("UA",phenotype),flag:="bottom"]
setnames(plotData,"rsID","SNP")
setnames(plotData,"chromosome","CHR")
setnames(plotData,"position","BP")
setorder(plotData,CHR,BP)
head(plotData)

#' # Add genes and novelty ####
#' ***
#' I add the information about candidate genes, novelty, and sex-specificity manually, as there is no other way right now. 
#'
#' * Candidate genes: taken from excel sheet table1_2022-10-15.xlsx (stored locally - not pretty!)
#' * Novelty: taken from excel sheet table1_2022-10-15.xlsx (stored locally - not pretty!)
#' * IA 
myTab = fread("../results/01_Locus_Definitions.txt")
myTab[,CandidateGenes := c("FAM9B","CDKL5","CDK16","EDA2R, AR","BRWD3","TSPAN6","ARMCX4, \nBEX4",
                           "BEX, \nTCEAL","CLDN2","ACSL4","SLC25A43","DCAF12L1","MST4","HPRT1",
                           "DUSP9","EDA2R, \nAR","CITED1, \nPIN4","ARMCX, \nADH4","BEX, \nTCEAL","DCAF12L1","HPRT1","DUSP9")]
myTab[,Novelty := c(F,F,F,T,T,T,F,F,T,T,F,F,F,F,F,T,F,F,F,F,F,F)]
myTab[,SexInteraction := c("male-specific",NA,NA,"male-specific",NA,NA,NA,
                           NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)]

stopifnot(sum(is.element(myTab$rsID,plotData$SNP))==22)
stopifnot(sum(is.element(myTab$rsID,IA$rsID))==22)

plotData[SNP %in% myTab$rsID]
plotData[,matchID := paste(SNP,phenotype,sep="::")]
myTab[,matchID := paste(rsID,phenotype,sep="::")]

matched1 = match(plotData$matchID,myTab$matchID)
plotData[,region_num := myTab[matched1,region]]
plotData[,region_chr := paste("region",myTab[matched1,region])]

plotData[,candidateGene :=myTab[matched1,CandidateGenes]]
plotData[,Novelty :=myTab[matched1,Novelty]]
plotData[,SexInteraction := myTab[matched1,SexInteraction]]

plotData[SNP == "rs149995096:100479327:C:T",region_num := 7]
plotData[SNP == "rs149995096:100479327:C:T",candidateGene := "DRP2"]
plotData[SNP == "rs149995096:100479327:C:T",Novelty := T]
plotData[SNP == "rs149995096:100479327:C:T",SexInteraction := "female-specific"]

#' # Plot ####
#' ***
mytitle = paste0("Miami Plot; top: eGFR, buttom: UA") 

ymaxpar1 = ifelse(plotData[flag=="top",max(logP,na.rm=T)] <7, 8,plotData[flag=="top",max(logP,na.rm=T)]+1)
ymaxpar2 = ifelse(plotData[flag=="bottom",max(logP,na.rm=T)] <7, 8,plotData[flag=="bottom",max(logP,na.rm=T)]+1)

plot1 = miamiPlot(x=plotData,
                 ymax = ymaxpar1,
                 ymin = -ymaxpar2,
                 title = mytitle,
                 xlabel = "",
                 ylabel=expression(paste("UA: ",log[10](p),"; eGFR: ",-log[10](p))),
                 hline1=-log10(5e-8),hline2=log10(5e-8),
                 sugline1=-log10(1e-6),sugline2=log10(1e-6),
                 highlight=T, diffsize = T,num_breaks_y=10,
                 plotGenes=T,
                 out_name="Fig1_MiamiPlot_BasePosition.pdf",
                 returnObject = T,
                 overall_max = 40, overall_min = -40,useBasePosition = T)

plot1

plot2 = miamiPlot(x=plotData,
                  ymax = ymaxpar1,
                  ymin = -ymaxpar2,
                  title = mytitle,
                  xlabel = "",
                  ylabel=expression(paste("UA: ",log[10](p),"; eGFR: ",-log[10](p))),
                  hline1=-log10(5e-8),hline2=log10(5e-8),
                  sugline1=-log10(1e-6),sugline2=log10(1e-6),
                  highlight=T, diffsize = T,num_breaks_y=10,
                  plotGenes=T,
                  out_name="Fig1_MiamiPlot_IndexPosition.pdf",
                  returnObject = T,
                  overall_max = 40, overall_min = -40,useBasePosition = F)

plot2

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


