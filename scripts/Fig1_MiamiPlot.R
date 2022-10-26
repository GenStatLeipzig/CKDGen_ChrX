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

setwd(projectpath_main)

#' # Load data ####
#' ***
eGFR_male   = fread("../data/CKDGen_ChrX_sumStat_eGFR_MALE.gz")
eGFR_female = fread("../data/CKDGen_ChrX_sumStat_eGFR_FEMALE.gz")
eGFR_all    = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
UA_male     = fread("../data/CKDGen_ChrX_sumStat_UA_MALE.gz")
UA_female   = fread("../data/CKDGen_ChrX_sumStat_UA_FEMALE.gz")
UA_all      = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")

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
head(eGFR)
eGFR = eGFR[!duplicated(rsID),]
eGFR[,table(phenotype)]

UA<-rbind(UA_male,UA_female,UA_all)
setorder(UA,"P")
head(UA)
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
myTab = fread("../results/01_Locus_Definitions.txt")
myTab

stopifnot(sum(is.element(myTab$rsID,plotData$SNP))==22)
plotData[SNP %in% myTab$rsID]
plotData[,matchID := paste(SNP,phenotype,sep="::")]
myTab[,matchID := paste(rsID,phenotype,sep="::")]

matched1 = match(plotData$matchID,myTab$matchID)
plotData[,region_num := myTab[matched1,region]]
plotData[,region_chr := paste("region",myTab[matched1,region])]

plotData[,candidateGene := paste("region",region_num)]
plotData[is.na(region_num),candidateGene := NA]

# plotData[,novelty := F]

#' # Plot ####
#' ***
mytitle = paste0("Miami Plot; top: eGFR, buttom: UA") 

ymaxpar1 = ifelse(plotData[flag=="top",max(logP,na.rm=T)] <7, 8,plotData[flag=="top",max(logP,na.rm=T)]+1)
ymaxpar2 = ifelse(plotData[flag=="bottom",max(logP,na.rm=T)] <7, 8,plotData[flag=="bottom",max(logP,na.rm=T)]+1)

plot = miamiPlot(x=plotData,
                 ymax = ymaxpar1,
                 ymin = -ymaxpar2,
                 title = mytitle,
                 xlabel = "",
                 ylabel=expression(paste("UA: ",log[10](p),"; eGFR: ",-log[10](p))),
                 hline1=-log10(5e-8),hline2=log10(5e-8),
                 sugline1=-log10(1e-6),sugline2=log10(1e-6),
                 highlight=T, diffsize = T,num_breaks_y=10,
                 plotGenes=T,
                 out_name="../results/Fig1_MiamiPlot.pdf",
                 returnObject = T,
                 overall_max = 40, overall_min = -40)

plot

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


