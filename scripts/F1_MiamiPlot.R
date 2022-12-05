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

setwd(paste0(projectpath,"scripts/"))

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
#' I use Markus annotation table ...
#'  
#' I add the information about candidate genes, novelty, and sex-specificity manually, as there is no other way right now. 
#'
#' * Candidate genes: taken from excel sheet table1_2022-10-15.xlsx (stored locally - not pretty!)
#' * Novelty: taken from excel sheet table1_2022-10-15.xlsx (stored locally - not pretty!)
#' * IA 

myTab = data.table(read_excel("../tables/MainTable1_ms.xlsx",sheet=1))
myTab[,...23:=NULL]
myTab[,...24:=NULL]
myTab[,...28:=NULL]

myTab[,CandidateGenes := genes]
myTab[,CandidateGenes := gsub(", ",",\n",CandidateGenes)]

stopifnot(sum(is.element(myTab$indexSNP,plotData$SNP))==23)

plotData[SNP %in% myTab$indexSNP]
plotData[,matchID := paste(SNP,phenotype,sep="::")]
myTab[,matchID := paste(indexSNP,bestSetting,sep="::")]

matched1 = match(plotData$matchID,myTab$matchID)
table(is.na(matched1))
plotData[,region_num := myTab[matched1,region]]
plotData[,region_chr := paste("region",myTab[matched1,region])]

plotData[,candidateGene :=myTab[matched1,CandidateGenes]]
plotData[,NoveltySexIA :=myTab[matched1,novelty]]
table(plotData$NoveltySexIA)
plotData[!is.na(NoveltySexIA)]
plotData[grepl("sex",NoveltySexIA)]
plotData[grepl("sex",NoveltySexIA),sexIA:=c("male","male","female","male","female","male")]

save(plotData,file = "../results/F1_MiamiPlot_PlotData.RData")

plotData[region_num %in% c(7,8) & is.na(sexIA), NoveltySexIA:=NA]

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
                 out_name="../figures/MainFigure1_MiamiPlot.pdf",
                 returnObject = T,
                 overall_max = 40, overall_min = -40,useBasePosition = T)

plot1


tiff(filename = "../figures/MainFigure1_MiamiPlot.tiff", 
     width = 4750, height = 3000, res=300, compression = 'lzw')
plot1
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


