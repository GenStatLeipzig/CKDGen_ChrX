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

source("../SourceFile_aman.R")
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

myTab[,CandidateGenes := `candidate genes`]
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
plotData[grepl("sex",NoveltySexIA),sexIA:=c("male","male","female","male","male","female")]

save(plotData,file = "../results/F1_MiamiPlot_PlotData.RData")
load("../results/F1_MiamiPlot_PlotData.RData")

#' # Plot ####
#' ***
#mytitle = paste0("Miami Plot; top: eGFR, buttom: UA") 

ymaxpar1 = ifelse(plotData[flag=="top",max(logP,na.rm=T)] <7, 8,plotData[flag=="top",max(logP,na.rm=T)]+1)
ymaxpar2 = ifelse(plotData[flag=="bottom",max(logP,na.rm=T)] <7, 8,plotData[flag=="bottom",max(logP,na.rm=T)]+1)

plot1 = miamiPlot(x=copy(plotData),
                 ymax = ymaxpar1,
                 ymin = -ymaxpar2,
                 title = "",
                 xlabel = "",
                 ylabel=expression(paste("UA: ",log[10](p),"                              eGFR: ",-log[10](p))),
                 hline1=-log10(5e-8),hline2=log10(5e-8),
                 sugline1=NULL,sugline2=NULL,
                 highlight=T, diffsize = T,num_breaks_y=10,
                 plotGenes=T,
                 out_name="../figures/MainFigure1_MiamiPlot_231121.pdf",
                 returnObject = T,
                 overall_max = 40, overall_min = -40,useBasePosition = T)

# raw miami plot
plot1

# add cytoband
cyto = copy(myTab)
cyto = cyto[!duplicated(cytoband)]
cyto = cyto[,c(2,5)]
setorder(cyto,position)
cyto[,cytoband := gsub("X","",cytoband)]
cyto =rbind(cyto,cyto[c(8,14)])
cyto[16,cytoband := "q22.1-3"]
cyto[16,position := 103029136]
cyto[17,cytoband := "q26.2-3"]
cyto[17,position := 132524288]
cyto = cyto[c(1:6,10:12,15:17)]
setorder(cyto,position)

plot2 <- plot1 + geom_text(data=cyto,aes(x=position, y=0, label=cytoband), colour = "black",size=4)
plot2

# add gene names
genes = copy(plotData)
genes = genes[!is.na(candidateGene)]

plot3 <- plot2 + 
  
  # top: novel & female-specific hits
  geom_label_repel(data = subset(genes, NoveltySexIA=="yes_sexia" & sexIA=="female"),
                   aes(x=BP, y=logP, label = candidateGene),
                   ylim = c(31,40),fontface = 'bold.italic',color = "#D73027") + 
  
  # top: known & female-specific hits
  geom_text_repel(data = subset(genes, NoveltySexIA=="no_sexia" & sexIA=="female"),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(30,40), xlim = c(Inf, 135000000), fontface = 'bold.italic',color = "#D73027") + 
  
  # top: known & male-specific hits
  geom_text_repel(data = subset(genes, NoveltySexIA=="no_sexia"& sexIA=="male" & flag == "top"),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(33,40),fontface = 'bold.italic',color = "#4575B4") + 
  
  # bottom: known & male-specific hits
  geom_text_repel(data = subset(genes, NoveltySexIA=="no_sexia"& sexIA=="male" & flag == "bottom"),
                  aes(x=BP, y=-logP, label = candidateGene),
                  ylim = c(-20,-15),fontface = 'bold.italic',color = "#4575B4") + 
  
  # top: novel & sex-unspecific hits
  geom_label_repel(data = subset(genes, NoveltySexIA=="yes" & flag == "top" & candidateGene != "TSPAN6"),
                   aes(x=BP, y=logP, label = candidateGene),
                   ylim = c(25,Inf),xlim = c(110000000,Inf),fontface = 'bold') + 
  geom_label_repel(data = subset(genes, NoveltySexIA=="yes" & flag == "top" & candidateGene == "TSPAN6"),
                   aes(x=BP, y=logP, label = candidateGene),
                   ylim = c(25,Inf),xlim = c(-Inf,99100000),fontface = 'bold') + 
  
  # top: known & sex-unspecific hits
  geom_text_repel(data = subset(genes, NoveltySexIA=="no" & flag == "top" & 
                                  candidateGene %nin% c( "ARMCX2,\nARMCX4" , "MORF4L2,\nTCEAL3","DCAF12L1")),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(15,Inf)) + 
  geom_text_repel(data = subset(genes, NoveltySexIA=="no" & flag == "top" & candidateGene == "ARMCX2,\nARMCX4"),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(12,25),xlim =c(-Inf,98000000) ) + 
  geom_text_repel(data = subset(genes, NoveltySexIA=="no" & flag == "top" & candidateGene == "MORF4L2,\nTCEAL3"),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(7.3,15),xlim =c(-Inf,100000000) ) + 
  geom_text_repel(data = subset(genes, NoveltySexIA=="no" & flag == "top" & candidateGene == "DCAF12L1"),
                  aes(x=BP, y=logP, label = candidateGene),
                  ylim = c(7.3,Inf),nudge_x= -10000, nudge_y= 3 ) + 
  
  # bottom: known & sex-unspecific hits
  geom_text_repel(data = subset(genes, NoveltySexIA=="no" & flag == "bottom"),
                  aes(x=BP, y=-logP, label = candidateGene),
                  ylim = c(-18,-13))

plot3

# add legend for labels
dummy = data.table(lab = c("novel \nloci","sex \ninteraction"),
                   yaxis = c(-11,-15),
                   xaxis = c(35092665,35092665))

plot4 = plot3 + annotate("rect", 
                           xmin = 10012628, xmax = 41482665,
                           ymin = -25, ymax = -8,
                           fill = "#F2F2F2")+   
  geom_text(data = subset(dummy, lab=="sex \ninteraction"),
            aes(x=xaxis, y=yaxis, label = lab),
            fontface = 'bold.italic') +   
  geom_label(data = subset(dummy, lab!="sex \ninteraction"),  
             aes(x=xaxis, y=yaxis, label = lab),
             fontface = 'bold') 

plot4

# save plot
message("Create PDF")

pdf_from_png(code2parseOrPlot = plot4, 
             pdf_filename = "../figures/MainFigure1_MiamiPlot_231121.pdf",
             weite = 12,
             laenge = 8,
             einheiten = "in",
             resolution = 150)


tiff(filename = "../figures/MainFigure1_MiamiPlot_231121.tiff", 
     width = 4800, height = 2440, res=300, compression = 'lzw')
plot4
dev.off()

pdf(file = "../figures/Figure1.pdf", width = 12, height = 8)
plot4
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")


