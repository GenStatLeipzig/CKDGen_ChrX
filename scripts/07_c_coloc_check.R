#' ---
#' title: "Co-localization Part 3: Check output"
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
#' *Co-localization analysis of gene-expression quantitative trait loci*: 
#' 
#' We tested for overlapping causal variants between kidney trait associations and gene-expression quantitative trait loci (eQTLs).
#' 
#' Used eQTL databases:
#' 
#' * GTEx Analysis V8, 49 tissues, hg38-> liftover using the GTEx SNP Annotation file
#' * NEPTUNE, 2 tissues, hg19
#' 
#' We analyzed all tissues, but laid main focus in the following tissues: 
#' 
#' * GTEx adrenal gland
#' * GTEx kidney cortex
#' * GTEx muscle skeletal
#' * GTEx whole blood
#' * NEPTUNE glomerulus
#' * NEPTUNE tubulointerstitial  
#' 
#' Gene selection: 
#' 
#' * nearby genes (+/- 250 kb)
#' * cis eQTL genes (LD r^2>=0.3)
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(paste0(projectpath,"scripts/"))
source("../helperFunctions/colocPlot.R")

#' # Load data ####
#' ***
load("../temp/07_allGenes.RData")
load("../results/07_b_coloc_eQTLs.RData")
coloc = copy(ColocTable)

tab6 = fread("../results/05_b_Cojo_Select_Results.txt")
tab6[,phenotype := gsub("_all","_ALL",phenotype)]
tab6[,phenotype := gsub("_male","_MALE",phenotype)]
tab6[,phenotype := gsub("_female","_FEMALE",phenotype)]

tab6[,leadSNP := rsID]
tab6[grepl("rs111884516:",rsID),leadSNP := "rs4328011:152898261:G:A"]
tab6[grepl("rs7056552:",rsID),leadSNP := "rs202138804:133799101:AGT:A"]
tab6[grepl("rs111410539:",rsID),leadSNP := "rs181497961:106168067:G:A"]

#' # Check data ####
#' ***
#' ## Check 1: nsnps ###
#' ***
#' number of used SNPs should be >100
check1 = coloc[,min(nsnps)]
check1
coloc[nsnps == min(nsnps)]
hist(coloc$nsnps)
hist(coloc[grepl("Kidney_Cortex_Glomerular",trait2) | grepl("Kidney_Cortex_Tubulointerstitial",trait2),nsnps],
     main = "Histogram NEPTUNE",
     xlab = "NEPTUNE overlap")
hist(coloc[!grepl("Kidney_Cortex_Glomerular",trait2) & !grepl("Kidney_Cortex_Tubulointerstitial",trait2),nsnps],
     main = "Histogram GTEx",
     xlab = "GTEx overlap")

#' ## Check 2: PP.H0 ####
#' ***
#' should not be too high 
check2 = copy(coloc)
check2 = check2[PP.H0.abf>=0.75,]
check2[,dumID := paste(region, trait1,gene,sep="::")]
check2 = check2[!duplicated(dumID),]

#' What could go wrong here? 
#' 
#' * GTEx: range +/- 1 mb around gene
#' * NEPTUNE: range +/- 500 kb around gene
#' * our selection: index SNP +/- 250 kb or eQTL (max dist 1 mb)
#' 
#' * Option 1: index SNP not in eQTL data --> regions overlap, but no signal for CKD traits (might be the case for NEPTUNE data as there are less SNPs available (filtered for MAF>3% per default and lower range than GTEx and our eQTL definition))
#' * Option 2: sex-specific association (check with Main Table 1)
#' 
dummy = check2[trait1 == "eGFR_ALL"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_ALL",]

#' **eGFR_ALL**: all outlier are in NEPTUNE data (see option 1) or not for the ALL phenotype
#' 
dummy = check2[trait1 == "eGFR_MALE"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_MALE",]

#' **eGFR_MALE**: all outlier are in NEPTUNE data (see option 1) or not for the MALE phenotype
#' 
dummy = check2[trait1 == "eGFR_FEMALE",]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_FEMALE",]

#' **eGFR_FEMALE**: all outlier are in NEPTUNE data (see option 1) or not for the FEMALE phenotype
#' 
dummy = check2[trait1 == "UA_ALL"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_ALL",]

data_GWAS = fread("../data/CKDGen_ChrX_sumStat_UA_ALL.gz")
load("../temp/07_coloc/GTEx_v8_filtered_Adipose_Subcutaneous.RData")
data1 = copy(data0)
data1 = data1[cyto == "Xq13.1",]
data_GWAS = data_GWAS[P<1e-6,]
tab6[phenotype=="UA_ALL" & region == 17,]
data_GWAS = data_GWAS[position <= tab6[phenotype=="UA_ALL" & region == 17,bp]+500000,]
data_GWAS = data_GWAS[position >= tab6[phenotype=="UA_ALL" & region == 17,bp]-500000,]
data_GWAS
table(is.element(data_GWAS$position,data1$pos_b37))

#' **UA_ALL**: 
#' 
#' * Region 22: outlier is in NEPTUNE data (see option 1)
#' * Region 17: Index SNP & the two other associated SNPs with suggestive significance are not available in GTEx. I check the GTEx Portal and searched for the rsID but found no matches. 
#' 
dummy = check2[trait1 == "UA_MALE"]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_MALE",]

#' **UA_MALE**: all outlier are in NEPTUNE data (see option 1) or not for the MALE phenotype
#' 
dummy = check2[trait1 == "UA_FEMALE",]
dummy[,c(1:7)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_FEMALE",]

#' **UA_FEMALE**: all outlier are in NEPTUNE data (see option 1) or not for the FEMALE phenotype
#' 
#' ## Check 3: PP.H2 ####
#' ***
#' should not be too high 
check3 = copy(coloc)
check3 = check3[PP.H2.abf>=0.75,]
check3[,dumID := paste(region, trait1,gene,sep="::")]
check3 = check3[!duplicated(dumID),]

dummy = check3[trait1 == "eGFR_ALL"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_ALL",]

#' **eGFR_ALL**: no problems
#' 
dummy = check3[trait1 == "eGFR_MALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_MALE",]

#' **eGFR_MALE**: no problems
#' 
dummy = check3[trait1 == "eGFR_FEMALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "eGFR_FEMALE",]

#' **eGFR_FEMALE**: no problems
#' 
dummy = check3[trait1 == "UA_ALL"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_ALL",]

#' **UA_ALL**: Region 17: Index SNPs not available
#' 
dummy = check3[trait1 == "UA_MALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_MALE",]

#' **UA_MALE**: no problems
#' 
dummy = check3[trait1 == "UA_FEMALE"]
dummy[,c(1:6,9)]
dummy = dummy[!grepl("Kidney_",trait2),]
tab7[genes %in% dummy$gene & phenotype == "UA_FEMALE",]

#' **UA_FEMALE**: no problems
#' 
#' ## Check 4: Coloc ####
#' ***
#' 
check4 = copy(coloc)
check4 = check4[PP.H4.abf>=0.75,]
length(unique(check4$gene))
length(unique(check4$region))

check4[,dumID := paste(gene,trait1,sep="::")]

dumTab6<-dcast(check4,
               formula = dumID ~ trait2,
               value.var = c("PP.H4.abf"),
               sep = "_")
names(dumTab6)
dumTab6[!is.na(`GE in Kidney_Cortex_Tubulointerstitial `),c(1,28)]
dumTab6[!is.na(`GE in Muscle_Skeletal`),c(1,31)]
dumTab6[!is.na(`GE in Whole_Blood`),c(1,45)]

candidates = dumTab6[!is.na(`GE in Kidney_Cortex_Tubulointerstitial `) | !is.na(`GE in Muscle_Skeletal`) | !is.na(`GE in Whole_Blood`) , dumID]
candidates = gsub("::.*","",candidates)
candidates = unique(candidates)
res = tab7[genes %in% candidates,]
res

candidates2 = dumTab6[, dumID]
candidates2 = gsub("::.*","",candidates2)
candidates2 = unique(candidates2)
res2 = tab7[genes %in% candidates2,]
res2

#' # Summary ####
#' ***
#' 
#' **SNP numbers**: all good
#' 
#' **PP H0, PP H2**: all good or explainable
#' 
#' **Co-localization**:
#' 
#' Best candidates in relevant tissues (kidney, muscle, blood):  
#' 
#' * ARMCX2 (UA and eGFR in ALL and MALE) - region 7 and 18
#' * CDKL5 (eGFR in ALL, MALE, and FEMALE) - region 2
#' * SLC25A5 (eGFR in ALL and MALE) - region 11
#' * NDUFB11 (eGFR in ALL) - region 3
#' * TCEAL3 (eGFR and UA in MALE) - region 8 and 19
#' * ACSL4 (eGFR in ALL, MALE and FEMALE) - region 10
#' * MORF4L2 (eGFR and UA in MALE) - region 8 and 19
#' 
#' In total, there are 38 unique genes with positive coloc in at least one tissue, corresponding to 14 of the 22 regions. 
#' 
#' # Plot ####
#' *** 
#' We want to plot the 7 interesting genes in the three relevant tissues for ALL, MALE, FEMALE and two additional genes that have kidney relevant function
#' 
#' * rows: genes_tissues
#' * columns: phenotype_setting
#' 
#' Feedback Markus (30.01.2023): exclude NDUFB11!
#' 
#' --> create a matrix with 8x6 entries!
#' 
plotData1 = copy(ColocTable)
myTab = data.table(read_excel("../tables/MainTable1_ms.xlsx",sheet=1))
myTab[,...23:=NULL]
myTab[,...24:=NULL]
myTab[,...28:=NULL]

candidate2 = myTab$`Coloc-genes`
candidate2 = candidate2[!is.na(candidate2)]
candidate2 = unlist(strsplit(candidate2,", "))
candidate2
candidate2 = candidate2[-4]
plotData1 = plotData1[gene %in% candidate2,]
plotData1[,phenotype := gsub("_.*","",trait1)]
plotData1[,setting := gsub(".*_","",trait1)]

plotData1[,dumID1 := paste(gene,trait2,sep="_")]
plotData1[,dumID2 := paste(trait1,sep="_")]
plotData1[,dumID1 := gsub("GE in ","",dumID1)]
plotData1[,dumID1 := gsub(" ","",dumID1)]

plotData2<-dcast(plotData1,
                 formula = dumID1 ~ dumID2,
                 value.var = c("PP.H4.abf"),
                 sep = "_")
names(plotData2)
M4<-as.matrix(plotData2[,-1])

plotData3<-dcast(plotData1,
                 formula = dumID1 ~ dumID2,
                 value.var = c("PP.H3.abf"),
                 sep = "_")
names(plotData3)
M3<-as.matrix(plotData3[,-1])

x1 = dim(M4)[1]
x2 = dim(M4)[2]

M<-matrix(0,x1,x2)
for (i in 1:x1){
  for (j in 1:x2){
    m4<-M4[i,j]
    m3<-M3[i,j]
    
    if(is.na(m3)==T){
      M[i,j] = NA
    }else if(m4>m3){
      M[i,j]<-m4
    }else{
      M[i,j]<- -m3
    }
  }
}
rownames(M)<-plotData3$dumID1
colnames(M)<-names(plotData3)[-1]

plotData4 = as.data.frame(M)
plotData4 = cbind(plotData3$dumID1,plotData4)
names(plotData4)[1] = "dumID"

dummy<-pmax(plotData4$UA_ALL,plotData4$UA_FEMALE,plotData4$UA_MALE,
            plotData4$eGFR_ALL,plotData4$eGFR_FEMALE,plotData4$eGFR_MALE, na.rm = T)
setDT(plotData4)
plotData4[,maxPPH4 := dummy]
filt = dummy>=0.75 
plotData5 = plotData4[filt,]
plotData5[,gene := gsub("_.*","",dumID)]
setorder(plotData5,-maxPPH4)
plotData6 = copy(plotData5)
plotData6 = plotData6[gene %in% c("CDK16","USP11") & !duplicated(gene) ]

plotData7 = copy(plotData5)
plotData7 = plotData7[gene %nin% c("CDK16","USP11"), ]
plotData7 = plotData7[grepl("Kidney_Cortex_Tubulointerstitial",dumID) | 
                        grepl("Muscle_Skeletal",dumID) | 
                        grepl("Whole_Blood",dumID),]

plotData8 = rbind(plotData7,plotData6)
for(i in 1:8){
  #i=1
  gene = plotData8[i,gene]
  x = grep(gene,myTab$`Coloc-genes`)
  plotData8[i,locus := myTab[x,region]]
}
setorder(plotData8,-locus)

plotData8[,dumID := gsub("_Kidney_Cortex_Tubulointerstitial"," (TI)",dumID)]
plotData8[,dumID := gsub("_Muscle_Skeletal"," (MS)",dumID)]
plotData8[,dumID := gsub("_Whole_Blood"," (WB)",dumID)]
plotData8[,dumID := gsub("_Thyroid"," (Thy)",dumID)]
plotData8[,dumID := gsub("_Stomach"," (Sto)",dumID)]
plotData8[,dumID := paste0(locus,": ",dumID)]

names(plotData8) = gsub("_"," (",names(plotData8))
names(plotData8) = gsub("ALL","ALL)",names(plotData8))
names(plotData8) = gsub("MALE","MALE)",names(plotData8))
plotData8

#update locus label for locus 7: -> 7A:
plotData8$dumID[plotData8$dumID == "7: ARMCX2 (TI)"] = "7A: ARMCX2 (TI)"

#' Change row order manually (ordered by alphabet as default in decast): ALL - MALE - FEMALE
#' 
setDF(plotData8)
colocPlot(x = plotData8[,c(1,5,7,6,2,4,3)],title = "coloc plot")

myPlot1 = colocPlot(x = plotData8[,c(1,5,7,6,2,4,3)],title = "")

tag = format(Sys.time(), "%Y-%m-%d")
tag = gsub("2023-","23-",tag)
tag = gsub("-","",tag)

tiff(filename = paste0("../figures/MainFigure5_ColocPlot_",tag,".tiff"),
     width = 1350, height = 1350, res=250, compression = 'lzw')
myPlot1
dev.off()

pdf_from_png(code2parseOrPlot = myPlot1, 
             pdf_filename = paste0("../figures/MainFigure5_ColocPlot_",tag,".pdf"),
             weite = 8,
             laenge = 8,
             einheiten = "in",
             resolution = 250)

#' # Sessioninfo ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

 