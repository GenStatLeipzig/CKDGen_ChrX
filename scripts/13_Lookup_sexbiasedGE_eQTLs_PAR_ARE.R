#' ---
#' title: "Lookup of PAR, ARE, sex biased eQTLs and sex biased GE"
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
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(projectpath_main)

#' # Check PAR ####
#' ***
#' There should be no SNPs in PAR1 and PAR2 (already filtered in data generating script, first step!)
#' 
#' In hg19 (GRCh37), taken from https://www.ncbi.nlm.nih.gov/grc/human
#' 
#' - PAR1: from base pair 60,001 to	2,699,520
#' - PAR2: from base pair 154,931,044 to 155,260,560
#' 
data1 = fread("../data/CKDGen_ChrX_sumStat_eGFR_ALL.gz")
table(data1$invalid_assoc)

PAR1_end = 2699520
PAR2_start = 154931044

dim(data1[position <= PAR1_end,])[1]
dim(data1[position >= PAR2_start,])[1]

data1[position <= PAR1_end,c(2,8,10,15:17,20)]
data1[position >= PAR2_start,c(2,8,10,15:17,20)]

#' ARE ####
#' ***

ARETab = data.table(read_xlsx("../../../06_lit/2016_Wilson_AndrogenResponseElements_SupTables.xlsx",sheet=4))
names(ARETab)
names(ARETab) = c("chr","start","stop","Fimo_pval","ARE_seq","Tier","GenomLoc","gene","distance","tba","regulation")
ARETab = ARETab[chr=="chrX",]
ARETab[,tba:=NULL]
ARETab[,length := stop-start]
ARETab[,table(length)]
ARE_sum = sum(ARETab$length)
length_nPAR = max(data1$position) - min(data1$position)
ARE_sum
length_nPAR
ARE_sum/length_nPAR
ARE_sum/length_nPAR * dim(data1[invalid_assoc==F,])[1]
dim(ARETab)[1]/length_nPAR * dim(data1[invalid_assoc==F,])[1]

data1[,ARE:=""]

for(i in 1:dim(data1)[1]){
  #i=1
  pos_i = data1[i,position]
  dummy = copy(ARETab)
  dummy = dummy[start>=pos_i,]
  dummy = dummy[stop<=pos_i,]
  if(dim(dummy)[1]==0){
    data1[i,ARE:="no"]
  }else{
    data1[i,ARE:=dummy$Tier[1]]
  }
}

table(data1$ARE)

#' Result: not one of all analyzed 325,770 SNPs is within an ARE!
#' 
#' This assumes that Wilson et al. used hg19 as well!
#' 
#' # SB eQTLs ####
#' ***
sb_eQTLs_complete = fread(paste0(path_GTExv8_sexStrat,"GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt"))
sb_eQTLs_complete = sb_eQTLs_complete[grepl("chrX",variant_id),]
sb_eQTLs = sb_eQTLs_complete[pvals.corrected<0.05,]
sb_eQTLs[,.N,by=hugo_gene_id]

#' candidate genes copy pasted from MainTable1_ms.xlsx (05/12/2022)
candidateGenes = c("FAM9B","CDKL5","CDK16", "USP11","AR", "EDA2R","BRWD3","TSPAN6","DRP2",
                   "ARMCX2", "ARMCX4","MORF4L2", "TCEAL3","CLDN2","ACSL4","SLC25A5","DCAF12L1",
                   "STK26","HPRT1","DUSP9", "FAM58A","AR", "EDA2R","CITED1", "PIN4","ARMCX2", 
                   "ARMCX4","MORF4L2", "TCEAL3","DCAF12L1","HPRT1","DUSP9", "FAM58A")
candidateGenes = unique(candidateGenes)

sb_eQTLs = sb_eQTLs[hugo_gene_id %in%  candidateGenes,]
setorder(sb_eQTLs,variant_id)
sb_eQTLs[,c(1:12)]

tab13_0 = sb_eQTLs[,.N,by=c("hugo_gene_id")]
tab13_1 = sb_eQTLs[Tissue == "Kidney_Cortex",.N,by=c("hugo_gene_id")]

#' # SB gene expression ####
#' ***
sb_GE = fread(paste0(path_GTExv8_sexStrat,"GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt"))
sb_GE = sb_GE[gene %in% sb_eQTLs_complete$ensembl_gene_id]

matched = match(sb_GE$gene,sb_eQTLs_complete$ensembl_gene_id)
sb_GE[, gene_id := sb_eQTLs_complete[matched,hugo_gene_id]]

sb_GE = sb_GE[gene_id %in% candidateGenes,]
sb_GE[,.N,by=gene_id]
sb_GE[,sign:= sign(effsize)]
sb_GE[sign == 1, higherIn := "females"]
sb_GE[sign == -1, higherIn := "males"]
tab13_1_1 = sb_GE[higherIn == "females",.N,by=c("gene_id","higherIn")]
tab13_1_2 = sb_GE[higherIn == "males",.N,by=c("gene_id","higherIn")]
tab13_2 = sb_GE[tissue == "Kidney_Cortex",.N,by=c("gene_id","higherIn")]

#' # Combine SB stuff ####
#' ***
myGenes = unique(c(tab13_0$hugo_gene_id,tab13_1_1$gene_id,tab13_1_2$gene_id,tab13_2$gene_id))

tab13 = data.table(genes = myGenes)
matched1 = match(tab13$genes,tab13_0$hugo_gene_id)
tab13[,sex_biased_eQTLs := tab13_0[matched1,N]]

matched2 = match(tab13$genes,tab13_1_1$gene_id)
tab13[,sex_biased_GE_higherInFemales := tab13_1_1[matched2,N]]
matched3 = match(tab13$genes,tab13_1_2$gene_id)
tab13[,sex_biased_GE_higherInMales := tab13_1_2[matched3,N]]
matched4 = match(tab13$genes,tab13_2$gene_id)
tab13[,kidneyCortex := tab13_2[matched4,higherIn]]

tab13[is.na(sex_biased_eQTLs),sex_biased_eQTLs := 0]
tab13[is.na(sex_biased_GE_higherInFemales),sex_biased_GE_higherInFemales := 0]
tab13[is.na(sex_biased_GE_higherInMales),sex_biased_GE_higherInMales := 0]

tab13[genes == "STK26",genes := "MST4"]
fwrite(tab13,file="../results/13_lookup_sexbiased_GE_eQTLs.txt")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
