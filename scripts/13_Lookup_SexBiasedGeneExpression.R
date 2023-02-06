#' ---
#' title: "Lookup of sex-biased gene expression"
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
#' I want to check the candidate genes for sex-biased gene expression and sex-biased eQTLs.
#' 
#' Data taken from [Oliva et al.](https://www.science.org/doi/10.1126/science.aba3066)
#' 
#' # Initialize ####
#' ***
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_angmar.R")

setwd(projectpath_main)

#' Load data ####
#' ***
sexBiasedGenes = fread(paste0(path_GTExv8_sexStrat,"/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt"))
sexBiasedGenes
sexBiasedeQTLs = fread(paste0(path_GTExv8_sexStrat,"/GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt"))
sexBiasedeQTLs

load("../results/07_a_usedGenes.RData")
myGenTab

myTab = data.table(read_excel("../tables/MainTable1_ms.xlsx",sheet=1))
myTab[,...23:=NULL]
myTab[,...24:=NULL]
myTab[,...28:=NULL]
myTab

#' # Filter for genes tested in coloc ####
#' ***
sexBiasedGenes[,gene2 := substr(gene,1,15)]
sexBiasedeQTLs[,gene2 := substr(ensembl_gene_id,1,15)]

sexBiasedGenes = sexBiasedGenes[gene2 %in% myGenTab$ensg]
sexBiasedeQTLs = sexBiasedeQTLs[gene2 %in% myGenTab$ensg]

matched = match(sexBiasedGenes$gene2,myGenTab$ensg)
sexBiasedGenes[,gene3 := myGenTab[matched,genename]]
matched = match(sexBiasedeQTLs$gene2,myGenTab$ensg)
sexBiasedeQTLs[,gene3 := myGenTab[matched,genename]]

sexBiasedeQTLs_sig = copy(sexBiasedeQTLs)
sexBiasedeQTLs_sig = sexBiasedeQTLs_sig[pvals.corrected <0.05]

x1 = sexBiasedGenes[,.N,by=gene3]
hist(x1$N)
x1[N>20]

y1 = sexBiasedeQTLs_sig[,.N,by=gene3]
hist(y1$N)
y1[N>2]

#' # Filter for candidate genes ####
#' ***
candidateGenes = myTab$`candidate genes`
candidateGenes = unlist(strsplit(candidateGenes,", "))
candidateGenes = unique(candidateGenes)

sexBiasedGenes2 = copy(sexBiasedGenes)
sexBiasedGenes2 = sexBiasedGenes2[gene3 %in% candidateGenes]

sexBiasedeQTLs_sig2 = copy(sexBiasedeQTLs_sig)
sexBiasedeQTLs_sig2 = sexBiasedeQTLs_sig2[gene3 %in% candidateGenes]

x2 = sexBiasedGenes2[,.N,by=gene3]
hist(x2$N)
setorder(x2,-N)
x2

y2 = sexBiasedeQTLs_sig2[,.N,by=gene3]
hist(y2$N)
setorder(y2,-N)
y2

#' # Filter for candidate genes with sexIA ####
#' ***
candidateGenes2 = myTab[grepl("sex",novelty),`candidate genes`]
candidateGenes2 = unlist(strsplit(candidateGenes2,", "))
candidateGenes2 = unique(candidateGenes2)
candidateGenes2

sexBiasedGenes3 = copy(sexBiasedGenes)
sexBiasedGenes3 = sexBiasedGenes3[gene3 %in% candidateGenes2]

sexBiasedeQTLs_sig3 = copy(sexBiasedeQTLs_sig)
sexBiasedeQTLs_sig3 = sexBiasedeQTLs_sig3[gene3 %in% candidateGenes2]

sexBiasedGenes3[,sign :=sign(effsize)]
x3 = sexBiasedGenes3[,.N,by=c("gene3","sign")]
setorder(x3,-N)
x3

y3 = sexBiasedeQTLs_sig3[,.N,by=gene3]
setorder(y3,-N)
y3
sexBiasedeQTLs_sig3

sexBiasedGenes3[tissue == "Kidney_Cortex",]

#' # Summary ####
#' ***
#' Of all 189 genes tested in coloc, there are 106 with sex-biased gene expression. 
#' 
#' Of all 23 candidate genes, there are 17 with sex-biased gene expression. 
#' 
#' Of the 7 genes with lead SNP sex interactions, there are 5 with sex-biased gene expression. 
#' 
#' Regarding eQTLs, comparison is a bit tricky, as only the MASH independent SNPs are reported (not all cis-SNPs per gene). However, 12 of 23 and 2 of 7 genes are reported with at least one sex-biased eQTL. 
#' 
#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
