#' ---
#' title: "Plot Pvalues of top hits in all phenotypes"
#' subtitle: "A Heatmap Plot"
#' author: "Katrin Horn"
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
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
rm(list=ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(projectpath_main)

#' # Preparations ####
#' ***
#' ##  ####
#' 
#' ## Load data of top hits 
load("../results/04_lookup_TopHits_allOtherTraits.RData")

#' ## Build data object for plotting
toPlot = copy(WideTable)
toPlot = toPlot[grepl("eGFR_", topPheno), ]
cols2keep = c("rsID", "P_eGFR_ALL", "P_eGFR_FEMALE", "P_eGFR_MALE","P_oneSided_BUN_ALL", "P_oneSided_BUN_FEMALE", "P_oneSided_BUN_MALE", 
              "P_oneSided_CKD_ALL", "P_oneSided_CKD_FEMALE", "P_oneSided_CKD_MALE")
colsOut = setdiff(colnames(toPlot), cols2keep)
toPlot[, get("colsOut") := NULL]
setcolorder(toPlot, cols2keep)

#remove "P_" from column names
setnames(toPlot, c("rsID", "P_eGFR_ALL", "P_eGFR_FEMALE", "P_eGFR_MALE", "P_oneSided_BUN_ALL", "P_oneSided_BUN_FEMALE", 
                   "P_oneSided_BUN_MALE", "P_oneSided_CKD_ALL", "P_oneSided_CKD_FEMALE", "P_oneSided_CKD_MALE"), 
         c("SNP", "eGFR_ALL", "eGFR_FEMALE", "eGFR_MALE", "BUN_ALL", "BUN_FEMALE", "BUN_MALE", "CKD_ALL", "CKD_FEMALE", "CKD_MALE"))

#short rsIDs
shorts = apply(toPlot, 1, function(x) return(unlist(strsplit(x[1], split =":"))[1]))
shorts[shorts == "chr23"] = "chr23:152898260"
toPlot[, SNP := shorts]

#get data in long format
toPlot2 = melt(toPlot, id.vars = 1)

#change scale (four different levels of significance)
toPlot2[, significance := "none"]
toPlot2[-log10(value) >= -log10(0.01), significance := "nominal"]
toPlot2[-log10(value) >= -log10(1*10^-6), significance := "suggestive"]
toPlot2[-log10(value) >= -log10(5*10^-8), significance := "genome-wide"]

#change to data.frame and set factor order
toPlot2 = as.data.frame(toPlot2)
toPlot2$significance = factor(toPlot2$significance, levels=c("genome-wide", "suggestive", "nominal", "none"))

#same order of SNPs as in main table 1
toPlot2$nr.SNP = NA
toPlot2$nr.SNP[toPlot2$SNP=="rs139036121"] = 16
toPlot2$nr.SNP[toPlot2$SNP=="rs5909184"] = 15
toPlot2$nr.SNP[toPlot2$SNP=="rs72616719"] = 14
toPlot2$nr.SNP[toPlot2$SNP=="rs189618857"] = 13
toPlot2$nr.SNP[toPlot2$SNP=="rs2063579"] = 12
toPlot2$nr.SNP[toPlot2$SNP=="rs1802288"] = 11
toPlot2$nr.SNP[toPlot2$SNP=="rs3850318"] = 10
toPlot2$nr.SNP[toPlot2$SNP=="rs149995096"] = 9
toPlot2$nr.SNP[toPlot2$SNP=="rs11092455"] = 8
toPlot2$nr.SNP[toPlot2$SNP=="rs181497961"] = 7
toPlot2$nr.SNP[toPlot2$SNP=="rs5942852"] = 6
toPlot2$nr.SNP[toPlot2$SNP=="rs16275"] = 5
toPlot2$nr.SNP[toPlot2$SNP=="rs5931180"] = 4
toPlot2$nr.SNP[toPlot2$SNP=="rs5933079"] = 3
toPlot2$nr.SNP[toPlot2$SNP=="rs5933443"] = 2
toPlot2$nr.SNP[toPlot2$SNP=="chr23:152898260"] = 1
toPlot2$SNP = factor(toPlot2$SNP, levels=(toPlot2$SNP)[order(unique(toPlot2$nr.SNP))])
head(toPlot2)


#' ## Plot Heatmap
col_fun = colorRamp2(c(0, 3), c("red", "white"))
myColors = col_fun(seq(0, 3)) 

p = ggplot(data = toPlot2, aes(x=variable, y=SNP, fill=significance)) + geom_tile(colour = "grey50") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 7, angle = 50, hjust = 1), axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7)) +
  scale_fill_manual(values = myColors) + theme(legend.position = "top") + 
  theme(legend.text = element_text(size = 7), legend.title = element_text(size = 7)) + 
  theme(plot.margin = margin(0,0.35,0,0.1, "cm")) + theme(legend.key.size = unit(0.3, "cm"), legend.key.width = unit(0.5,"cm"))
p

tiff(file="../figures/MainFigure2_Heatmap_eGFR_hits.tiff", width = 1400, height = 800, res = 300, compression = 'lzw')
p
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

