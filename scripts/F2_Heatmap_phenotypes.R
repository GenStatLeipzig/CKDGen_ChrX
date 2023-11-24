#' ---
#' title: "Plot Pvalues of top hits in phenotypes eGFR, UA, BUN, CKD, UACR and MA"
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
#' This is a script to visualise the connection between the phenotypes eGFR, UA, BUN, CKD, UACR and MA. This is done with a heatmap plot of 
#' p values for the 23 index SNPs of the aforementioned phenotypes. The p values for eGFR and UA are two-sided from the Meta-GWAS results. 
#' For all other phenotypes one-sided p values are used. The best setting (ALL, MALE or FEMALE) is chosen for all phenotypes.
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
table(WideTable[, topPheno])

#' ## Build data object for plotting
#ALL
ALL = copy(WideTable)
ALL = ALL[grepl("_ALL", topPheno), ]
cols2keep = c("region", "Cytoband", "rsID", "beta_eGFR_ALL", "P_eGFR_ALL", "beta_UA_ALL", "P_UA_ALL", "beta_BUN_ALL", "P_oneSided_BUN_ALL", 
              "beta_CKD_ALL", "P_oneSided_CKD_ALL", "beta_UACR_ALL", "P_oneSided_UACR_ALL", "beta_MA_ALL", "P_oneSided_MA_ALL")
colsOut = setdiff(colnames(ALL), cols2keep)
ALL[, get("colsOut") := NULL]
setcolorder(ALL, cols2keep)
setnames(ALL, cols2keep, c("Locus", "Cytoband", "rsID", "beta_eGFR", "eGFR", "beta_UA", "UA", "beta_BUN", "BUN", "beta_CKD", "CKD", "beta_UACR", "UACR", "beta_MA", "MA"))
ALL[, Locus2 := as.character(Locus)]
ALL[Locus2 == "7", Locus2 := "7A"]
ALL[, Cytoband := paste0(Locus2, ": ALL - ", Cytoband)]
ALL[, Locus2 := NULL]

#FEMALE
FEMALE = copy(WideTable)
FEMALE = FEMALE[grepl("_FEMALE", topPheno), ]
cols2keep = c("region", "Cytoband", "rsID", "beta_eGFR_FEMALE", "P_eGFR_FEMALE", "beta_UA_FEMALE", "P_UA_FEMALE", "beta_BUN_FEMALE", "P_oneSided_BUN_FEMALE", 
              "beta_CKD_FEMALE", "P_oneSided_CKD_FEMALE", "beta_UACR_FEMALE", "P_oneSided_UACR_FEMALE", "beta_MA_FEMALE", "P_oneSided_MA_FEMALE")
colsOut = setdiff(colnames(FEMALE), cols2keep)
FEMALE[, get("colsOut") := NULL]
setcolorder(FEMALE, cols2keep)
setnames(FEMALE, cols2keep, c("Locus", "Cytoband", "rsID", "beta_eGFR", "eGFR", "beta_UA", "UA", "beta_BUN", "BUN", "beta_CKD", "CKD", "beta_UACR", "UACR", "beta_MA", "MA"))
FEMALE[, Locus2 := as.character(Locus)]
FEMALE[Locus2 == 7, Locus2 := "7B"]
FEMALE[, Cytoband := paste0(Locus2, ": FEMALE - ", Cytoband)]
FEMALE[, Locus2 := NULL]
  
#MALE
MALE = copy(WideTable)
MALE = MALE[grepl("_MALE", topPheno), ]
cols2keep = c("region", "Cytoband", "rsID", "beta_eGFR_MALE", "P_eGFR_MALE", "beta_UA_MALE", "P_UA_MALE", "beta_BUN_MALE", "P_oneSided_BUN_MALE", 
              "beta_CKD_MALE", "P_oneSided_CKD_MALE", "beta_UACR_MALE", "P_oneSided_UACR_MALE", "beta_MA_MALE", "P_oneSided_MA_MALE")
colsOut = setdiff(colnames(MALE), cols2keep)
MALE[, get("colsOut") := NULL]
setcolorder(MALE, cols2keep)
setnames(MALE, cols2keep, c("Locus", "Cytoband", "rsID", "beta_eGFR", "eGFR", "beta_UA", "UA", "beta_BUN", "BUN", "beta_CKD", "CKD", "beta_UACR", "UACR", "beta_MA", "MA"))
MALE[, Cytoband := paste0(Locus, ": MALE - ", Cytoband)]

#combine to one object
toPlot = rbindlist(list(ALL, FEMALE, MALE), use.names = T)
setkey(toPlot, "Locus")

#short rsIDs
shorts = apply(toPlot, 1, function(x) return(unlist(strsplit(x[3], split =":"))[1]))
shorts[shorts == "chr23"] = "chr23:152898260"
toPlot[, rsID := paste(Cytoband, " - ", shorts)]
toPlot[, Cytoband := NULL]
toPlot[, Locus := NULL]

#get data in long format
toPlot.beta = melt(data = toPlot, id.vars = 1, measure.vars = c("beta_eGFR", "beta_UA", "beta_BUN", "beta_CKD", "beta_UACR", "beta_MA"))
toPlot.p = melt(data = toPlot, id.vars = 1, measure.vars = c("eGFR", "UA", "BUN", "CKD", "UACR", "MA"))
setnames(toPlot.p, c("rsID", "variable", "value"), c("SNP", "phenotype", "Pvalue"))
setnames(toPlot.beta, c("variable", "value"), c("phenotype2", "Beta"))
toPlot2 = cbind(toPlot.p, toPlot.beta)

#a little check
table(toPlot2[, rsID] == toPlot2[, SNP])   #all identical
table(toPlot2[, phenotype], toPlot2[, phenotype2])

#remove columns not needed anymore
toPlot2[, rsID := NULL]
toPlot2[, phenotype2 := NULL]

#change scale (three different levels of significance)
toPlot2[, signif := "missing"]
toPlot2[!is.na(Pvalue), signif := "none"]
toPlot2[-log10(Pvalue) >= -log10(0.05), signif := "nominal"]
toPlot2[-log10(Pvalue) >= -log10(5*10^-8), signif := "genome-wide"]

#define colour depending on significance and Beta
toPlot2[, significance := "missing"]  
toPlot2[signif == "none" & Beta > 0, significance := "not significant, beta > 0"]
toPlot2[signif == "none" & Beta < 0, significance := "not significant, beta < 0"]
toPlot2[signif == "nominal" & Beta > 0, significance := "nominal, beta > 0"]
toPlot2[signif == "nominal" & Beta < 0, significance := "nominal, beta < 0"]
toPlot2[signif == "genome-wide" & Beta > 0, significance := "genome-wide, beta > 0"]
toPlot2[signif == "genome-wide" & Beta < 0, significance := "genome-wide, beta < 0"]
table(toPlot2[, significance])

#change to data.frame and set factor order
toPlot2 = as.data.frame(toPlot2)
toPlot2$significance = factor(toPlot2$significance, levels=c("genome-wide, beta > 0", "nominal, beta > 0", "not significant, beta > 0", 
                                                             "genome-wide, beta < 0", "nominal, beta < 0", "not significant, beta < 0", 
                                                             "missing"))
toPlot2$phenotype = factor(toPlot2$phenotype, levels=c("eGFR", "UA", "BUN", "CKD", "UACR", "MA"))

#add main phenotype of top hit to data object
toPlot2$topPheno = "eGFR"
UA_SNPs = c("16: ALL - Xq12  -  rs6625094", "17: ALL - Xq13.1  -  rs34687188", "18: ALL - Xq22.1  -  rs34884874", 
            "19: ALL - Xq22.1  -  rs34815154", "20: ALL - Xq25  -  rs112708523", "21: ALL - Xq26.3  -  rs202138804", 
            "22: ALL - Xq28  -  rs4328011")
toPlot2$topPheno[is.element(toPlot2$SNP, UA_SNPs)] = "UA"

#same order of SNPs as in main table 1
toPlot2$nr.SNP = NA
toPlot2$nr.SNP[grep("rs139036121", toPlot2$SNP)] = 23
toPlot2$nr.SNP[grep("rs5909184", toPlot2$SNP)] = 22
toPlot2$nr.SNP[grep("rs72616719", toPlot2$SNP)] = 21
toPlot2$nr.SNP[grep("rs189618857", toPlot2$SNP)] = 20
toPlot2$nr.SNP[grep("rs2063579", toPlot2$SNP)] = 19
toPlot2$nr.SNP[grep("rs1802288", toPlot2$SNP)] = 18
toPlot2$nr.SNP[grep("rs3850318", toPlot2$SNP)] = 17
toPlot2$nr.SNP[grep("rs149995096", toPlot2$SNP)] = 16
toPlot2$nr.SNP[grep("rs11092455", toPlot2$SNP)] = 15
toPlot2$nr.SNP[grep("rs181497961", toPlot2$SNP)] = 14
toPlot2$nr.SNP[grep("rs5942852", toPlot2$SNP)] = 13
toPlot2$nr.SNP[grep("rs16275", toPlot2$SNP)] = 12
toPlot2$nr.SNP[grep("rs5931180", toPlot2$SNP)] = 11
toPlot2$nr.SNP[grep("rs5933079", toPlot2$SNP)] = 10
toPlot2$nr.SNP[grep("rs5933443", toPlot2$SNP)] = 9
toPlot2$nr.SNP[grep("chr23:152898260", toPlot2$SNP)] = 8
toPlot2$nr.SNP[grep("rs6625094", toPlot2$SNP)] = 7
toPlot2$nr.SNP[grep("rs34687188", toPlot2$SNP)] = 6
toPlot2$nr.SNP[grep("rs34884874", toPlot2$SNP)] = 5
toPlot2$nr.SNP[grep("rs34815154", toPlot2$SNP)] = 4
toPlot2$nr.SNP[grep("rs112708523", toPlot2$SNP)] = 3
toPlot2$nr.SNP[grep("rs202138804", toPlot2$SNP)] = 2
toPlot2$nr.SNP[grep("rs4328011", toPlot2$SNP)] = 1
toPlot2$SNP = factor(toPlot2$SNP, levels=(toPlot2$SNP)[order(unique(toPlot2$nr.SNP))])
head(toPlot2)

#' ## Plot Heatmap
col_red = colorRamp2(c(0, 3), c("coral1", "white"))
col_blue = colorRamp2(c(0, 3), c("dodgerblue", "white"))
MyColors = c(col_red(seq(0, 2)), col_blue(seq(0, 3)))

p = ggplot(data = toPlot2, aes(x=phenotype, y=SNP, fill=significance)) + geom_tile(colour = "grey40") + 
  facet_grid(rows = "topPheno", scales = "free_y", space = "free_y") +
  theme_bw() +  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12, angle = 50, hjust = 1), 
        axis.ticks.x = element_blank(), axis.title.y = element_text(size = 12)) + 
  theme(axis.text.y = element_text(size = 12), strip.text = element_text(size = 10)) +
  scale_fill_manual(values = MyColors) + theme(legend.position = "right") + 
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 12)) + 
  theme(plot.margin = margin(0.1, 0.15, 0.1 ,0.1, "cm")) + theme(legend.key.size = unit(0.3, "cm"), legend.key.width = unit(0.5,"cm")) +
  coord_cartesian(clip = "off") + theme(plot.tag.position = c(.8, .3), plot.tag = element_text(size = 7))

p

pdf(file="../figures/Figure4.pdf", width = 9, height = 6)
p
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

