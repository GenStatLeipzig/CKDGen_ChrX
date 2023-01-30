#' ---
#' title: "Replication in HUNT study"
#' subtitle: "CKDGen - Chr X"
#' author: "Janne Pott, Katrin Horn"
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
#' This is a script to look up the top hits of eGFR ALL from our CKDGen meta-analysis in the data of the HUNT study. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(paste0(projectpath,"scripts/"))

#' # Loading data ####
#' ***
#' ## HUNT study ###
#' The HUNT data is stored locally and the path to it is provided in the parameter *path_HUNT*. 
#' 
data_HUNT = fread(path_HUNT)
data_HUNT
data_HUNT[, markerID.1 := paste0("chr", chr, ":", position, ":", coded_all, ":", noncoded_all)]
data_HUNT[, markerID.2 := paste0("chr", chr, ":", position, ":", noncoded_all, ":", coded_all)]

#' ## GWAMA hits ###
#' We use the results of script 01_locus_definition.R. The replication will only be done for the 15 eGFR loci, not the 7 UA loci
#' 
data_CKDGen = fread("../results/01_Locus_Definitions.txt")
data_CKDGen
dummy = strsplit(data_CKDGen[, rsID], split = ":")
dummy = as.data.table(dummy)
dummy = t(dummy)
data_CKDGen[, ourSNPs := paste0("chr23:", dummy[, 2], ":", dummy[, 3], ":", dummy[, 4])]
data_CKDGen = data_CKDGen[grepl("eGFR",phenotype),]

data_CKDGen_indep = fread("../results/05_b_Cojo_Select_Results.txt")
dummy = strsplit(data_CKDGen_indep[, rsID], split = ":")
dummy = as.data.table(dummy)
dummy = t(dummy)
data_CKDGen_indep[, ourSNPs := paste0("chr23:", dummy[, 2], ":", dummy[, 3], ":", dummy[, 4])]
data_CKDGen_indep = data_CKDGen_indep[grepl("eGFR_a",phenotype),]
setorder(data_CKDGen_indep,p)
data_CKDGen_indep = data_CKDGen_indep[!duplicated(rsID),]
setorder(data_CKDGen_indep,region)

#' # Match data ####
#' ***
matched1 = match(data_CKDGen[, ourSNPs], data_HUNT[, markerID.1])
sum(is.na(matched1))   
matched2 = match(data_CKDGen[, ourSNPs], data_HUNT[, markerID.2])
sum(is.na(matched2))  
matched = pmax(matched1, matched2, na.rm = T)
matched
data_CKDGen[is.na(matched),]

matched3 = match(data_CKDGen_indep[, ourSNPs], data_HUNT[, markerID.1])
sum(is.na(matched3))   
matched4 = match(data_CKDGen_indep[, ourSNPs], data_HUNT[, markerID.2])
sum(is.na(matched4))  
matched5 = pmax(matched3, matched4, na.rm = T)
matched5
data_CKDGen_indep[is.na(matched5),]

#' These male specific hits cannot be matched to HUNT! I will remove these SNPs
data_HUNT2 = copy(data_HUNT)
data_HUNT2 = data_HUNT2[matched, ]
data_HUNT2 = data_HUNT2[!is.na(markerID.2), ]
data_HUNT2

data_HUNT3 = copy(data_HUNT)
data_HUNT3 = data_HUNT3[matched5, ]
data_HUNT3 = data_HUNT3[!is.na(markerID.2), ]
data_HUNT3

data_CKDGen2 = copy(data_CKDGen)
data_CKDGen2 = data_CKDGen2[grepl("ALL",phenotype), ]
data_CKDGen2

data_CKDGen3 = copy(data_CKDGen_indep)
data_CKDGen3 = data_CKDGen3[ourSNPs %in% data_HUNT3$markerID.1 | ourSNPs %in% data_HUNT3$markerID.2, ]
data_CKDGen3

#' Check effect alleles
table(data_CKDGen2$effect_allele == data_HUNT2$coded_all)
table(data_CKDGen2$other_allele == data_HUNT2$noncoded_all)
table(data_CKDGen3$refA == data_HUNT3$coded_all)

plot(data_CKDGen2$EAF, data_HUNT2$AF_coded_all)
abline(0,1)
plot(data_CKDGen3$freq, data_HUNT3$AF_coded_all)
abline(0,1)

#' Merge data: We want a data object with SNP, chr, position, effect and other allele, statistics for both CKDGen and HUNT (EAF, beta, se, p-value, log(p-value), sample size, and number of studies in case of CKDGen)
names(data_CKDGen2)[c(6,1,5,7,8,15,16,12,17:20,10,9)]
data_HUNT2[,logP:=-log10(Pval)]
data_HUNT3[,logP:=-log10(Pval)]
names(data_HUNT2)[c(9,6:8,15,10)]

data_Merge2 = cbind(data_CKDGen2[,c(6,1,5,7,8,15,16,12,17:20,10,9)],
                    data_HUNT2[,c(9,6:8,15,10)])
data_Merge2
myNames = c("phenotype","region","SNP","chr","pos","EA","OA","EAF_CKDGen","beta_CKDGen","SE_CKDGen",
            "pval_CKDGen","logP_CKDGen","NSamples_CKDGen","NStudies_CKDGen",
            "EAF_HUNT","beta_HUNT","SE_HUNT","pval_HUNT","logP_HUNT","NSamples_HUNT")
dim(data_Merge2)
length(data_Merge2)
names(data_Merge2) = myNames
data_Merge2

data_CKDGen3[,logP:=-log10(p)]
names(data_CKDGen3)[c(1,2,18,4,6:11,20,12)]

data_Merge3 = cbind(data_CKDGen3[,c(1,2,18,4,6:11,20,12)],
                    data_HUNT3[,c(9,6:8,15,10)])
data_Merge3
myNames = c("phenotype","region","SNP","chr","pos","EA","EAF_CKDGen","beta_CKDGen","SE_CKDGen",
            "pval_CKDGen","logP_CKDGen","NSamples_CKDGen",
            "EAF_HUNT","beta_HUNT","SE_HUNT","pval_HUNT","logP_HUNT","NSamples_HUNT")
dim(data_Merge3)
length(data_Merge3)
names(data_Merge3) = myNames
data_Merge3

#' # Test one-sidedness ####
#' ***
#' For replication, we want the one sided p-value. Hence, we check the given p-value and add the one sided pvalue and log(p-value) for HUNT
#' 
#' ## For index SNPs ###
data_Merge2[,pval_HUNT_oneSided := pnorm(-abs(beta_HUNT/SE_HUNT))]
data_Merge2[,logP_HUNT_oneSided := -log10(pval_HUNT_oneSided)]

data_Merge2[,plot(logP_HUNT,logP_HUNT_oneSided)]
abline(0,1)
data_Merge2[,plot(pval_HUNT,pval_HUNT_oneSided)]
abline(0,1)

data_Merge2[,table(pval_HUNT<0.05,pval_HUNT_oneSided<0.05)]
data_Merge2[pval_HUNT_oneSided<0.05,]
data_Merge2[pval_HUNT_oneSided<0.05 & sign(beta_HUNT)==sign(beta_CKDGen),]

#' Using the one-sided p-value, we can replicate 9 of 13 SNPs in HUNT (eGFR ALL). All SNPs have the same effect direction in CKDGen and HUNT. 
#' 
data_Merge2[pval_HUNT_oneSided<0.05 & sign(beta_HUNT)==sign(beta_CKDGen),replicated := T]
data_Merge2[pval_HUNT_oneSided>0.05 ,replicated := F]

#' ## For independent SNPs ###
data_Merge3[,pval_HUNT_oneSided := pnorm(-abs(beta_HUNT/SE_HUNT))]
data_Merge3[,logP_HUNT_oneSided := -log10(pval_HUNT_oneSided)]

data_Merge3[,plot(logP_HUNT,logP_HUNT_oneSided)]
abline(0,1)
data_Merge3[,plot(pval_HUNT,pval_HUNT_oneSided)]
abline(0,1)

data_Merge3[,table(pval_HUNT<0.05,pval_HUNT_oneSided<0.05)]
data_Merge3[pval_HUNT_oneSided<0.05,]
data_Merge3[pval_HUNT_oneSided<0.05 & sign(beta_HUNT)==sign(beta_CKDGen),]

#' Using the one-sided p-value, we can replicate 11 of 15 SNPs in HUNT. All SNPs have the same effect direction in CKDGen and HUNT. 
#' 
data_Merge3[pval_HUNT_oneSided<0.05 & sign(beta_HUNT)==sign(beta_CKDGen),replicated := T]
data_Merge3[pval_HUNT_oneSided>0.05 ,replicated := F]

#' # Plotting ####
#' ***
plotData1 = copy(data_Merge2)
plotData2 = copy(data_Merge3)

#' ## EAF-EAF Plot ####
myPlot0 = ggplot(plotData1, aes(x=EAF_CKDGen, y=EAF_HUNT,color=replicated)) + 
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size = 3) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="EAF in CKDGen", 
       y = "EAF in HUNT",
       color = "Successful \nreplication")

myPlot0

myPlot0 = ggplot(plotData2, aes(x=EAF_CKDGen, y=EAF_HUNT,color=replicated)) + 
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size = 3) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="EAF in CKDGen", 
       y = "EAF in HUNT",
       color = "Successful \nreplication")

myPlot0

#' ## Beta-Beta Plot ####
c_CKDGen = qnorm(0.975,lower.tail = T)
c_HUNT = qnorm(0.95,lower.tail = T)

myPlot1 = ggplot(plotData1, aes(x=beta_CKDGen, y=beta_HUNT,color=replicated)) + 
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = beta_CKDGen-c_CKDGen*SE_CKDGen, xmax = beta_CKDGen + c_CKDGen*SE_CKDGen)) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="Effect on eGFR ALL in CKDGen", 
       y = "Effect on eGFR ALL in HUNT",
       color = "Successful \nreplication")

myPlot2 = myPlot1 + 
  geom_errorbar(data = subset(plotData1, beta_CKDGen>0),
                aes(ymin = beta_HUNT-c_HUNT*SE_HUNT, ymax = beta_HUNT)) +
  geom_errorbar(data = subset(plotData1, beta_CKDGen<0),
                aes(ymin = beta_HUNT, ymax = beta_HUNT+c_HUNT*SE_HUNT))
myPlot2

plotData2 = plotData2[!grepl("rs111410539:",SNP),]
myPlot3 = ggplot(plotData2, aes(x=beta_CKDGen, y=beta_HUNT,color=replicated)) + 
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed", size=1.25)+
  geom_hline(yintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed", size=1.15)+
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = beta_CKDGen-c_CKDGen*SE_CKDGen, xmax = beta_CKDGen + c_CKDGen*SE_CKDGen)) +
  theme_bw(base_size = 10)+
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=12,face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text = element_text(size=12,face="bold"))+
  labs(x="Effect on eGFR ALL in CKDGen", 
       y = "Effect on eGFR ALL in HUNT",
       color = "Successful \nreplication")

myPlot4 = myPlot3 + 
  geom_errorbar(data = subset(plotData2, beta_CKDGen>0),
                aes(ymin = beta_HUNT-c_HUNT*SE_HUNT, ymax = beta_HUNT)) +
  geom_errorbar(data = subset(plotData2, beta_CKDGen<0),
                aes(ymin = beta_HUNT, ymax = beta_HUNT+c_HUNT*SE_HUNT))
myPlot4

#' # Save results ####
#' ***
#' ## Save plot ####

# tiff(filename = "../figures/SupplementalFigure_BetaBeta_ReplicationHunt_indexSNPs_eGFR_ALL.tiff", 
#      width = 1800, height = 1350, res=250, compression = 'lzw')
# myPlot2
# dev.off()

tiff(filename = "../figures/SupplementalFigure_BetaBeta_ReplicationHunt_indepSNPs_eGFR.tiff", 
     width = 1800, height = 1350, res=250, compression = 'lzw')
myPlot4
dev.off()

#' ## Save data ####
#' 
#' * General information: region, SNP, chr, pos, EA, replicated
#' * CKDGen: phenotype, sample size, EAF, beta, SE, logP, CIs
#' * HUNT: sample size, EAF, beta, SE, logP one sided, CIs one sided
#' 
data_Merge3[, CIlower_CKDGen := round(beta_CKDGen-c_CKDGen*SE_CKDGen, digits = 6)]
data_Merge3[, CIupper_CKDGen := round(beta_CKDGen+c_CKDGen*SE_CKDGen, digits = 6)]
data_Merge3[, CIlower_HUNT := round(beta_HUNT-c_HUNT*SE_HUNT, digits = 6)]
data_Merge3[, CIupper_HUNT := round(beta_HUNT+c_HUNT*SE_HUNT, digits = 6)]
data_Merge3[beta_CKDGen>0,CIupper_HUNT:=Inf]
data_Merge3[beta_CKDGen<0,CIlower_HUNT:=-Inf]
data_Merge3[, ts_power := abs(beta_CKDGen/SE_CKDGen*sqrt(NSamples_HUNT/NSamples_CKDGen))]
data_Merge3[, power := 1-pnorm(c_HUNT,mean=ts_power)]
data_Merge3[, ts_power := NULL]

result = copy(data_Merge3)
myNames = names(result)[c(2:6,21,1,12,7:9,11,22,23,18,13:15,20,24,25,26)]
colsOut<-setdiff(colnames(result),myNames)
result[,get("colsOut"):=NULL]
setcolorder(result,myNames)
result

cor.test(result$beta_CKDGen,result$beta_HUNT)

data.description = data.table(column = names(result),
                              description = c("region according to locus definition (see Table 1)",
                                              "SNP ID according to 1000 Genomes phase 5 version 3",
                                              "chromosome number, 23 corresponds to chromosome X",
                                              "position according to hg19",
                                              "effect allele, also known as coding allele",
                                              "TRUE/FALSE vector indicating successful replication in HUNT",
                                              "subgroup analyzed in CKDGen (ALL, MALE, or FEMALE)",
                                              "number of samples included in CKDGen Meta analysis for eGFR",
                                              "effect allele frequency in CKDGen Meta analysis for eGFR",
                                              "beta estimate in CKDGen Meta analysis for eGFR",
                                              "standard error in CKDGen Meta analysis for eGFR",
                                              "-log10 transformed p-value in CKDGen Meta analysis for eGFR",
                                              "lower bound of the 95% confidence interval in CKDGen Meta analysis",
                                              "upper bound of the 95% confidence interval in CKDGen Meta analysis",
                                              "number of samples included in HUNT study",
                                              "effect allele frequency in HUNT study",
                                              "beta estimate in HUNT study",
                                              "standard error in HUNT study",
                                              "-log10 transformed one-sided p-value in HUNT study",
                                              "lower bound of the one-sided 95% confidence interval in HUNT study",
                                              "upper bound of the one-sided 95% confidence interval in HUNT study",
                                              "power for replication, given the CKDGen statistics and the square-root of HUNT and CKDGen sample size ratio"))

WriteXLS(x = c("result","data.description"), 
         ExcelFileName = "../results/09_replication_HUNT.xlsx",
         SheetNames = c("ReplicationHunt_results","ReplicationHunt_description"),
         AutoFilter=T, 
         BoldHeaderRow=T)
save(data_Merge3, file = "../results/09_replication_HUNT.RData")
save(result,data.description, file = "../results/09_replication_HUNT_summary.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
