#' ---
#' title: "Generate Forest plots of MR-Mega results"
#' subtitle: ""
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
#' ## Load data object with MR-Mega hits
SNPsToPlot = fread("../results/12_MR_Mega_Hits_including_SNP_from_table1.txt")
unique(SNPsToPlot[, setting])

#' ## Load data of Meta-GWAS with study specific values
eGFR_MALE = fread("../temp/12_MR_Mega/meta_results/GWASMA_eGFR_overall_ChrX_M_2021-12-16_12-10-52.gz")
UA_ALL = fread("../temp/12_MR_Mega/meta_results/GWASMA_uric_acid_overall_2021-06-23_11-41-05.gz")
UA_FEMALE = fread("../temp/12_MR_Mega/meta_results/GWASMA_uric_acid_overall_ChrX_F_2022-01-11_15-59-15.gz")
BUN_MALE = fread("../temp/12_MR_Mega/meta_results/GWASMA_BUN_ChrX_M_2022-01-27_14-46-00.gz")

#' ## Define function for Forest plots
#' It expects the SNP and the data object to use as input.
myForestPlot = function(snp, studydata) {
  #use data object specified
  setting = studydata
  studydata = get(studydata)
  
  #find row of SNP
  rowNR = which(studydata[, markerID] == snp)
  
  #find columns with beta, SE or samplesize values
  posBeta = colnames(studydata)[which(substring(colnames(studydata),1,5) == "beta.")]
  posSE = colnames(studydata)[which(substring(colnames(studydata),1,3) == "se.")]
  posN = colnames(studydata)[which(substring(colnames(studydata),1,2) == "n.")]
  posEAF = colnames(studydata)[which(substring(colnames(studydata),1,4) == "eaf.")]
  
  #get beta, SE values for single studies
  fp.mean = as.numeric(studydata[rowNR , get("posBeta"), with = F])
  myFilt = !is.na(fp.mean)
  fp.mean = fp.mean[!is.na(fp.mean)]
  fp.se = as.numeric(studydata[rowNR, get("posSE"), with = F])
  fp.se = fp.se[!is.na(fp.se)]
  fp.lower = fp.mean - fp.se
  fp.upper = fp.mean + fp.se
  fp.N = as.numeric(studydata[rowNR , get("posN"), with = F])
  fp.N = fp.N[!is.na(fp.N)]
  fp.EAF = as.numeric(studydata[rowNR , get("posEAF"), with = F])
  fp.EAF = fp.EAF[!is.na(fp.EAF)]
  fp.EAF = round(fp.EAF, digits = 2)
  
  #extract study names
  studyname = posBeta[myFilt]
  studyname = substring(studyname, 6, nchar(studyname))
  studyname = gsub("TRAILS_Pop", "TRAILS", studyname) #fix studyname for plot routine
  dummy = strsplit(studyname, split = "_")
  studyname = sapply(dummy, function(x) return(paste0(x[[1]], "_", x[[2]])))
  fp.labeltext = as.matrix(data.frame(study = studyname))
  fp.labeltext = paste0(fp.labeltext, " (N=", fp.N, ", EAF=", fp.EAF, ")")
  
  #extract ethnicity
  dummy = strsplit(studyname, split = "_")
  fp.eth = sapply(dummy, function(x) return(x[2]))
  
  #add result from MetaGWAS
  fp.mean = c(fp.mean, studydata[rowNR, betaFEM])
  fp.lower = c(fp.lower, studydata[rowNR, betaFEM] - studydata[rowNR, seFEM])
  fp.upper = c(fp.upper, studydata[rowNR, betaFEM] + studydata[rowNR, seFEM])
  fp.N = c(fp.N, studydata[rowNR, totalN])
  fp.labeltext = c(fp.labeltext, paste0("MetaGWAS", " (N=", studydata[rowNR, totalN], ", I2=", round(studydata[rowNR, I2], digits = 2), ")"))
  fp.eth = c(fp.eth, "Mixed")
  
  #build object for plotting
  plotData = data.frame(DV = paste0(setting, " - ", snp), study = fp.labeltext, beta = fp.mean, lower = fp.lower, upper = fp.upper, 
                        N = fp.N, ethnicity = fp.eth)
  
  #set order of rows in plot data object
  dummy = as.numeric(plotData$beta[-nrow(plotData)])
  plotOrder = order(dummy) 
  plotOrder = c(plotOrder, nrow(plotData))
  plotData$study = factor(plotData$study, levels = rev(plotData$study[plotOrder]))
  
  #set vectors for shape, size and colour of dots and lines in plot
  myShape = rep(15, times = nrow(plotData))
  myShape[grepl("MetaGWAS", plotData$study)] = 18  
  
  mySize = rep(1, times = nrow(plotData))
  mySize[plotData$N > 1000 & plotData$N < 10000] = 2
  mySize[plotData$N > 10000 & plotData$N < 100000] = 3
  mySize[plotData$N > 100000] = 4
    
  myCol = rep(1, times = nrow(plotData))
  myCol[plotData$ethnicity == "Mixed"] = 2
  myCol[plotData$ethnicity == "AA"] = 3
  myCol[plotData$ethnicity == "EAS"] = 4
  myCol[plotData$ethnicity == "HIS"] = 5
  
  myMin = min(plotData$lower)
  myMax = max(plotData$upper)
  
  #plot
  ggplot(data=plotData, aes(x=study, y=beta, ymin=lower, ymax=upper)) +
    geom_pointrange(colour = "white") + # Makes range for ggplot values based on the data and AES specified in first line
    geom_hline(yintercept = 0, lty=1, size = 0.5, color = "gray70") +  # add a dotted line at x=0 after flip
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.5, cex=0.5, colour = myCol) + # Makes whiskers on the range (more aesthetically pleasing)
    facet_wrap(~DV) + # Makes DV header (Can handle multiple DVs)
    coord_flip() + # flip coordinates (puts labels on y axis)
    geom_point(shape = myShape, size = mySize, colour = myCol) + # specifies the size and shape of the geompoint
    ggtitle("") + # Blank Title for the Graph
    xlab("") + # Label on the Y axis (flipped specification do to coord_flip)
    ylab("Beta") + # Label on the X axis (flipped specification do to coord_flip)
    scale_y_continuous(limits = c(myMin, myMax))  + # limits and tic marks on X axis (flipped specification do to coord_flip)
    theme(line = element_line(colour = "black", size = 2), # My personal theme for GGplots
          strip.background = element_rect(fill="gray95"), 
          legend.position ="none", 
          axis.line.x = element_line(colour = "black"), 
          axis.line.y = element_blank(), 
          panel.border= element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.ticks = element_blank(),
          axis.title.x = element_text(colour = "Black", margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(colour = "Black", margin = margin(t = 0, r = 10, b = 0, l = 0)),
          plot.title = element_text(colour = "Black", margin = margin(t = 0, r = 0, b = 10, l = 0)),
          axis.text=element_text(size=12, color = "Black"), 
          text=element_text(size=20), plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "mm"))
}


#' # Generate Forest Plot of table 1 index SNP with low pvalue anc_het
pdf(file = "../figures/SupplementalFigure_Forest_plot_table1_indexSNP.pdf", width = 8, height = 12)

foreach(r = 1) %do% {
  mySNP = SNPsToPlot[r, MarkerName]
  myData = SNPsToPlot[r, setting]
  
  myForestPlot(snp = mySNP, studydata = myData)
}

#' # Generate Forest Plots of MR-MEGA hits
pdf(file = "../figures/SupplementalFigure_Forest_plots_MR_MEGA_results.pdf", width = 8, height = 12)

foreach(r = c(2:nrow(SNPsToPlot))) %do% {
  mySNP = SNPsToPlot[r, MarkerName]
  myData = SNPsToPlot[r, setting]
  
  myForestPlot(snp = mySNP, studydata = myData)
}

dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
