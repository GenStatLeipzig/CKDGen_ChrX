#' ---
#' title: "Locus Definition"
#' subtitle: "Define loci for all future analysis"
#' author: "Katrin Horn, Janne Pott"
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
#' Definition of the genome wide significant loci of all phenotypes and settings. 
#' This important loci are the basis of the paper and all secondary analyses.
#' 
#' # Initialize ####
#' ***
#' Update server & set working directory
time0 = Sys.time()

source("../SourceFile_forostar.R")

setwd(projectpath_main)

#' # Analysis ####
#' ***
#' ## get genome-wide sig. SNPs ####
#' 
#' 1. Step: get genome-wide significant SNPs of all phenotypes and settings 
myFiles = list.files(path = path_data, pattern = "CKDGen_ChrX_sumStat") 

result.1 = foreach(f = myFiles) %do% {
  input_fn = paste0(path_data, f)
  dat = fread(input_fn)
  
  #filter invalid associations and for genome-wide significance
  dat = dat[(invalid_assoc == FALSE) & (P < 5*10^-8), ]
  dat
}
result.1 = rbindlist(result.1)
table(result.1[, phenotype])

#' ## Shrink SNP list by 500 kb range ###
#' 
#' 2. Step: shrink list to best SNPs in phenotype and setting (loop: take best SNP and remove SNPs in +/- 500 kb range until no SNPs left to remove from range of a better SNP)
subsets = unique(result.1[, phenotype])

# helper function to get minimal distance between significant SNPs --> do we have to collapse or not?
getSmallestDist = function(x) {
  y = c(x[2:length(x)], max(x)+1000000)
  return(min(y-x))
}

result.2 = foreach(s = subsets) %do% {
  # s = subsets[1]
  subdata = copy(result.1)
  subdata = subdata[phenotype == s, ]
  setkey(subdata, position)
  subdata[, keep := NA]
  subdata[, NR_SNPs := as.numeric(NA)]
  
  smallestDist = getSmallestDist(subdata[, position])
  while(smallestDist < 500000) {
    minP = min(subdata[is.na(keep), P])
    subdata[minP == P, keep := T]
    myPos = subdata[minP == P, position]
    stopifnot(length(myPos)==1)
    #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
    myFilt = (subdata[, position] < (myPos - 500000)) | 
             (subdata[, position] > (myPos + 500000)) | 
              subdata[, keep] 
    myFilt[is.na(myFilt)] = FALSE
    subdata = subdata[myFilt == TRUE, ]
    #NR_SNPs sums the filtered SNPs for a locus
    subdata[minP == P, NR_SNPs := sum(myFilt==F)]
    smallestDist = getSmallestDist(subdata[, position])
  }
  
  stopifnot(sum(is.na(subdata[,keep])) <= 1)
  subdata[is.na(keep), NR_SNPs := 0]
  subdata[is.na(keep), keep := TRUE]
  subdata
}
result.2 = rbindlist(result.2)

#add regions
result.2[, region_start := position - 500000]
result.2[, region_end := position + 500000]
result.2

#' ## Shrink SNP list by phenotype ###
#' 
#' 3. Step: again shrink list to best SNPs of phenotypes eGFR and uric acid (UA), combine regions of settings if needed (eGFR had a male-specific hit, hence eGFR_MALE is included in this step) and keep best SNP per new region

result.2 = result.2[phenotype %in% c("eGFR_ALL", "eGFR_MALE","UA_ALL")]
subsets = c("eGFR", "UA")

result.3 = foreach(s = subsets) %do% {
  # s = subsets[2]
  subdata = copy(result.2)
  subdata = subdata[grepl(s, phenotype), ]
  setkey(subdata, position)
  subdata[, keep := NA]
  
  #assign new region number to all lines of data.table
  subdata[, region := as.numeric(NA)]
  subdata[1, region := 1]
  foreach(l = c(2:nrow(subdata))) %do% {
    if (subdata[l, region_start] <= subdata[l-1, region_end]) {
      subdata[l, region := subdata[l-1, region]]
    } else { subdata[l, region := subdata[l-1, region] + 1]}
  }
  
  #keep best SNPs per region
  allRegions = unique(subdata[, region])
  
  subresult = foreach(r = allRegions) %do% {
    subdata.2 = copy(subdata)
    subdata.2 = subdata.2[region == r, ]
    minP = min(subdata.2[, P])
    subdata.2[P == minP, region_start := min(subdata.2[, region_start])]
    subdata.2[P == minP, region_end := max(subdata.2[, region_end])]
    subdata.2[P == minP, NR_SNPs := sum(subdata.2[,NR_SNPs]) + nrow(subdata.2) - 1]
    subdata.2[P == minP, ]
  }
  subresult = rbindlist(subresult)
  subresult
}
result.3 = rbindlist(result.3)
result.3[, region := c(1:nrow(result.3))]
result.3

cols2keep = c("region", "region_start", "region_end", "ID_Meta", "rsID", "phenotype", "chromosome", "position", "numberOfStudies", 
              "N", "I2", "EAF", "MAF", "infoScore", "effect_allele", "other_allele", "beta", "SE", "P", "logP", "NR_SNPs")
colsOut = setdiff(colnames(result.3), cols2keep)
result.3[, get("colsOut") := NULL]
setcolorder(result.3, cols2keep)

write.table(result.3, file = "../results/01_Locus_Definitions.txt", col.names = T, row.names = F, quote = F, sep = ",", dec = ".")

#' ## Shrink SNP list overall ###
#' 
#' 4. Step: again shrink list to best SNPs of either phenotype
#' 
#'  Same procedure as step 3. See how regions collapse over phenotypes.
subdata = copy(result.2)
setkey(subdata, position)

#assign new region number to all lines of data.table
subdata[, region := as.numeric(NA)]
subdata[1, region := 1]
done = foreach(l = c(2:nrow(subdata))) %do% {
  if (subdata[l, region_start] <= subdata[l-1, region_end]) {
    subdata[l, region := subdata[l-1, region]]
  } else { subdata[l, region := subdata[l-1, region] + 1]}
}

#keep best SNPs per region
allRegions = unique(subdata[, region])
  
result.4 = foreach(r = allRegions) %do% {
  subdata.2 = copy(subdata)
  subdata.2 = subdata.2[region == r, ]
  minP = min(subdata.2[, P])
  subdata.2[P == minP, region_start := min(subdata.2[, region_start])]
  subdata.2[P == minP, region_end := max(subdata.2[, region_end])]
  subdata.2[P == minP, NR_SNPs := sum(subdata.2[,NR_SNPs]) + nrow(subdata.2) - 1]
  subdata.2[P == minP, phenotypes_in_region := paste(unique(subdata.2[, phenotype]), collapse = " | ")]
  subdata.2[P == minP, ]
}
result.4 = rbindlist(result.4)

cols2keep = c("region", "region_start", "region_end", "ID_Meta", "rsID", "phenotype", "chromosome", "position", "numberOfStudies", 
              "N", "I2", "EAF", "MAF", "infoScore", "effect_allele", "other_allele", "beta", "SE", "P", "logP", "NR_SNPs", "phenotypes_in_region")
colsOut = setdiff(colnames(result.4), cols2keep)
result.4[, get("colsOut") := NULL]
setcolorder(result.4, cols2keep)
result.4

write.table(result.4, file = "../results/01_Loci_collapsed_overall.txt", col.names = T, row.names = F, quote = F, sep = ",", dec = ".")


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

