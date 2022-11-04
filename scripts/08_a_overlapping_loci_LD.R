#' ---
#' title: "Calculate LD of SNPs in overlapping regions of eGFR and UA"
#' subtitle: "LD of SNPs in overlapping regions"
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
#' ## Load loci defintion
loci = fread("../results/01_Locus_Definitions.txt")
cols2keep = c("region", "region_start", "region_end", "ID_Meta", "rsID", "phenotype")
colsOut = setdiff(colnames(loci), cols2keep)
loci[, get("colsOut") := NULL]

# Add eGFR FEMALE hit for region 7 to this analysis
femaleHit = list(7.1, 100438892 , 101438892, "chr23:100479327:T:C", "rs149995096:100479327:C:T", "eGFR_FEMALE")
loci = rbindlist(list(loci, femaleHit))
setkey(loci, region)

#' ## Find overlapping regions
regions = loci[grepl("eGFR_", phenotype), region]
toDo = foreach(r = regions) %do% {
  myStart = loci[region == r, region_start]
  myEnd = loci[region == r, region_end]
  
  #four different types of overlaps
  a = myStart >= loci[, region_start] & myEnd <= loci[, region_end]
  b = myStart <= loci[, region_start] & myEnd > loci[, region_start] & myEnd <= loci[, region_end]
  c = myStart > loci[, region_start] & myStart <= loci[, region_end] & myEnd > loci[, region_end]
  d = myStart <= loci[, region_start] & myEnd >= loci[, region_end]
  
  e = as.numeric(a) + as.numeric(b) + as.numeric(c) + as.numeric(d)
  overlapping = loci[e > 0, region]
  overlapping = overlapping[!is.element(overlapping, r)]
  #remove overlap of regions 7 and 7.1
  if (r == 7 | r == 7.1) {
    overlapping = overlapping[overlapping > 8]
  }
  if (identical(overlapping, numeric(0))) { overlapping = NA }
  
  list(Region = r, Overlap = overlapping)
}
toDo = rbindlist(toDo)
toDo = toDo[!is.na(Overlap)]
toDo

#' ## Add UKBB ID to loci definition
step25 = fread("../results/05_b_step25_Matching_Table_UKBB.txt.gz")
matched = match(loci[, ID_Meta], step25[, markerID])
loci[, UKBB.ID := step25[matched, UKBB.ID]]


#' # Calculate LD of SNPs in overlapping regions
data.UKBB = "../results/05_a_ukb20272_imp_chrX_v3_s486631"

LD = foreach(l = c(1:nrow(toDo))) %do% {
  r1 = toDo[l, Region]
  r2 = toDo[l, Overlap]
  
  mySNP1 = loci[region == r1, UKBB.ID]
  mySNP2 = loci[region == r2, UKBB.ID]
  
  outFile = paste0("../temp/LD_calculation/08_a_", gsub(":", "_", mySNP1), "_", gsub(":", "_", mySNP2))
  plinkCall = paste0(path_plink2, " --bfile ", data.UKBB, " --ld ", mySNP1, " ", mySNP2, " --out ", outFile)
  system(plinkCall)
  
  #files have different length
  linenumber = system(paste0("grep -n r^2 ", outFile, ".log"), intern = T)
  linenumber = as.numeric(substring(linenumber, 1, 2)) - 1
  output = fread(paste0(outFile, ".log"), skip = linenumber, nrows = 1)
  dummy = list(eGFR.region = r1, eGFR.SNP = loci[mySNP1 == UKBB.ID, rsID], UA.region = r2, UA.SNP = loci[mySNP2 == UKBB.ID, rsID], var1 = output[, V1], r2 = output[, V3], var2 = output[, V4], Dprime = output[, V6])
  dummy
}
LD = rbindlist(LD)
LD
LD[, var1 := NULL]
LD[, var2 := NULL]
LD

#save results
write.table(LD, file = "../results/08_a_LD_overlapping_regions.txt", col.names = T, row.names = F, quote = F)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

