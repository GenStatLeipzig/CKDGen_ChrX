#' @title Function to perform a co-localization test
#' @description This function takes two locus-wide summary statistic sets and matches them by ID as provided. Then, it tests for co-localization using beta, se, sample size, MAF, and ID. If the plotting option is TRUE, then a co-localization plot is generated and the coloc result is included in the sub title. Both data sets need to have the same columns!
#' @param tab1 data set for trait 1
#' @param tab2 data set for trait 2
#' @param trait1 phenotype name for trait 1, Default: 'trait1'
#' @param trait2 phenotype name for trait 1, Default: 'trait2'
#' @param locus cytoband of tested region; will be used to label the X axis of the coloc plot, Default: 'cytoband'
#' @param locus_name gene name of tested region (eg if some specific eQTL-coloc is done); will be used in title of the coloc plot, Default: 'gene'
#' @param plotting should a plot be created?, Default: F
#' @param col_SNPID column name for matching, Default: 'markername'
#' @param col_pos column name for base position, Default: 'pos'
#' @param col_beta column name for beta estimate, Default: 'beta'
#' @param col_se column name for standard error, Default: 'se'
#' @param col_P column name for p-value (used for plotting), Default: 'pval'
#' @param col_N column name for sample size, Default: 'n_samples'
#' @param col_MAF column name for minor allele frequency, Default: 'maf'
#' @return first entry of coloc.abf result: vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
#' @details see [coloc](https://cran.r-project.org/web/packages/coloc/index.html) package for more information:
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{setcolorder}}
#'  \code{\link[coloc]{coloc.abf}}
#' @rdname colocFunction_jp
#' @export
#' @importFrom data.table setcolorder
#' @importFrom coloc coloc.abf
colocFunction_jp<-function(tab1,tab2,trait1 = "trait1",trait2 = "trait2",locus = "cytoband",locus_name = "gene",
                           plotting=F, col_SNPID = "markername", col_pos = "pos", col_beta = "beta", col_se = "se",
                           col_P = "pval", col_N = "n_samples", col_MAF = "maf",col_effect ="effect_allele", 
                           col_other = "other_allele" ){
  # tab1 = copy(data_GWAS1)
  # tab2 = copy(data_eQTL1)
  # trait1 = myPheno
  # trait2 = myTissue2
  # locus = moreInfo$cytoband
  # locus_name = moreInfo$gene
  # plotting=F
  # col_SNPID = "chrPos_b37"
  # col_pos = "pos_b37"
  # col_beta = "beta"
  # col_se = "se"
  # col_P = "pval"
  # col_N = "n_samples"
  # col_MAF = "maf"
  # col_effect = "effect_allele"
  # col_other = "other_allele"
  
  # Step 1: check if all provided column names are there; maf allowed to be missing in external data
  myNames = c(col_SNPID, col_pos, col_beta, col_se, col_P, col_N,col_effect,col_other,col_MAF)
  stopifnot(sum(myNames %in% colnames(tab1))==9)
  stopifnot(sum(myNames %in% colnames(tab2))>=8)
  
  if(sum(myNames %in% colnames(tab2))==8){
    stopifnot(sum(myNames[1:8] %in% colnames(tab2))==8)
    tab2[,maf:=0.25]
    updateMAF = T
  }else{
    updateMAF = F
  }
  
  # Step 2: reduce date sets to relevant columns and order them
  colsOut<-setdiff(colnames(tab1),myNames)
  tab1[,get("colsOut"):=NULL]
  data.table::setcolorder(tab1,myNames)
  
  colsOut<-setdiff(colnames(tab2),myNames)
  tab2[,get("colsOut"):=NULL]
  data.table::setcolorder(tab2,myNames)
  
  # Step 3: rename columns
  myNames2<-c("markername","pos","beta","se","pval","n_samples","EA","OA","maf")
  names(tab1)<-myNames2
  names(tab2)<-myNames2
  
  # Step 4: check SNP overlap and match data set 2 to data set 1
  table(is.element(tab1$markername,tab2$markername))
  table(is.element(tab2$markername,tab1$markername))
  
  dumTab1<-tab1[is.element(tab1$markername,tab2$markername),]
  dumTab2<-tab2[is.element(tab2$markername,tab1$markername),]
  
  matched<-match(dumTab1$markername,dumTab2$markername)
  dumTab2<-dumTab2[matched,]
  stopifnot(dumTab1$markername==dumTab2$markername)
  
  # Step 5: check if the same allele were used
  filt1 = dumTab1$EA == dumTab2$EA & dumTab1$OA == dumTab2$OA
  filt2 = dumTab1$EA == dumTab2$OA & dumTab1$OA == dumTab2$EA
  table(filt1,filt2)
  filt = filt1 | filt2
  dumTab1<-dumTab1[filt,]
  dumTab2<-dumTab2[filt,]
  stopifnot(dumTab1$markername==dumTab2$markername)
  
  
  # Step 6: check for NA in beta estimates (should have been filtered before!)
  filt<-is.na(dumTab1$beta) | is.na(dumTab2$beta)
  dumTab1<-dumTab1[!filt,]
  dumTab2<-dumTab2[!filt,]
  
  # Step 6: perform coloc
  if(updateMAF==T){
    my_res<- coloc::coloc.abf(dataset1=list(beta=dumTab1$beta,
                                            varbeta=(dumTab1$se)^2,
                                            N=dumTab1$n_samples,
                                            snp=dumTab1$markername,
                                            type="quant"),
                              dataset2=list(beta=dumTab2$beta,
                                            varbeta=(dumTab2$se)^2,
                                            N=dumTab2$n_samples,
                                            snp=dumTab2$markername,
                                            type="quant"),
                              MAF=dumTab1$maf)
  }else{
    my_res<- coloc::coloc.abf(dataset1=list(beta=dumTab1$beta,
                                            varbeta=(dumTab1$se)^2,
                                            N=dumTab1$n_samples,
                                            snp=dumTab1$markername,
                                            MAF=dumTab1$maf,
                                            type="quant"),
                              dataset2=list(beta=dumTab2$beta,
                                            varbeta=(dumTab2$se)^2,
                                            N=dumTab2$n_samples,
                                            snp=dumTab2$markername,
                                            MAF=dumTab2$maf,
                                            type="quant"))
  }

  
  my_res2<-my_res[[1]]
  
  # Step 7: plot if plotting == T
  if(plotting==T){
    myXlab<-locus
    
    dumTab1[,plot(pos, -log10(pval),col = rgb(0,0,1,0.3),pch=19,
                  xlab = myXlab,
                  ylab = "")]
    axis(side = 2, col = "blue", col.ticks = "blue",col.axis="blue")
    mtext(side = 2, line = 2.1,
          bquote(.(trait1) ~ -log[10](italic(p-value))), col = "blue", cex = 0.7)
    
    par(new = T)
    dumTab2[,plot(pos, -log10(pval),axes=F, xlab=NA, ylab=NA,
                  col = rgb(1,0,0,0.2), pch = 17)]
    axis(side = 4, col = "red", col.ticks = "red",col.axis="red")
    mtext(side = 4, line = 2.1,
          bquote(.(trait2) ~ -log[10](italic(p-value))), col = "red", cex = 0.7)
    
    #legendtext = paste0("APOB gene\nmore GWAS-risk,\nmore expression")
    #legend("topright",legend=legendtext,  text.col=c("black"), bty="n")
    
    title(paste0(trait1," vs. ",trait2," at ", locus_name), cex.main = 1)
    
    mtext(paste0("H0 | H1 | ... | H4:  ", paste(formatC(round(my_res2[2:6],3), digits=3, format="f" ), collapse = " | ")),cex = 0.7)
  }
  
  # Step 8: return coloc result
  return(my_res2)
}
