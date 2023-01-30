#' ---
#' title: "Heritability and Genetic Correlation"
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
#' We want to estimate the heritability for eGFR and UA on chromosome X and compare it to the autosomale heritability of previous CKDGen publications. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0 = Sys.time()

source("../SourceFile_aman.R")

setwd(paste0(projectpath,"scripts/"))

#' # Get ToDo ####
#' ***
ToDoList = data.table(results = list.files(path = "../temp/10_GCTA_Heritability/Step5_SlurmOutput/",pattern = "hsq"))
ToDoList[,path_fn := paste0("../temp/10_GCTA_Heritability/Step5_SlurmOutput/",results)]
dummy = ToDoList$results
dummy2 = unlist(strsplit(dummy,"_"))
ToDoList[,num := dummy2[seq(1,length(dummy2),4)]]
ToDoList[,num := as.numeric(num)]
ToDoList[,phenotype := dummy2[seq(2,length(dummy2),4)]]
ToDoList[,sex := dummy2[seq(3,length(dummy2),4)]]
ToDoList[,adjustment := dummy2[seq(4,length(dummy2),4)]]
ToDoList[,adjustment := gsub(".hsq","",adjustment)]
ToDoList[,c(1,3:6)]

#' # Load results ####
#' ***

dumTab = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  tab = readLines(myRow[,path_fn])
  tab
  
  tab2 = tab[5]
  tab2 = unlist(strsplit(tab2,"\t"))
  stopifnot(tab2[1] == "V(G)/Vp")

  tab3 = tab[6]
  tab3 = unlist(strsplit(tab3,"\t"))
  stopifnot(tab3[1] == "logL")
  
  tab4 = tab[11]
  tab4 = unlist(strsplit(tab4,"\t"))
  stopifnot(tab4[1] == "n")
  
  myRow[,h2 := tab2[2]]
  myRow[,h2_se := tab2[3]]
  myRow[,logL := tab3[2]]
  myRow[,n := tab4[2]]
  myRow
}
dumTab = rbindlist(dumTab)
dumTab[,c(1,3:10)]

dumTab[,h2 := as.numeric(h2)]
dumTab[,h2_se := as.numeric(h2_se)]

#' # Create Barplot ####
#' ***
plotdata = copy(dumTab)
plotdata[,lowerbound := h2-1.96*h2_se]
plotdata[,upperbound := h2+1.96*h2_se]
plotdata[lowerbound<0,lowerbound:=0]
plotdata[upperbound>1,upperbound:=1]

ggplot(plotdata, aes(x=phenotype, y=h2, fill=sex)) + 
  facet_wrap(~adjustment,scales = "free_x",ncol = 1) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lowerbound, ymax=upperbound), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(values=c("#B2182B","#2166AC"),labels = c("females","males"))+
  labs(x="phenotype",
       y = "heritability",
       fill = "sex") +
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

plotdata2 = copy(plotdata)
plotdata2 = plotdata2[adjustment == "Adj"]
plotdata2 = rbind(plotdata2,plotdata[c(1,3)])
plotdata2[c(5,6),sex := "combined"]

ggplot(plotdata2, aes(x=phenotype, y=h2, fill=sex)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=lowerbound, ymax=upperbound), width=.2,
                position=position_dodge(.9)) + 
  scale_fill_manual(values=c("#82B446","#B2182B","#2166AC"),labels = c("combined","females","males"))+
  labs(x="phenotype",
       y = "X-heritability",
       fill = "sex") +
  theme_bw(base_size = 10) + 
  theme(plot.title = element_text(hjust = 0, size=22,face="bold"),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"),
        strip.text.x = element_text(size = 12,face="bold"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#' # Save results ####
#' ***
myTab = copy(dumTab)
myTab[,path_fn := NULL]
myTab[,results := NULL]

write.table(myTab, file = "../results/10_GCTA_heritability.txt")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
