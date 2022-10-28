miamiPlot<-function(x,ymax=NULL,ymin=NULL, 
                               title="Miami Plot (GWAS, TWAS)",xlabel="chromosome",
                               ylabel=expression(paste("TWAS: ",log[10](p)," and GWAS: ",-log[10](p))),
                               hline1=NULL,hline2=NULL,sugline1=NULL,sugline2=NULL, highlight=F, 
                               returnObject = F,suppress19_21label =F, diffsize = FALSE , num_breaks_y=5, 
                               mypointsize=1,plotGenes=T,
                               showcolors =F, out_name,overall_max = 25, overall_min = -25,
                    useBasePosition=F) {
  
  # debug
  # x=plotData
  # ymax = ymaxpar1
  # ymin = -ymaxpar2
  # title = mytitle
  # xlabel = ""
  # ylabel=expression(paste("UA: ",log[10](p),"; eGFR: ",-log[10](p)))
  # hline1=-log10(5e-8)
  # hline2=log10(5e-8)
  # sugline1=-log10(1e-6)
  # sugline2=log10(1e-6)
  # highlight=T
  # diffsize = T
  # num_breaks_y=10
  # plotGenes=T
  # out_name="../results/Fig1_MiamiPlot_withGenes.pdf"
  # returnObject = T
  # overall_max = 40
  # overall_min = -40
  # suppress19_21label =F
  # showcolors =F
  # mypointsize=1

  if(suppress19_21label ==T) {x_labels = c(1:18, 20, 22:23)} else {x_labels= c(1:13, 15, 17, 19, 21, 23)}
  
  #load needed package
  library("ggplot2")
  library("scales")
  library(RColorBrewer)
  localenv <- environment()
  
  message("Preparing data")
  #check column names of object x
  check1<-match(c("CHR","BP","P"),colnames(x))
  if (sum(is.na(check1)) > 0) {
    stop("ERROR: Check the column names of your data object!")
  }
  
  #check if NAs are in x, if so -> remove them
  x<-x[!is.na(x$P),]
  
  #sort data set, right order is CHR:BP  
  x<-x[with(x, order(CHR, BP)), ]
  
  #recalculate given p values
  myY<- -log10(x$P)
  
  #calculate ymax
  if (is.null(ymax)) {
    ymax<- ceiling(max(myY))
    myText<-sprintf("ymax is set automatically to: %d",ymax)
    print(myText)
  }
  
  # set ymin to 0 if not provided
  if(is.null(ymin)) ymin = 0
  
  # calculate x coordinates for plot
  if(useBasePosition==T){
    myX = x$BP
  }else{
    myX<-seq(from=1,to=length(myY),by=1)
  }
  
  
  # choose color labels for chromosomes
  if(showcolors ==F) {
    myEven<-(x$CHR %% 2) == 0
    myColor<-vector(mode="character",length=dim(x)[1])
    myColor[myEven]<-2
    myColor[!myEven]<-1
  } else 	myColor = factor(x$CHR)
  
  
  # mylabel um spaeter noch die SNP - positionen rueckschlussfolgern zu koennen
  myLabel = x$SNP
  myFlag = x$flag
  myPhenotype = x$phenotype
  
  #build plot object for ggplot
  if(plotGenes==T){
    myGene = x$candidateGene
    myNovelty = x$Novelty
    mySexIA = x$SexInteraction
    myPD<-data.table(myX,myY,myColor, myLabel,myFlag,myPhenotype,myGene,myNovelty,mySexIA)
  }else{
    myPD<-data.table(myX,myY,myColor, myLabel,myFlag,myPhenotype)
  }
  
  myPD[myFlag=="bottom",myY:= myY*(-1)]
  
  if (highlight==T) {
    # I want colors per setting?
    # eGFR females - all - males: 
    myPD[,myPhenotype := gsub(".*_","",myPhenotype)]
    
    myPD[myY>hline1 & myPhenotype=="FEMALE" ,myColor:=3]
    myPD[myY>hline1 & myPhenotype=="ALL" ,myColor:=4]
    myPD[myY>hline1 & myPhenotype=="MALE" ,myColor:=5]
    
    myPD[myY<hline2 & myPhenotype=="FEMALE" ,myColor:="#D1E5F0"]
    myPD[myY<hline2 & myPhenotype=="ALL" ,myColor:=6]
    myPD[myY<hline2 & myPhenotype=="MALE" ,myColor:=7]
    
  }
  
  #recalculate position and values of x axis ticks and labels
  xAL<-unique(x$CHR)
  xAP<-sapply(xAL,function(chr) quantile(which(x$CHR==chr))[3])
  
  #recalculate position and values of y axis ticks and labels
  if (ymax < 50) { yAL<-seq(from=2,to=ymax,by=2) }else { yAL<-seq(from=5,to=ymax,by=5) }
  yAP<-yAL
  
  ## nun noch die x werte angezeigter art auf x achse ausduennen
  xAL = as.character(xAL)
  print(xAL)
  xAL[  xAL %in% as.character(x_labels) ==F] = ""
  print(xAL)
  
  
  myPD$pointsize = mypointsize
  if(diffsize == T & is.null(hline1)==F  & is.null(hline2)==F){ myPD[myY>hline1 | myY<hline2,pointsize:=2]}
  myPD$pointsize = factor(myPD$pointsize)
  
  # filter data to limit y-axis
  message("Filter data to limit y-axis")
  if(max(myPD$myY)>overall_max){
    sub1<-paste0("y axis limited to ", 
                 overall_max, " (max. original y-value: ", 
                 round(max(myPD$myY, na.rm = T),1), ")")
    myPD[myY>overall_max,myY:=overall_max]
    ymax<-overall_max
  } else sub1<-NULL
  
  if(min(myPD$myY)<overall_min){
    sub2<-paste0("y axis limited to ", 
                 overall_min, " (min. original y-value: ", 
                 round(min(myPD$myY, na.rm = T),1), ")")
    myPD[myY<overall_min,myY:=overall_min]
    ymin<-overall_min
  } else sub2<-NULL
  
  if(!is.null(sub1) & !is.null(sub2)){
    mySubtitle<-paste0(sub2, " and dummy",sub1)
    mySubtitle<-gsub("dummyy axis limited to","",mySubtitle)
    mySubtitle
  }
  if(!is.null(sub1) & is.null(sub2)) mySubtitle<-sub1
  if(is.null(sub1) & !is.null(sub2)) mySubtitle<-sub2
  if(is.null(sub1) & is.null(sub2)) mySubtitle<-""
  
  #print the manhattan plot 
  message("Start plot ...")
  myPlot <- ggplot() 
  # data=myPD, aes(x=myX, y=myY, colour=myColor, label = myLabel)
  myPlot <- myPlot + 
    geom_point(data=myPD, aes(x=myX, y=myY, colour=myColor, label = myLabel, size = pointsize)) + 
    scale_size_discrete(range=c(1,2.5)) + 
    guides(size="none")
  if(is.null(highlight) & showcolors!= F) {
    myPlot <- myPlot + scale_colour_manual(values= rep(rev(brewer.pal(5, "Dark2")), 23))
  }else { 
    myPlot <- myPlot + scale_colour_manual(values=c("#000000","#FDDBC7","#EF8A62","#B2182B","#67A9CF","#2166AC"),
                                           labels=c("non sig.","eGFR FEMALES","eGFR ALL","eGFR MALES","UA ALL","UA MALES"))
  }
  if(useBasePosition == T){
    myPlot <- myPlot + coord_cartesian(xlim=c(min(myX)-1,max(myX)+1),ylim=c(ymin-1,ymax+1), expand = FALSE)
  }else{
    myPlot <- myPlot + coord_cartesian(xlim=c(0,length(myX)+1),ylim=c(ymin-1,ymax+1))
  }
  
  myPlot <- myPlot + scale_x_continuous(name=xlabel,labels=xAL, breaks=xAP)
  myPlot <- myPlot + scale_y_continuous(name=ylabel,breaks=pretty_breaks(n=num_breaks_y)) 
  myPlot <- myPlot + ggtitle(title,subtitle = mySubtitle)
  myPlot <- myPlot + labs(color = "Legend")
  myPlot <- myPlot + theme(legend.position = c(0.8, 0.9)) 
  # myPlot <- myPlot + guides(color="none")
  myPlot <- myPlot + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
  myPlot <- myPlot + theme(axis.text.y = element_text(colour="black",size = 12),axis.line.y = element_line(colour="black"))
  myPlot <- myPlot + theme(panel.background = element_blank())
  myPlot <- myPlot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #myPlot <- myPlot + theme(axis.line = element_line(colour="black"))
  myPlot <- myPlot + geom_hline(yintercept = 0)
  if (!is.null(hline1)) {myPlot <- myPlot + geom_hline(yintercept = hline1,colour="#990000", linetype="dashed")}
  if (!is.null(hline2)) {myPlot <- myPlot + geom_hline(yintercept = hline2,colour="#990000", linetype="dashed")}
  if (!is.null(sugline1)) {myPlot <- myPlot + geom_hline(yintercept = sugline1,colour="#000099", linetype="dotted")}
  if (!is.null(sugline2)) {myPlot <- myPlot + geom_hline(yintercept = sugline2,colour="#000099", linetype="dotted")}
  
  # x axis line
  # myPlot <- myPlot + geom_segment(x = xAP[1], xend = xAP[1], y = 0, yend = 0 + 0.5, colour = "black")
  # myPlot <- myPlot + geom_text(x=xAP[1], y=-0.5, label=xAL[1], colour = "black")
  #dumTab<-data.table(label=xAL,num=1:23,breaks=xAP)
  
  #myPlot <- myPlot + geom_segment(data=dumTab, aes(x = breaks, xend = breaks, y = 0, yend = 0 + 0.5), colour = "black")
  #myPlot <- myPlot + geom_text(data=dumTab,aes(x=breaks, y=-0.75, label=label), colour = "black")
  
  # add genenames 
  myPlot2 <-myPlot + 
    geom_label_repel(data = subset(myPD, myFlag=="top" & myY>=7.3 & myNovelty==T),
                     aes(x=myX, y=myY, label = myGene),
                     ylim = c(ymax-5,ymax+5)) + 
    geom_label_repel(data = subset(myPD, myFlag=="top" & myY>=7.3 & myNovelty==F),
                     aes(x=myX, y=myY, label = myGene),color = "#666666",
                     ylim = c(ymax-15,ymax-10)) +
    geom_label_repel(data = subset(myPD, myFlag=="bottom" & myY<=-7.3 & myNovelty==T),
                     aes(x=myX, y=myY, label = myGene),
                     ylim = c(ymin-5,ymin+5)) + 
    geom_label_repel(data = subset(myPD, myFlag=="bottom" & myY<=-7.3 & myNovelty==F),
                     aes(x=myX, y=myY, label = myGene),color = "#666666",
                     ylim = c(ymin+10,ymin+5))
  
  
  myPlot2
  
  # check time for saving plot!
  timeplot0 = Sys.time()
  message("Create PDF")
  
  pdf_from_png(code2parseOrPlot = myPlot2, 
               pdf_filename = out_name,
               weite = 12,
               laenge = 8,
               einheiten = "in",
               resolution = 150)
  
  message("\nTOTAL TIME : " ,formateTimediff(Sys.time()-timeplot0))
  
  if(returnObject) return(myPlot2)
  if(returnObject ==F) print(myPlot2)
  
}
