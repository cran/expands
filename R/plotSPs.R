plotSPs<-function(dm, sampleID=NA,cex=0.5){
  
  keep=which(!is.na(dm[,"SP"]),); dm=dm[keep,];
  ia=order(dm[,"%maxP"]);dm=dm[ia,];
  ia=order(dm[,"SP"],decreasing = TRUE);dm=dm[ia,];
  
  maxPloidy=max(dm[,"PM"],na.rm=TRUE);
  yticklab=c(1:maxPloidy,seq(0,1,by=0.1));
  at=c(sort(c(1:maxPloidy)*-0.1),seq(0,1,by=0.1));
  
  par(xpd=T, cex=cex, cex.axis=1/cex,cex.lab=1/cex, cex.main=1/cex,mar=par()$mar+c(0,0.5,0,4.2));
  #   plot.new(); 
  par(xpd=FALSE)
  
  plot(c(1:length(at)),at, col="white",pch=8,xlim=c(0,nrow(dm)),
       yaxt="n", bty="L", main=sampleID, xlab="Mutation", 
       ylab="Copy-number <---> Allele-frequency and SP size");
  axis(2, at, yticklab) 
  
  ia=order(dm[,"chr"]);dm=dm[ia,];
  ia=order(dm[,"SP"],decreasing = TRUE);dm=dm[ia,];
  
  legend1=.plotSPPerChr(dm,8,0,0);
  
  x=gray.colors(100)
  norm=1/length(x);
  for (k in 1:nrow(dm)){
    ci=max(1,ceil(dm[k,"%maxP"]/norm));
    matpoints(k,dm[k,"SP"],pch=15,col=x[ci]);
    if (k==1){
      legend1$text[length(legend1$text)+1]="SP";
      legend1$col[length(legend1$text)]=x[ci];
      par(xpd=TRUE)
      legend("topright",legend1$text,fill=legend1$col,inset=c(-0.1,-0.02),cex=0.75/cex,bty = "n")
    }
  }
  .plotSPPerChr(dm,17,1,0);
  .plotSPPerChr(dm,8,0,1);
  .plotSPPerChr(dm,8,1,1);
  lines(c(0,nrow(dm)),c(0,0),col="black");
  
  #ylim([-0.12*maxPloidy,1.1]);
}


.plotSPPerChr<-function(dm,lineType,lohFlag,cnFlag){
  x=rainbow(40);
  legend1=list("text"=c(),"col"=c());
  maxploidy=max(dm[,"PM"],na.rm=TRUE)+1;
  for (i in 1:22){
    idx=which(dm[,"chr"]==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
      if (!cnFlag){
        matpoints(idx,dm[idx,"AF_Tumor"],pch=lineType,col=x[i]);
      }else{
        if (any(!is.na(dm[idx,"PM"]))){
          matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
        }
      }
    }
    legend1$text[i]=paste("chr",i);
    legend1$col[i]=x[i];
  }
  return(legend1);
}
