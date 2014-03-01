buildPhylo<-function(ploidy,outF,treeAlgorithm="bionjs"){
  #library(ape)
  ii=grep("SP",colnames(ploidy));
  cnv=ploidy[,ii];
  if(is.null(ncol(cnv))){
    print("Less than two SPs coexist in this tumor. Aborting phylogeny reconstruction");
    return(NULL);
  }
  
  print(paste("Building phylogeny using ",treeAlgorithm," algorithm",sep=""))
  print("Pairwise SP distances calculated as: % segments with identical copy number");
  ##add ancestor cell - ploidy:=consensus at all positions - if not already added
  if (!any(colnames(cnv)=="Consensus_SP",na.rm=T)){
    cnv=cbind(cnv,matrix(matrix(NaN,nrow(cnv),1),nrow=nrow(cnv), ncol=1, dimnames=list(rownames(cnv),"Consensus_SP")));
    cnv[,"Consensus_SP"]=round(colMeans(t(cnv),na.rm=T));
  }
 
  toRm=c(); 
  ##distance matrix from pairwise alignments
  cols=gsub(" ","",colnames(cnv));
  D=matrix(matrix(1,ncol(cnv),ncol(cnv)), nrow=ncol(cnv), ncol=ncol(cnv), dimnames=list(cols,cols));
  for (i in 1:ncol(cnv)){
    for (j in i:ncol(cnv)){
      ii=which(!is.na(cnv[,i]) & !is.na(cnv[,j]));
      if (length(ii)==0){
        D[i,j]<-D[j,i]<-NA;
        next;
      }
      x=cnv[ii,i];      y=cnv[ii,j];
      dd=length(which(x!=y))/length(ii);
      if (i!=j){
        dd=dd+0.3;
      }
      D[i,j]<-D[j,i]<-dd
    }
    if (any(is.na(D[i,1:i]))){
      toRm=cbind(toRm,i);#remove NAs
    }
  }
  #remove NAs
  if (nrow(D)-length(toRm)<2){
    print("No two SPs found between which distance could be calculated. Aborting phylogeny reconstruction");
    return(NULL);
  }
  if (length(toRm)>0){
    print(paste("Insufficient copy number segments for ",rownames(D)[toRm],". SP excluded from phylogeny",sep=""))
    D=D[-toRm,];
    D=D[,-toRm];
  }
  
  D=D*100;
  write.table(D,paste(outF,".dist",sep=""),quote=F,sep="\t");
  print(paste("distance-matrix saved under ",outF,".dist",sep=""));
  tr <- njs(D);
  if (treeAlgorithm=="bionjs"){
    tr <- bionjs(D);
  }
  tr$root.edge <- 0; ## adds root dummy
  write.tree(tr, file = paste(outF,".tree",sep=""));
  print(paste("tree saved under ",outF,".tree",sep=""));
  
  return(tr);
}
