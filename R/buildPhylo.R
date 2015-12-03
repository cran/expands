buildPhylo<-function(ploidy,outF,treeAlgorithm="bionjs",dm=NA){
  #library(ape)
  out=list("tree"=NULL,"dm"=dm);
  ii=grep("SP",colnames(ploidy));
  cnv=ploidy[,ii];
  if(is.null(ncol(cnv))){
    print("Less than two SPs coexist in this tumor. Aborting phylogeny reconstruction");
    return(out);
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
  if (length(toRm)>0){
    print(paste("Insufficient copy number segments for ",rownames(D)[toRm],". SP excluded from phylogeny",sep=""))
    D=D[-toRm,];
    D=D[,-toRm];
  }
  
  if(is.null(nrow(D)) || nrow(D)-length(grep("Consensus",rownames(D)))<2){
    print("No two SPs found between which distance could be calculated. Aborting phylogeny reconstruction");
    return(out);
  }
  
  D=D*100;
  #D=D[-1*which(rownames(D)=="Consensus_SP"),]
  #D=D[,-1*which(colnames(D)=="Consensus_SP")]
  write.table(D,paste(outF,".dist",sep=""),quote=F,sep="\t");
  print(paste("distance-matrix saved under ",outF,".dist",sep=""));
  tr =c();
  if (treeAlgorithm=="bionjs"){
    tr <- bionjs(D);
  }else{
    tr <- njs(D);
  }
  tr$root.edge <- 0; ## adds root dummy
  
  outF=paste(outF,".tree",sep="");
  write.tree(tr, file = outF);
  print(paste("tree saved under ",outF,sep=""));
  
  out$tree=tr;
  if(!is.na(dm)){
    dm1=try(.assignSNVsToMultipleSPs(dm,outF),silent=FALSE)
    if(class(dm1)!="try-error"){
      dm=dm1;
    }
    out$dm=dm;
  }
  return(out);
}

.assignSNVsToMultipleSPs <-function(dm,outF){
  if (!requireNamespace("phylobase", quietly = TRUE)) {
    print("Package \'phylobase\' needed for assigning SNVs to Multiple SPs. Please install it.")
    return(dm);
  }
  tr=phylobase::readNewick(outF,check.names=F);
  print("Assigning SNVs to SPs...")
  SPs = sort(unique(c(dm[, "SP_cnv"],dm[,"SP"])))
  spSizes = unique(round(SPs * 1000)/1000)
  spNames= paste("SP_", as.character(spSizes), sep = "");
  x = colnames(dm)
  x[(length(x) + 1):(length(x) + length(SPs))] =spNames
  dm = matrix(cbind(dm, matrix(0, nrow(dm), length(SPs))), nrow = nrow(dm), ncol = length(x), dimnames = list(1:nrow(dm), x))
  
  for (k in 1:nrow(dm)) {
    if(is.na(dm[k,"SP"])){
      next;
    }
    if (mod(k, 100) == 0) {
      print(paste("Assigning SPs for SNV", k, 
                  "out of ", nrow(dm), "..."))
    }
    thisSP = paste("SP ", as.character(round(dm[k,"SP"] * 1000)/1000),sep="");
    dm[k,gsub(" ","_",thisSP)]=1; #dm[k,"PM_B"]; binary assignment for now
    dm=.propagateSNVToMultipleSPs(thisSP,dm,k,tr,spSizes)
  }
  
  for (i in 1:(length(spNames)-1)){
    iPhylo=which(sum(t(dm[,spNames]!=0))>length(spNames)-i)
    print(paste(length(iPhylo), " SNVs assigned to >",length(spNames)-i," SPs"))
  }
  return(dm)
}

.propagateSNVToMultipleSPs <-function(thisSP,dm,k,tr, spSizes){
  xx=phylobase::getNode(tr,type="tip");
  ii_This=match(thisSP,names(xx)); ##Node representing SP which harbors this SNV
  sibl=phylobase::siblings(tr,xx[ii_This]); ##Siblings of SP with this SNV
  i_toRM=grep("Consensus",names(sibl));
  if(!isempty(i_toRM)){
    sibl=sibl[-i_toRM];
  }
  ij=match(names(sibl),paste("SP", as.character(spSizes), sep = " "));  ij=ij[!is.na(ij)];
  if (isempty(ij)){ ##No tree tips among Siblings
    return(dm);
  }
  
  for (s in 1:length(sibl)){
    if(is.na(names(sibl[s]))){
      next;
    }
    otherSP=gsub(" ","_",names(sibl[s]));
    if (dm[k,otherSP]==0 && phylobase::nodeDepth(tr,sibl[s])> 2*phylobase::nodeDepth(tr,xx[ii_This])){ ##This SP is likely ancestor of SP assigned as sibling. TODO: find tree reconstruction algorithm that can assign "living populations" as common ancestors
      dm[k,otherSP]=1; #dm[k,"PM_B"];  ##Other SP inherits mutation of this SP. 
      ##TODO: what if the ploidy of the otherSP in this region is < B-allele ploidy of this SP! need to check this and set to minimum. Temporary solution --> binary assignment (above)
      dm=.propagateSNVToMultipleSPs(gsub("_"," ",otherSP),dm,k,tr,spSizes)
    }
  }
  return(dm);
}

