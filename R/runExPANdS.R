runExPANdS<-function(SNV, CBS, maxScore=2.5, max_PM=6, precision=NA, plotF=2,snvF="out.expands",maxN=8000,region=NA){
if (!exists("SNV") || !exists("CBS")){
    print("Input-parameters SNV and/or CBS missing. Please add the paths to tabdelimited files containing the SNVs and copy numbers.");
    return();
}
if (is.na("maxScore")){
    maxScore=1;
}

if (is.na("max_PM")){
    max_PM=6;
}

if (is.na("plotF")){
    plotF=1;
}

nullResult=list("finalSPs"=NULL,"dm"=NULL,"densities"=NULL,"ploidy"=NULL);
dirF=getwd();
##SNVs
if (is.character(SNV) && file.exists(SNV)){
	print(paste("Running ExPANdS on: ",SNV))
	dm=read.table(SNV,sep="\t",header=TRUE);
	dm <- dm[ as.character(dm[,"chr"]) %in% as.character(seq(22)), ];
	dm=data.matrix(dm);
	print("Only SNVs with autosomal coordinates included.")
	snvF=fileparts(SNV);
	dirF=dirname(SNV); snvF=snvF$name;
}else if (is.matrix(SNV)){
	dm=SNV;
	if (!is.numeric(dm)) {
		print("SNV matrix has to be numeric. Likely cause: only mutations detected on autosomes accepted for ExPANdS model. Remove SNVs with allosomal and mitochondrial coordinates before you proceed.")
		return();
	}
}else{
	print("No SNVs provided. Aborting ExPANdS.");
	return();
}

##CBS
if (is.character(CBS) && file.exists(CBS)){
	copyNumber=as.matrix(read.table(CBS,sep="\t",header=TRUE))
}else if (is.matrix(CBS)){
	copyNumber=CBS;
}else{
	print("No copy number information provided. Aborting ExPANdS.")
	return();
}

dm=assignQuantityToMutation(dm,copyNumber,"CN_Estimate");
ii=which(is.na(dm[,"CN_Estimate"]));
if (length(ii)>0){
	print(paste(length(ii), " SNV(s) excluded due to unavailable copy number in that region."));
	dm=dm[-ii,];
}
ii=which(dm[,"CN_Estimate"]<1);
if (length(ii)>0){
        print(paste(length(ii), " SNV(s) excluded due to homozygous deletions within that region."));
        dm=dm[-ii,];
}
ii=which(dm[,"AF_Tumor"]*dm[,"CN_Estimate"]<0.1);
if (length(ii)>0){
        print(paste(length(ii), " SNV(s) excluded due to AF*CN below 0.1 (SNV can't be explained by an SP present in 10% or more of the sample)."));
        dm=dm[-ii,];
}


if (is.na(precision)){
   precision=0.1/log(nrow(dm)/7);
}
if(nrow(dm)<20){
   print("Not enough mutations provided. Minimum 20 SNVs required to attempt a run.")
   return(nullResult);
}

##Regions of interest
if (is.character(region) && file.exists(region)){
  roi=read.table(region,sep="\t",header=TRUE);
  roi <- roi[ as.character(roi[,"chr"]) %in% as.character(seq(22)), ];
  roi=data.matrix(roi);
}else if (is.matrix(region)){
  roi=region;
}

idx_R=c(1:size(dm,1));
if (size(dm,1)>maxN){
  if (is.na(region)){
    print(paste("Input contains more than ",maxN," SNVs and parameter <region> not set. Using default regional boundary (SureSelectExome_hg19)"));
  }

  if (!is.numeric(roi)) {
    print("Region matrix has to be numeric.");
    return(nullResult);
  }
  print(paste("Gathering SNVs within regions of interest ..."))
  ##Find SNVs within regions of interest
  idx_R=c();
  so=sort(roi[,'start'],index.return=T);
  roi=roi[so$ix,];
  for(i in 1:nrow(dm)){
    ii=which(roi[,'chr']==dm[i,'chr'])
    ix=which.min(abs(roi[ii,'start']-dm[i,'startpos']))
    ii=ii[c(max(0,(ix - 2)):min(length(ii),(ix + 2)))];
    if(i %% 100 ==0){
      print(paste(i," out of ",nrow(dm), " SNVs tested"))
    }
    for (j in ii){
      if (roi[j,'start']<=dm[i,'startpos'] && roi[j,'end']>=dm[i,'startpos']){
        idx_R=cbind(idx_R,i);
        break;
      }
    }
  }
}else{
  if (!is.na(region)){
    print(paste("Input contains less than ",maxN," SNVs. Parameter <region> ignored."));
  }
}

cfd=computeCellFrequencyDistributions(dm, max_PM, precision);
toUseIdx=which(apply(is.finite(cfd$densities),1,all) );
toUseIdx=intersect(toUseIdx,idx_R);
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, plotF-1,label=snvF);
if(is.null(SPs) || size(SPs,1)==0){
  print("No SPs found.")
  result=list("finalSPs"=NULL,"dm"=cfd$dm,"densities"=cfd$densities);
  return(result);
}
aM= assignMutations( cfd$dm, SPs,cfd$densities);
dm=aM$dm; finalSPs=aM$finalSPs;
if(is.null(dim(finalSPs))){
  finalSPs=finalSPs[finalSPs["score"]<=maxScore];
}else{
	finalSPs=finalSPs[finalSPs[,"score"]<=maxScore,];
}
if(!is.null(dim(finalSPs))){
	ia=order(finalSPs[,"Mean Weighted"]);
	finalSPs=matrix(finalSPs[ia,],nrow=nrow(finalSPs), ncol=ncol(finalSPs), dimnames=list(1:nrow(finalSPs),colnames(finalSPs)));
}
if(size(finalSPs,1)==0){
	print(paste("No SPs found below score:",maxScore))
	result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities);
	return(result);
}
aM = assignMutations( dm, finalSPs,cfd$densities);
dm=aM$dm; finalSPs=aM$finalSPs;

##ia=order(dm[,"%maxP"]);
##dm=dm[ia,]; cfd$densities=cfd$densities[ia,];

if (plotF>0){
    #jpgF=paste(dirF, .Platform$file.sep,snvF,".sps.jpg",sep="")
    #jpeg(filename = jpgF, width = 1080, height = 680); #X11();
    plotSPs(dm[toUseIdx,],snvF);
    ##dev.copy(jpeg,filename=paste(snvF,".sps.jpg"))
    #dev.off();
}
if (plotF>2){
    if(!require(rgl)){
	message("Plot supressed: Package rgl required for 3D plot of subpopulation clusters. Load this package before using this option.")
    }else{
    ##plot probability distributions
    cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
    kk=ceil(nrow(finalSPs)/2); par(mfcol=c(2,kk));
    for (i in 2:nrow(finalSPs)){
        idx=which(dm[toUseIdx,"SP"]==finalSPs[i,"Mean Weighted"]);
        open3d();
	        persp3d(as.numeric(1:length(idx)),as.numeric(cfd$freq),t(cfd$densities[idx,]),col=cols[mod(i,length(cols))+1],aspect=c(1, 1, 0.5), add=FALSE,xlab="Mutation", ylab="cell-frequency", zlab="Probability");
	        title3d(paste("SP_", round(finalSPs[i,"Mean Weighted"], digits=2),"\noo",sep=""));
        	play3d(spin3d(axis=c(0,0,1), rpm=10), duration=5)
	    }
    }
}

output=paste(dirF, .Platform$file.sep, snvF,".sps",sep="");
write.table(dm,file = output, quote = FALSE, sep = "\t", row.names=FALSE);
print(paste("Output saved under ",output));

##phylogeny
aQ=try(assignQuantityToSP(copyNumber, dm),silent=FALSE);
tr=NULL;
if(class(aQ)=="try-error" || is.null(ncol(aQ))){
	print("Error encountered while reconstructing phylogeny")
}else {
	output=paste(dirF, .Platform$file.sep, gsub("\\.","_",snvF), sep="");
	tr=buildPhylo(aQ,output);
}

if (plotF>1 && !is.null(tr)){
	plot(tr);
}

result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities,"ploidy"=aQ,"tree"=tr);
return(result);
}
