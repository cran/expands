runExPANdS<-function(SNV, CBS, maxS=0.7, max_PM=6, min_CF=0.1, p=NA, ploidy=2, 
                     nc=1, plotF=2, snvF=NULL, maxN=8000, region=NA, verbose=T){
  MINSNVS=15
  
  if (!exists("SNV") || !exists("CBS")){
    print("Input-parameters SNV and/or CBS missing. Please add the paths to tabdelimited files containing the SNVs and copy numbers.");
    return();
  }

  dirF=getwd();
  ##Read point mutation and copy number input
  tmp=.readSNVandCBS(SNV,CBS,max_PM=max_PM, min_CF=min_CF, snvF=snvF, verbose=verbose);
  dm=tmp$dm; copyNumber=tmp$copyNumber; snvF=tmp$snvF;
  
  ii=which(dm[,"CN_Estimate"]<0.5);
  if (length(ii)>0){
    print(paste(length(ii), " SNV(s) excluded due to homozygous deletions within that region."));
    dm=dm[-ii,];
  }
  
  if(is.null(dm) || is.null(dim(dm)) || nrow(dm)<MINSNVS){
    print("Not enough mutations provided. Minimum",MINSNVS,"SNVs required to attempt a run.")
    return(list("finalSPs"=NULL,"dm"=NULL,"densities"=NULL,"sp_cbs"=NULL));
  }
  
  ##Regions of interest
  idx_R=.indexOfSNVsWithin(region, dm, maxN, verbose=verbose)
  if (is.na(p)){
    p=0.1/log(length(idx_R)/7);
  }

  #############################
  ###Input parsing ends here###
  cfd=computeCellFrequencyDistributions(dm = dm, max_PM = max_PM, p = p, 
                                        ploidy=ploidy, min_CF=min_CF, nc=nc, v=verbose);
  toUseIdx=which(apply(is.finite(cfd$densities),1,all) );
  toUseIdx=intersect(toUseIdx,idx_R);
  SPs=clusterCellFrequencies(densities = cfd$densities[toUseIdx,], p = p, min_CF=min_CF,verbose=verbose);
  if(is.null(SPs) || size(SPs,1)==0 || all(is.na(SPs)) || sum(SPs[,"score"]<=maxS)==0 ){
    print("No SPs found.")
    result=list("finalSPs"=NULL,"dm"=cfd$dm,"densities"=cfd$densities);
    return(result);
  }
  finalSPs=SPs[SPs[,"score"]<=maxS,,drop=F];
  ia=order(finalSPs[,"Mean Weighted"]);
  finalSPs=matrix(finalSPs[ia,,drop=F],nrow=nrow(finalSPs), ncol=ncol(finalSPs), dimnames=list(1:nrow(finalSPs),colnames(finalSPs)));

  if(size(finalSPs,1)==0){
    print(paste("No SPs found below score:",maxS))
    result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities);
    return(result);
  }
  ###########################
  ###Assigning SNVs to SPs###
  aM = assignMutations( dm = cfd$dm, finalSPs = finalSPs, max_PM=max_PM, ploidy=ploidy,verbose=verbose);
  dm=aM$dm; SPs=aM$finalSPs;
  
#   ########################################################################
#   ##Remove SPs with no single SNV/CNV cooccurence and reassign mutations##
#   iF=which(SPs[,'snv_cnv_Co']==0);
#   while(length(iF)>0 && nrow(SPs)>2){ ##Remove only if enough SPs are left
#     iF=iF[which.max(SPs[iF,'score'])]
#     .notifyUser(paste("Excluding subpopulation",SPs[iF,'Mean Weighted'],"with 0 SNV/CNV cooccurence incidence."),verbose=verbose)
#     SPs=SPs[-iF,,drop=F]
#     aM = assignMutations( dm = cfd$dm, finalSPs = SPs, max_PM=max_PM, ploidy=ploidy,verbose=verbose);
#     dm=aM$dm; SPs=aM$finalSPs;
#     iF=which(SPs[,'snv_cnv_Co']==0);
#   }
  
  ###########################
  ##Choose between doublets##
  continueprune=TRUE;
  while (!is.null(dim(SPs)) && nrow(SPs)>1 && continueprune){
    x=sort(SPs[,'Mean Weighted'],index.return=TRUE);   SPs=SPs[x$ix,];
    toRm=c();
    for  (sp in 1:nrow(SPs)){
      ii=sp; ##sp_i * x != sp_j for all x in 2:6 and all SP pairs (i,j)
      x=2; ##Check for doublets only
      maxDev=SPs[sp,'precision']* 2/3 ;
      i_=which(abs(SPs[sp,'Mean Weighted']*x-SPs[,'Mean Weighted'])<maxDev);
      ii=union(ii,i_); 
      if( length(ii)>1 ){
        ii=ii[which(SPs[ii,'nMutations']<0.6*max(SPs[ii,'nMutations'],na.rm=T))]; ##Keep only SP of max kurtosis
        toRm=c(toRm,ii )
      }
    }
    
    toRm=unique(toRm)
    if (!is.null(toRm) && length(toRm)>0){
      iReassign=which(dm[,'SP'] %in% SPs[toRm,'Mean Weighted'] | dm[,'SP_cnv'] %in% SPs[toRm,'Mean Weighted'] );
      print(paste('Reassigning SNVs after pruning',length(toRm),'doublet subpopulation(s).'))
      .notifyUser("Pruned subpopulation(s):",verbose = verbose)
      .notifyUser(SPs[toRm,'Mean Weighted'],verbose = verbose)
      SPs=SPs[-toRm,, drop=FALSE];
      aM = assignMutations(dm = dm[iReassign,,drop=FALSE], finalSPs = SPs,max_PM=max_PM, ploidy=ploidy,verbose=verbose);
      dm[iReassign,]=aM$dm; 
      #Recount mutations after reassignment (absence of pruned SP_cnv may cause reassignment of mutations among remaining SPs):
      SPs[,"nMutations"]=0;
      tmp=count(dm[,"SP"]); ia=match(tmp$x,SPs[,"Mean Weighted"]);
      SPs[ia,"nMutations"]=tmp$freq
      finalSPs=SPs;
    }else{
      continueprune=FALSE;
      finalSPs=SPs;
    }
  }
  
  #if (plotF>2){
  #    if(!require(rgl)){
  #	message("Plot supressed: Package rgl required for 3D plot of subpopulation clusters. Load this package before using this option.")
  #    }else{
  #    ##plot probability distributions
  #    cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
  #    kk=ceil(nrow(finalSPs)/2); par(mfcol=c(2,kk));
  #    for (i in 2:nrow(finalSPs)){
  #        idx=which(dm[toUseIdx,"SP"]==finalSPs[i,"Mean Weighted"]);
  #        open3d();
  #	        persp3d(as.numeric(1:length(idx)),as.numeric(cfd$freq),t(cfd$densities[idx,]),col=cols[mod(i,length(cols))+1],aspect=c(1, 1, 0.5), add=FALSE,xlab="Mutation", ylab="cell-frequency", zlab="Probability");
  #	        title3d(paste("SP_", round(finalSPs[i,"Mean Weighted"], digits=2),"\noo",sep=""));
  #        	play3d(spin3d(axis=c(0,0,1), rpm=10), duration=5)
  #	    }
  #    }
  #}
  
  ##phylogeny
  finalSPs=.addColumn(finalSPs,'Ancestor',NA);
  finalSPs=.addColumn(finalSPs,'ClosestDescendant',NA);
  aQ=try(assignQuantityToSP(cbs = copyNumber, dm = dm,v=verbose),silent=FALSE);
  tr=NULL;
  if(class(aQ)=="try-error"){
    print("Error encountered while assigning subpopulation specific copy number to genomic segments. Phylogeny will not be inferred.")
  }else {
    .writeExpandsOutput(X=aQ, dirF,snvF,suffix=".sps.cbs", message="Subpopulation specific copy numbers")
    output=paste(dirF, .Platform$file.sep, snvF, sep="");    # output=paste(dirF, .Platform$file.sep, gsub("\\.","_",snvF), sep="");
    tr=try(buildPhylo(sp_cbs = aQ,outF = output,dm=dm,verbose=verbose),silent=FALSE);
    if(class(tr)!="try-error" ){
      if(class(tr$dm)!="try-error" && !is.na(tr$dm)){
        dm=tr$dm;
        ##Assign ancestor and closest descendant
        phySPs=colnames(tr$spRelations)
        for (sp in rownames(tr$spRelations)){
          desc=as.numeric(gsub('SP_','',c(sp,phySPs[tr$spRelations[sp,]==1])))
          if(length(desc)>1){
            iAnc=which.min(abs(finalSPs[,'Mean Weighted']-desc[1]))
            desc=desc[which.min(abs(desc[2:length(desc)]-desc[1]))+1]
            iDesc=which.min(abs(finalSPs[,'Mean Weighted']-desc))
            finalSPs[iDesc,'Ancestor']=finalSPs[iAnc,'Mean Weighted']
            finalSPs[iAnc,'ClosestDescendant']=finalSPs[iDesc,'Mean Weighted']
          }
        }
      }else{
        print("Error encountered while reconstructing phylogeny")
      }
      tr=tr$tree;
    }
  }
  
  ####################
  ###Write output#####
  .writeExpandsOutput(X=dm, dirF,snvF,suffix=".sps", message="Subpopulation specific point mutations")
  .writeExpandsOutput(X=finalSPs, dirF,snvF,suffix=".spstats", message="Summary file of detected subpopulations")
  
  if (plotF>0){
    tmp=try(plotSPs(dm = dm[toUseIdx,],sampleID = snvF),silent=FALSE);
  }
  if (plotF>1 && !is.null(tr) && class(tr)!="try-error"){
    try(plot(tr),silent=FALSE);
  }
  
  result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities,"sp_cbs"=aQ,"tree"=tr);
  return(result);
}


.indexOfSNVsWithin<-function(region, dm, maxN, verbose){
  if (is.character(region) && file.exists(region)){
    roi=read.table(region,sep="\t",header=TRUE,check.names = F,stringsAsFactors = F);
    roi <- roi[ as.character(roi[,"chr"]) %in% as.character(seq(100)), ,drop=F];
    if(nrow(roi)==0){
      warning(paste(region,"does not contain genomic regions on known chromosomes. Regions will be ignored. Correct format does not include \'chr\'-prefix."))
    }
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
      return( list("finalSPs"=NULL,"dm"=NULL,"densities"=NULL,"sp_cbs"=NULL) );
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
        .notifyUser(paste(i," out of ",nrow(dm), " SNVs tested"),verbose=verbose)
      }
      for (j in ii){
        if (roi[j,'start']<=dm[i,'startpos'] && roi[j,'end']>=dm[i,'startpos']){
          idx_R=cbind(idx_R,i);
          break;
        }
      }
    }
    print(paste('Found ',length(idx_R), ' SNVs within regions of interest',sep=""))
    idx_R=sample(idx_R,min(maxN,length(idx_R)),replace=FALSE);
    print(paste("Keeping ",length(idx_R), ' of these SNVs (randomly selected).',sep=""))
  }else{
    if (!is.na(region)){
      print(paste("Input contains less than ",maxN," SNVs. Parameter <region> ignored."));
    }
  }
  return(idx_R)
}
