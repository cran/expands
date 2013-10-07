runExPANdS<-function(SNV, CBS, maxScore=2.5, max_PM=6, precision=NA, plotF=1,snvF="out.expands"){
#RUNEXPANDS Computes coexistent subpopulations in a tumor from allele-frequency
#and copy number of mutated loci.
#     finalSPs = runExPANdS(SNV, CBS) assigns each SNV within SNV-input-file to a subpopulation. 
# Input-parameters SNV and CBS hold the paths to tabdelimited files containing the SNVs and the copy numbers respectively. 
# The files typically consists of column names in the first row, row names in the first column, and numeric data starting in the (2,2) position. 
# Columns in SNV must include:
#     > chr - the chrmosome on which each SNV is located
#     > startpos - the position of each SNV
#     > AF_Tumor - the allele-frequency of each SNV
#     > PN_B - binary variable: 1 if the SNV is a germline variant, 0 if
#     somatic.
# CBS is typically the output of a circular-binary-segmentation algorithm. Columns in CBS must include:
#     > chr - chromosome
#     > startpos - the first position of a copy number segment
#     > endpos - the last position of a copy number segment
#     > CN_Estimate - the copy number estimated for each segment
# An output file is saved under the SNV-input-file address with extension 'sps', in which the content of the input file 
# is extended by four additional columns:
#     > SP - the subpopulation to which each mutation has been assigned as
#     fraction of cells in the tumor bulk
#     > %maxP - the confidence with which each mutation has been assigned to
#     a subpopulation - between 0 (lowest) and 100 (highest).
#     > PM - the number of copies of each mutated locus in the assigned SP
#     > PM_B - the number of copies of each SNV in the assigned SP
# The detected subpopulations are returned as a DataMatrix object containing one SP per row.
#
#     [finalSPs,dm]=runExPANdS(SNV, CBS) also returns the contents of the
#     output file as a DataMatrix object.
#
#     [finalSPs,dm,densities]=runExPANdS(SNV, CBS) also returns the
#     probability distributions of cellular-frequencies for each mutation.
#
#     finalSPs = runExPANdS(SNV, CBS, maxScore) keeps only SPs identified
#     at a score below maxScore (default 1).
#
#     finalSPs = runExPANdS(SNV, CBS, maxScore, max_PM) limits the maximum
#     ploidy of any locus in a mutated cell to max_PM (default 6). Increasing the value of this variable 
#     is not recommanded unless extensive depth of coverage and physical coverage underly the measurements 
#     of copy numbers and allele frequencies.
#
#     finalSPs = runExPANdS(SNV, CBS, maxScore, max_PM, precision) the
#     precision with which SP-size is predicted, a small value reflects a high resolution and can trigger
#     a higher number of predicted SPs (default 0.1/log(n/7), where n=# SNVs).
#     
#
#     finalSPs = runExPANdS(SNV, CBS, maxScore, max_PM, precision, plotF) option for displaying
#     a visual representation of the identified SPs (0 - no display; 1 - display SP size; 2 - display SP size and cell-frequency probability clusters; default: 1)
#
#     Example: SNV='Data/TCGA-06-0211-01.snvs'
#              CBS='Data/TCGA-06-0211-01.cbs'  
#              finalSPs = runExPANdS(SNV, CBS, 1, 6, 0.02, 1);
#
#     See also plotSPs
  
  
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

dirF=getwd();
##ExomeSeq
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
}

if (is.character(CBS) && file.exists(CBS)){
	copyNumber=as.matrix(read.table(CBS,sep="\t",header=TRUE))
}else if (is.matrix(CBS)){
	copyNumber=CBS;
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
   return(list("finalSPs"=NULL,"dm"=NULL,"densities"=NULL, "spGrid"=NULL));
}
maxN=8000;
if (size(dm,1)>maxN){
        idx_R=sample(1:size(dm,1), maxN, replace=F);
        dm=dm[idx_R,];
        print(paste("Input contains more than ",maxN," SNVs.", size(dm,1), " random SNVs selected for SP size calculation."))
}
cfd=computeCellFrequencyDistributions(dm, max_PM, precision);
toUseIdx=which(apply(is.finite(cfd$densities),1,all) );
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, plotF-1,label=snvF);
if(is.null(SPs) || size(SPs,1)==0){
	print("No SPs found.")
	result=list("finalSPs"=NULL,"dm"=dm,"densities"=cfd$densities);
	return(result);
}
aM= assignMutations( dm, SPs,cfd$densities);
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
    plotSPs(dm,snvF);
    ##dev.copy(jpeg,filename=paste(snvF,".sps.jpg"))
    #dev.off();
}
if (plotF>1){
    if(!require(rgl)){
	message("Plot supressed: Package rgl required for 3D plot of subpopulation clusters. Load this package before using this option.")
    }else{
    ##plot probability distributions
    cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
    kk=ceil(nrow(finalSPs)/2); par(mfcol=c(2,kk));
    for (i in 2:nrow(finalSPs)){
        idx=which(dm[,"SP"]==finalSPs[i,"Mean Weighted"]);
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

result=list("finalSPs"=finalSPs,"dm"=dm,"densities"=cfd$densities);
return(result);
}
