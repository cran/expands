plotSPs<-function(dm,sampleID=NA,cex=0.5){
#PLOTSPS plots coexistent subpopulations determined by ExPANdS 
#plotSPs(DM,sampleID) For each mutation (x-axis) the function displays: 
#   > the SP to which the mutation has been assigned (squares)
#   > the ploidy of the locus in that SP 
#   > the allele frequency of the mutation. 
#Allele frequencies and ploidities are colored based on the chromosome on which the mutation is located (stars ? somatic SNVs, triangles - LOH). 
#Subpopulations are colored based on the confidence with which the mutation has been assigned to the corresponding subpopulation (black - highest, white - lowest).
#DM is a DataMatrix object returned by runExPANdS. Columns must include:
#     > chr - the chrmosome on which each SNV is located
#     > AF_Tumor - the allele-frequency of each SNV
#     > PN_B - binary variable: 1 if the SNV is a germline variant, 0 if
#     somatic.
#     > SP - the subpopulation to which each mutation has been assigned as
#     fraction of cells in the tumor bulk
#     > #maxP - the confidence with which each mutation has been assigned to
#     a subpopulation - between 0 (lowest) and 100 (highest).
#     > PM - the number of copies of each mutated locus in the assigned SP
#
#plotSPs(DM,sampleID,summaryFlag)  option for displaying
#     a summary representation of the identified SPs (0 - no display; 1 -
#     display; default: 0). Each circle within the summary plot represents a SP predicted to exist
#     in the given percentage (y-axis) of the corresponding sample (x-axis). The area of the circle 
#     is proportional to the number of SNVs assigned to a SP. 
#
#   See also runExPANdS

keep=which(!is.na(dm[,"SP"]),); dm=dm[keep,];
ia=order(dm[,"%maxP"]);dm=dm[ia,];
ia=order(dm[,"SP"],decreasing = TRUE);dm=dm[ia,];
           
   maxPloidy=max(dm[,"PM"]);
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

legend1=.plotSPPerChr(dm,8,0,0,5);

x=gray.colors(100)
norm=101/length(x);
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
.plotSPPerChr(dm,17,1,0,7);
.plotSPPerChr(dm,8,0,1);
.plotSPPerChr(dm,8,1,1);
lines(c(0,nrow(dm)),c(0,0),col="black");

#ylim([-0.12*maxPloidy,1.1]);
}


.plotSPPerChr<-function(dm,lineType,lohFlag,cnFlag,sizeM){
x=rainbow(40);
legend1=list("text"=c(),"col"=c());
maxploidy=max(dm[,"PM"])+1;
for (i in 1:22){
    idx=which(dm[,"chr"]==i & ((lohFlag & dm[,"PN_B"]==1) |(!lohFlag & dm[,"PN_B"]==0)) );
    if (!isempty(idx)){
        if (!cnFlag){
            matpoints(idx,dm[idx,"AF_Tumor"],pch=lineType,col=x[i]);
        }else{
            matpoints(idx,0.1*(dm[idx,"PM"]-maxploidy),pch=20,col=x[i]);
        }
    }
    legend1$text[i]=paste("chr",i);
    legend1$col[i]=x[i];
}
return(legend1);
}
