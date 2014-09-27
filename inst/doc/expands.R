### R code from vignette source 'expands.Rnw'

###################################################
### code chunk number 1: expands.Rnw:40-49
###################################################
library(expands)
##loading mutations
data(snv);
## use only a subset of all mutations (for performance reasons).
set.seed(6); idx=sample(1:nrow(snv), 130, replace=FALSE); snv=snv[idx,];
##loading copy number segments
data(cbs);
##assign copy number to mutations
dm=assignQuantityToMutation(snv,cbs,"CN_Estimate");


###################################################
### code chunk number 2: expands.Rnw:54-59
###################################################
##parameters
max_PM=6; maxScore=2.5; precision=0.018;
plotF=1; 
##the name of the sample
snvF="TCGA-06-0152-01";


###################################################
### code chunk number 3: expands.Rnw:65-67
###################################################
##compute the cell frequency probability distribution for each mutation
cfd=computeCellFrequencyDistributions(dm, max_PM, precision)


###################################################
### code chunk number 4: expands.Rnw:70-72
###################################################
##cluster mutations with valid distributions
toUseIdx=which(apply(is.finite(cfd$densities),1,all) )


###################################################
### code chunk number 5: expands.Rnw:78-79
###################################################
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision)


###################################################
### code chunk number 6: expands.Rnw:82-83
###################################################
print(SPs)


###################################################
### code chunk number 7: expands.Rnw:88-90
###################################################
##assign mutations to subpopulations
aM= assignMutations( dm, SPs,cfd$densities)


###################################################
### code chunk number 8: expands.Rnw:100-101
###################################################
plotSPs(aM$dm, snvF,cex=1)


###################################################
### code chunk number 9: expands.Rnw:112-114
###################################################
##assigning copy number to subpopulations
aQ=assignQuantityToSP(cbs, aM$dm)


###################################################
### code chunk number 10: expands.Rnw:117-119
###################################################
##building phylogeny
tr=buildPhylo(aQ,snvF)


###################################################
### code chunk number 11: expands.Rnw:124-125
###################################################
plot(tr,cex=2.5)


###################################################
### code chunk number 12: expands.Rnw:136-144
###################################################
#Patient and sample labels
patient='ID_MRD_001';
samples=c('_primPancreas','_metKidney','_metLung');
output=patient;
#The CBS files for each sample:
cbs=as.list(paste(patient, samples,'.cbs',sep=""));
#The SP files for each sample (previously calculated via runExPANdS-function):
sps=as.list(paste(patient, samples,'.sps',sep=""));


###################################################
### code chunk number 13: expands.Rnw:147-158
###################################################
sampleGroup=list(cbs=cbs,sps=sps,labels=samples)
tr=buildMultiSamplePhylo(sampleGroup,output,keepAmbigSeg = TRUE, plotF=0);
##Tree tip color labels according to sample origin of SPs:
jet <- colorRampPalette(c("#00007F", "blue", "#007FFF",
            "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colmap = jet( length(sampleGroup$labels) )
colors <- rep(colmap[1], each = length(tr$tip.label))
for (i in 1: length(sampleGroup$labels) ) {
     ii = grep(sampleGroup$labels[[i]], tr$tip.label)
     colors[ii] = colmap[i]
}


###################################################
### code chunk number 14: expands.Rnw:162-163
###################################################
plot(tr, tip.col = colors, cex = 1.6, type = "u")


