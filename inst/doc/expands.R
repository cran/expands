### R code from vignette source 'expands.Rnw'

###################################################
### code chunk number 1: expands.Rnw:43-52
###################################################
library(expands)
##load mutations:
data(snv);
## use only a subset of mutations (to reduce time required to run example):
set.seed(6); idx=sample(1:nrow(snv), 80, replace=FALSE); snv=snv[idx,];
##load copy number segments:
data(cbs);
##assign copy numbers to point mutations:
dm=assignQuantityToMutation(snv,cbs,"CN_Estimate");


###################################################
### code chunk number 2: expands.Rnw:57-62
###################################################
##parameters
max_PM=6; maxS=0.7; precision=0.018;
plotF=1; 
##the name of the sample
snvF="TCGA-06-0152-01";


###################################################
### code chunk number 3: expands.Rnw:68-70
###################################################
##calculate cell frequency probability distribution for each mutation
cfd=computeCellFrequencyDistributions(dm, max_PM, p=precision)


###################################################
### code chunk number 4: expands.Rnw:73-75
###################################################
##cluster mutations with valid distributions
toUseIdx=which(apply(is.finite(cfd$densities),1,all) )


###################################################
### code chunk number 5: expands.Rnw:81-83
###################################################
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], p=precision)
SPs=SPs[SPs[,"score"]<=maxS,]; ## exclude SPs detected at high noise levels


###################################################
### code chunk number 6: expands.Rnw:86-87
###################################################
print(SPs)


###################################################
### code chunk number 7: expands.Rnw:91-93
###################################################
##assign mutations to subpopulations:
aM= assignMutations( dm, SPs, verbose = F)


###################################################
### code chunk number 8: expands.Rnw:103-104
###################################################
o=plotSPs(aM$dm, snvF,cex=1)


###################################################
### code chunk number 9: expands.Rnw:115-117
###################################################
##assigning copy number to subpopulations
aQ=assignQuantityToSP(cbs, aM$dm, v=F)


###################################################
### code chunk number 10: expands.Rnw:120-122
###################################################
##building phylogeny
spPhylo=buildPhylo(aQ,snvF,add = NULL)


###################################################
### code chunk number 11: expands.Rnw:127-128
###################################################
plot(spPhylo$tree,cex=3,type = "c")


