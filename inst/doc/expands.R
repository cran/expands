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
SPs=clusterCellFrequencies(cfd$densities[toUseIdx,], precision, label=snvF)


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
### code chunk number 8: expands.Rnw:99-100
###################################################
plotSPs(aM$dm, snvF)


