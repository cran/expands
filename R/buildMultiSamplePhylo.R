buildMultiSamplePhylo<-function (samGr, out, treeAlgorithm = "bionjs", keepAmbigSeg = FALSE, plotF=1) 
{
  library(expands)
  cols = c("Count", "chr", "startpos", "endpos",  "CN_Estimate")
  dummySNVcols=c("Count","endpos");
  allCBS = c()
  allDM = c()
  dmPris=list();
  n_Samples = length(samGr$labels)
  for (i in 1:n_Samples) {
    cbsPri = read.table(samGr$cbs[[i]], sep = "\t", 
                        header = T)
    cbsPri[, c("startpos", "endpos")] = round(cbsPri[, c("startpos", 
                                                         "endpos")]/10000) * 10000
    dmPri = read.table(samGr$sps[[i]], sep = "\t", 
                       header = T) 
    for (j in 1:length( dummySNVcols)){
      if (!any(colnames(dmPri)==dummySNVcols[j])){
        tmp=colnames(dmPri);
        dmPri=cbind(dmPri,matrix(NA,nrow(dmPri),1));
        colnames(dmPri)=c(tmp,dummySNVcols[j]);
        if(dummySNVcols[j]=='endpos'){
          dmPri[,dummySNVcols[j]]=dmPri[,'startpos'];
        }
      }
    }
    dmPris[[i]]=dmPri;
    allCBS = as.matrix(rbind(allCBS, cbsPri))
    allCBS[, "Count"] = c(1:nrow(allCBS))
    allDM = as.matrix(rbind(allDM, dmPri))
    allDM[, "Count"] = c(1:nrow(allDM))
  }
  dupI = which(duplicated(allCBS[, c("chr", "startpos", "endpos")]))
  if (length(dupI) > 0) {
    allCBS = allCBS[-1 * dupI, ]
  }
  dupI = which(duplicated(allDM[, c("chr", "startpos")]))
  if (length(dupI) > 0) {
    allDM = allDM[-1 * dupI, ]
  }
  aqCBS = allCBS[, cols]
  aqDM = allDM[, cols]
  for (i in 1:n_Samples) {
    dmPri = dmPris[[i]];
    aQpriCBS = try(assignQuantityToSP(allCBS[, cols], dmPri, 
                                      keepAmbigSeg = keepAmbigSeg), silent = FALSE)
    dmPri[, "PM_B"] = sign(dmPri[, "PM_B"])
    aQpriDM = try(assignQuantityToSP(allDM[, cols], dmPri, 
                                     colName = "PM_B"), silent = FALSE)
    firstI = min(grep("SP", colnames(aQpriCBS)))
    aqCBS = cbind(aqCBS, aQpriCBS[, firstI:ncol(aQpriCBS)])
    aqDM = cbind(aqDM, aQpriDM[, firstI:ncol(aQpriDM)])
    nSPs = length(unique(dmPri[!is.na(dmPri[, "SP"]), "SP"]))
    lab = paste(samGr$labels[[i]], "_SP", sep = "")
    colns = colnames(aQpriDM)
    colnames(aqCBS) = c(colnames(aqCBS[, 1:(ncol(aqCBS) - 
                                              nSPs)]), gsub("SP", lab, colns[firstI:ncol(aQpriDM)]))
  }
  aQ = rbind(aqCBS, aqDM)
  tr = NULL
  if (class(aQ) == "try-error" || is.null(ncol(aQ))) {
    print("Error encountered while reconstructing phylogeny")
  }
  else {
    tr = buildPhylo(aQ, out, treeAlgorithm = treeAlgorithm)
    if (plotF>0){
   	jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                              "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
   	colmap = jet(n_Samples)
    	colors <- rep(colmap[1], each = length(tr$tip.label))
    	for (i in 1:n_Samples) {
    	  ii = grep(samGr$labels[[i]], tr$tip.label)
    	  colors[ii] = colmap[i]
    	}
    	plot(tr, tip.col = colors, cex = 1, type = "u")
    }
    return(tr)
  }
  return(NULL)
}
