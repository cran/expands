assignQuantityToSP<-function (cbs, dm, colName = "PM_cnv", keepAmbigSeg = FALSE) 
{
  print("Assigning copy number to SPs...")
  SPs = sort(unique(c(dm[, "SP_cnv"],dm[,"SP"])))
  spSize = unique(round(SPs * 1000)/1000)
  x = colnames(cbs)
  x[(length(x) + 1):(length(x) + length(SPs))] = paste("SP_", 
                                                       as.character(spSize), sep = "")
  out = matrix(cbind(cbs, matrix(NaN, nrow(cbs), length(SPs))), 
               nrow = nrow(cbs), ncol = length(x), dimnames = list(1:nrow(cbs), 
                                                                   x))
  toD = c()
  allAssigned=c();
  for (k in 1:nrow(out)) {
    if (mod(k, 100) == 0) {
      print(paste("Finding overlaps for CBS segment", k, 
                  "out of ", nrow(out), "..."))
    }
    idx = which(dm[, "chr"] == out[k, "chr"] & dm[, "startpos"] >= 
                  out[k, "startpos"] & dm[, "startpos"] <= out[k, "endpos"])
    if (length(idx) == 0) {
      next
    }
    for (j in 1:length(idx)) {
      dmx = dm[idx[j], ]
      ##First assign SP with abberant copy number
      if (is.na(dmx["SP_cnv"])) {
        next
      }
      sp = paste("SP_", as.character(round(dmx["SP_cnv"] * 
                                             1000)/1000), sep = "")
      thisAssigned=c(k,dmx["SP_cnv"] ,out[k,"chr"],out[k,"startpos"],out[k,"endpos"],dmx[colName]);
      allAssigned = rbind(allAssigned, thisAssigned)
      if (is.na(out[k, sp]) || out[k, sp] == as.double(dmx[colName])) {
        out[k, sp] = as.double(dmx[colName])
      } else {
        toD = rbind(toD, thisAssigned)
      }
      
      if (is.na(dmx["SP"]) || dmx["SP"]==dmx["SP_cnv"]) {
        next
      }
      ##Next assign SP with copy number, if co-occurence assumption violation applies 
      if (is.na(dmx["PM"]) ){
        next
      }
      sp = paste("SP_", as.character(round(dmx["SP"] *
                                             1000)/1000), sep = "")
      thisAssigned=c(k,dmx["SP"] ,out[k,"chr"],out[k,"startpos"],out[k,"endpos"],dmx["PM"]);
      allAssigned = rbind(allAssigned, thisAssigned)
      if (is.na(out[k, sp]) || out[k, sp] == dmx["PM"]) {
        out[k, sp] = as.double(dmx["PM"])
      } else {
        toD = rbind(toD, thisAssigned)
      }
      
    }
  }
  
  ##Either remove ambiguous segments or calculate their median ploidy based on SNV ploidy 
  printErr=FALSE;
  if (!is.null(nrow(toD))) {
    toD=matrix(toD,nrow=nrow(toD),ncol=ncol(toD),dimnames=list(
      rows=c(1:nrow(toD)),cols=c("Idx","SP","chr","startpos","endpos",colName)))
    colnames(allAssigned)=colnames(toD)  
    uD=unique(toD);
    for (i in 1:nrow(uD)){
      ii=which(as.numeric(uD[i,"SP"])==as.numeric(allAssigned[,"SP"]) & 
                 as.numeric(uD[i,"chr"])==as.numeric(allAssigned[,"chr"]) &
                 as.numeric(uD[i,"startpos"])==as.numeric(allAssigned[,"startpos"]) &
                 as.numeric(uD[i,"endpos"])==as.numeric(allAssigned[,"endpos"]));
      sp = paste("SP_", as.character(round(as.numeric(uD[i,"SP"] )* 1000)/1000), sep = "")
      tmpI=as.numeric(uD[i,"Idx"]);
      if (!keepAmbigSeg){
        out[tmpI,sp]=NA;
        printErr=TRUE;    
      }else{
        out[tmpI,sp]=round(median(as.numeric(allAssigned[ii,colName])));
      }
    }
  }else if (is.null(nrow(toD)) && length(toD) > 0) {
    if (!keepAmbigSeg){
      out[as.numeric(toD[1]),]=NA;
      printErr=TRUE;    
    }
  }
  
  if(printErr){
    print(paste("Ambiguous SP specific ploidies found for ",length(toD)," segment-SP pairs."));
    print("Ploidies not assigned for these segments in corresponding SPs.")
  } 
  
  print("... Done.")
  if (keepAmbigSeg){
    print("Warning: parameter <keepAmbigSeg> set to TRUE. Output includes segment-assignements where subpopulation specific ploidy is ambiguous.Recommend repeating circular binary segmentation with less stringent parameters instead, to reduce segment length and thus the prevalence of ambiguous assignements.")
  }
  
  return(out)
}
