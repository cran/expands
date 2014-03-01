assignQuantityToSP<-function(cbs, dm){
#library(matlab)
print("Assigning copy number to SPs...")
SPs=sort(unique(dm[,'SP']));
spSize=unique(round(SPs*1000)/1000);
x=colnames(cbs);
x[(length(x)+1):(length(x)+length(SPs))]=paste("SP_",as.character(spSize),sep="");
out=matrix(cbind(cbs,matrix(NaN,nrow(cbs),length(SPs))),nrow=nrow(cbs), ncol=length(x), dimnames=list(1:nrow(cbs),x))
toD=c();
##Assign copy numbers in cbs to SPs in dm
for (k in 1:nrow(out)){
    if (mod(k,100)==0){
        print(paste("Finding overlaps for CBS segment", k,"out of ",nrow(out),"..."));
    }
    idx=which(dm[,"chr"]==out[k,"chr"] & dm[,"startpos"]>=out[k,"startpos"] & dm[,"startpos"]<=out[k,"endpos"]);
    if (length(idx)==0){
        next;
    }
    for (j in 1:length(idx)){
	dmx=dm[idx[j],];
	if (is.na(dmx['SP'])){
	   next;
	}
	sp=paste("SP_",as.character(round(dmx['SP']*1000)/1000),sep="");
	if (is.na(out[k,sp]) || out[k,sp]==dmx["PM"] ){
		out[k,sp]=as.double(dmx["PM"]);
	}else{
		toD=rbind(toD, k);
		print(paste("Ambiguous ploidy of ", sp," on chr",as.character(out[k,"chr"]),":",as.character(out[k,"startpos"]),"-",as.character(out[k,"endpos"]),". Excluded"));
	}
    }
}
out=out[-toD,];
##impute missing values
#for (k in 1:nrow(out)){
#	ii= which(is.na(out[k,]));
#	for (ix in length(ii):1){
#		i=ii[ix];
#		if (i==ncol(out)){
#			out[k,i]=2;
#		}else{
#			out[k,i]=out[k,i+1];
#		}
#	}
#}

print("... Done.")
return(out);

}
