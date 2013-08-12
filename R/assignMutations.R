assignMutations<-function( dm, finalSPs,densities, max_PM=6){
if(is.null(dim(finalSPs))){
	precision=finalSPs["x_p"];
	spFreq=finalSPs["Mean Weighted"];
}else{
	spFreq=finalSPs[,"Mean Weighted"];
	precision=finalSPs[1,"x_p"];
}
freq=t(seq(0.1,1.1,by=precision/10));
##Assign mutations to SPs
idx=matrix(NA,1,size(finalSPs,1));
for (i in 1:size(finalSPs,1)){
    idx[i]=which.min(abs(freq-spFreq[i]));
}
    
addCols=c("%maxP","SP","PM","PM_B");
for (k in 1:length(addCols)){
    dm=.addColumn(dm,addCols[k],NA);
}
if (!any(colnames(dm)=="f")){
    dm=.addF(dm,  max_PM);
}
dm[,"SP"]=NA; ##delete any potentially existing SP info

for(k in 1:nrow(dm)){
    #     [a,ia]=max(densities(k,idx));
    #     p_fadj=densities(k,idx(ia));
    if (is.na(dm[k,"f"]) || dm[k,"f"]<0.05){
        next;
    }
    ia=which.min(abs(spFreq-dm[k,"f"]));
    p_fadj=densities[k,idx[ia]];
    
    pm=(dm[k,"CN_Estimate"]-(1-spFreq[ia])*2)/spFreq[ia];
    pm=max(1,pm); pm=min(max_PM,pm);
    pmb=(dm[k,"CN_Estimate"]*dm[k,"AF_Tumor"]-(1-spFreq[ia])*dm[k,"PN_B"])/spFreq[ia];
    pmb=min(pm,pmb);pmb=max(1,pmb);
    maxP=max(densities[k,]);
    xxx=100*p_fadj/maxP;
    dm[k,c("SP","%maxP")]=c(spFreq[ia],xxx);
    dm[k,c("PM","PM_B")]=c(round(pm),round(pmb));
}
dm[is.na(dm[,"%maxP"]),"%maxP"]=0;

for (j in 1:size(finalSPs,1)){
	if(is.null(dim(finalSPs))){
	    idx=which(dm[,"SP"]==finalSPs["Mean Weighted"]);
	    finalSPs["nMutations"]=length(idx);	  
	}else{
	    idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
	    finalSPs[j,"nMutations"]=length(idx);
	}
}
output=list("dm"=dm,"finalSPs"=finalSPs);
return(output);
}


.addF<-function (dm,  max_PM){
#.jaddClassPath("ExPANydS.jar")
.jinit(classpath="ExPANdS.jar")
#javaImport(packages = "core.analysis.ngs.algorithms.*")
dm=.addColumn(dm,"f",NA);

for (k in 1:nrow(dm)){
   expands <-try(.jnew("ExPANdS", as.double(dm[k,"AF_Tumor"]),as.double(dm[k,"CN_Estimate"]),
                  as.integer(dm[k,"PN_B"]),as.integer(max_PM)));
    if (class(expands)=="try-error"){
	print(expands);
	print(paste("At SNV ",k,": -->"));
	print(dm[k,]);
    }else{
	.jcall(expands,,"run")
    	dm[k,"f"]<-.jcall(expands,"D","getF");
    }
}
return(dm);
}
