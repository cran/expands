cellfrequency_pdf <-function(af,cnv,pnb,freq, max_PM=6){
.jaddClassPath("ExPANydS.jar")
.jinit(classpath="ExPANdS.jar")
javaImport(packages = "core.analysis.ngs.algorithms.*")
##Compute all possible solutions for variable PM, PM_B and e
expands <-.jnew("ExPANdS", as.double(af),as.double(cnv),as.integer(pnb),as.integer(max_PM));
.jcall(expands,,"run")
bestF<-.jcall(expands,"D","getF");
fits<-.jcall(expands,"Ljava/util/Collection;","solutions");
fits<-.jcall(fits,"Ljava/util/Iterator;","iterator");
results=c();
while (.jcall(fits,"Z","hasNext")){
    solutions=.jcall(fits,"Ljava/lang/Object;","next")
    solutions=.jcall(solutions,"Ljava/util/Iterator;","iterator");
    while (.jcall(solutions,"Z","hasNext")){
        rawR<-.jcall(solutions,"Ljava/lang/Object;","next")
        rawR<-.jcall(rawR,"[D","toDouble");
        results=rbind(results,rawR);
    }
}

##Keep only valid solutions
fit=matrix(results, nrow = nrow(results), ncol = ncol(results), 
           dimnames = list(paste(1:nrow(results)), .jfield(expands,,"SOLUTION_ENTITIES")))
fit=fit[fit[,"f"]>=0.1 & fit[,"f"]<=1.1,];

z=round(fit[,"f"]*100);
z1=sort(unique(z));
dm=matrix(nrow=length(z1), ncol=ncol(fit), 
          dimnames=list(paste(1:length(z1)),colnames(fit)));
tfit=t(fit);
for (i in 1:length(z1)){
    f=z1[i];
    similarFrequencies=t(tfit[,z==f]);
    ia=which.min(similarFrequencies[,"dev"]);
    dm[i,]=similarFrequencies[ia,];
}

##create frequency array weghted by deviation
 normdev=dm[,"dev"];
 normdev=round(-100*.sigmoid(normdev*50, 2, 4)+101);
# normdev=round(-1*log10(normdev/max(normdev)))+1;
 #normdev=1+round(100*(1-normdev));
# normdev=round(100/double(dm(:,"dev")));
f=matrix(NA,sum(normdev),1);
for (i in 1:nrow(dm)){
    if (i==1){
        start = 1;
    }else{
        start=sum(normdev[1:i-1])+1;
    }
    idx=start:sum(normdev[1:i]);
    f[idx]=dm[i,"f"];
}

nComponents=ceil(fit[1,"CN_Estimate"])+1; ##number of components in gaussian mixture model
p=matrix(NA,1,length(freq));

obj=suppressWarnings(densityMclust(f[!is.na(f)],G=max(1,nComponents-1):nComponents));
stowarn<-warnings();
for (w in names(stowarn)){
        if ( isempty(agrep("occurs at min or max",w,max.distance=0.3)) ){
                 warning(as.character(stowarn[[w]]),": ",w)
        }
}
p=predict(obj, t(freq));  p=p/sum(p)
output=list("p"=p,"bestF"=bestF,"normdev"=normdev,"fit"=fit,"errors"=NULL,"f"=f);
return(output)
}

.sigmoid<-function(x, t1, t2){
res = 1/(1 + exp(-t1*(x-t2)));
return(res);
}

