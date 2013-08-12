.addColumn<-function(M,newCol,initVal){
  if (!any(colnames(M)==newCol)){
    M=matrix(cbind(M,matrix(initVal,nrow(M),1)),nrow=nrow(M),ncol=ncol(M)+1,
             dimnames = list(rownames(M), c(colnames(M),newCol)));
  }
  return(M);
}
