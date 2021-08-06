getseqlistav=function(topclonotypes,clusters){
seqlist=list(1);av=list(1);for (m in 1:length(clusters)){if(length(clusters[[m]])>1) {seqlist[[m]]=topclonotypes[as.numeric(clusters[[m]]),]};av[[m]]=mean(unlist(log(seqlist[[m]]["cloneFreq"],10)))}

out=list("")
out[[1]]=seqlist
out[[2]]=av
return(out)}
