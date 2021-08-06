makecombinedshortlisttable=function(clusterlistfilename,focalsamples,pathtorawcountsfile,maxthreshold,alpha,chain){

#requires loadcountsandfreqs.R, findclusters.R, TCRdist.R

WD=paste0(WD,"workingFiles/")
setwd(WD)

overallclonotypelist=read.csv(file=paste0(WD,"/",clusterlistfilename),head=T)
loadcountsandfreqs(pathtorawcountsfile)

samplefilenames=paste0(focalsamples,"_top1000",chain,".csv")

selected=Reduce(intersect,list(which(as.numeric(as.character(overallclonotypelist[,"p.value"]))<alpha),which(overallclonotypelist$sample %in% samplefilenames)))
shortlistedclonotypes=unique(overallclonotypelist[selected,"descrR1"])
shortlist=freqs[which(freqs$descrR1 %in% shortlistedclonotypes),]

shortlist_chain=shortlist[which(shortlist$class==paste0("TR",chain)),]
chainfilename=paste0("shortlist_",chain,".csv")
write.csv(shortlist_chain,file=chainfilename,quote=F,row.names=F)
shortlistclusters=findclusters(TCRdist(chainfilename,chainfilename,chain),maxthreshold)
clusterID_pre=mat.or.vec(nrow(shortlist_chain),1); for (k in 1:length(shortlistclusters)){clusterID_pre[as.numeric(shortlistclusters[[k]])]=k}
shortlisttable=cbind(shortlist_chain[order(clusterID_pre),c("descrR1","class","nSeqCDR3","bestVHit","bestJHit","aaSeqCDR3",focalsamples)],sort(clusterID_pre))
colnames(shortlisttable)[ncol(shortlisttable)]="CET_id"
return(shortlisttable)}