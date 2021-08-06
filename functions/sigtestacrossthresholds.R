sigtestacrossthresholds=function(file,chain,thresholds,topclonotypes){

#requires extracttopclonotypes.R, TCRdist.R, findclusters.R, getseqlistav.R, CETsigtest.R

load(paste0(WD,"functions/randomtrials/100_randomtrials_chain=A.Rdat"))
load(paste0(WD,"functions/randomtrials/100_randomtrials_chain=B.Rdat"))
setwd(paste0(WD,"workingFiles/"))

###CALCULATE DISTS/CLUSTERS###
dists=TCRdist(file,file,chain)
outlist=as.data.frame(mat.or.vec(0,6))
colnames(outlist)=c("descrR1","sample","threshold","cluster_id","p-value","geomean")

for (threshold in thresholds){
clusters=findclusters(dists,threshold)

if (length(clusters)>0) {
seqlist=getseqlistav(topclonotypes,clusters)[[1]]

###TEST SIGNIFICANCE###
if (chain=="A"){randomtrialsvector=randomtrialsA[,paste0("threshold_",threshold)]}
if (chain=="B"){randomtrialsvector=randomtrialsB[,paste0("threshold_",threshold)]}
clusterlist=CETsigtest(topclonotypes,clusters,clustersizetrials=randomtrialsvector,numberoftrials=1000)
print(paste0(length(clusters)," clusters at threshold ",threshold," for ",file))


A=as.numeric(unlist(lapply(seqlist,function(X) X[,"descrR1"])))
B=rep(file,length(A))
C=as.numeric(rep(threshold,length(unlist(lapply(seqlist,function(X) X[,"descrR1"])))))
D=paste(B,C,rep(1:length(seqlist)[1],clusterlist[,"size"]),sep="_")
E=as.numeric(rep(clusterlist[,"p-value"],clusterlist[,"size"]))
F=as.numeric(rep(clusterlist[,"size"],clusterlist[,"size"]))
G=as.numeric(rep(clusterlist[,"geomean"],clusterlist[,"size"]))

newoutlist=as.data.frame(cbind(A,B,C,D,E,F,G))
colnames(newoutlist)=c("descrR1","sample","threshold","cluster_id","p-value","size","geomean")
outlist=rbind(outlist,newoutlist)}
}

return(outlist)}


