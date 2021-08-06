findCETs=function(infocallist,maxthresh=50,incounts=counts,infreqs=freqs,topn=1000,clusterlistfilename){

WD=paste0(WD,"workingFiles/")
setwd(WD)

# requires extractfocal.R, extracttopclonotypes.R, sigtestacrossthresholds.R

thresholds=seq(0,maxthresh,by=5)
overallclonotypelist=as.data.frame(mat.or.vec(0,7))
colnames(overallclonotypelist)=c("descrR1","sample","threshold","cluster_id","p.value","size","geomean")

for (chain in c("A","B")){
for (focal in infocallist){

###EXTRACT TOP CLONOTYPES###
if(length(grep("random",focal))==1){randomvalue=T;focaltoextract=gsub("_random","",focal)}
if(length(grep("random",focal))==0){randomvalue=F;focaltoextract=focal}
data=extractfocal(incounts,infreqs,focaltoextract)
topclonotypes=extracttopclonotypes(data,chain,topn,random=randomvalue)
if(length(grep("random",focal))==1){topclonotypes[,c("cloneCount","cloneFreq")]=extracttopclonotypes(data,chain,topn,random=F)[,c("cloneCount","cloneFreq")]}
file=paste0(focal,"_top",topn,chain,".csv");pathtofile=paste0(WD,file)
write.csv(topclonotypes,pathtofile,quote=F,row.names=F)

###LIST SIG CLONOTYPES ACROSS THRESHOLDS###
clonotypelist=sigtestacrossthresholds(file,chain,thresholds,topclonotypes)
if (nrow(clonotypelist)>0){overallclonotypelist=rbind(overallclonotypelist,clonotypelist)}
print(paste(focal,chain," completed"))
}}

###WRITE CLONOTYPELIST###W
write.csv(overallclonotypelist,paste0(WD,clusterlistfilename),quote=F,row.names=F)

return(overallclonotypelist)}