CETsigtest=function(topclonotypes,clusters,clustersizetrials,numberoftrials=10000){

#requires getseqlistav.R
	
#### MODELS #OCCURRENCES OF EACH SIZE PER TRIAL IN A SAMPLE-SPECIFIC WAY
outall=list(NULL)
for (trials in 1:numberoftrials){
	set=count(unlist(clustersizetrials[sample(1:(length(clustersizetrials)),1)]))

if (nrow(set)!=0){
	subsetfreqs=topclonotypes[,"cloneFreq"]
	out=mat.or.vec(nrow(set),2)
	for (row in 1:nrow(set)){
		values=mat.or.vec(1,set[row,2])
		for (samplenumber in 1:set[row,2]){values[samplenumber]=mean(log(sample(subsetfreqs,set[row,1]),10))}
	out[row,1]=set[row,1]
	out[row,2]=max(values)
}
outall[[trials]]=out}
}

av=getseqlistav(topclonotypes,clusters)[[2]]

clusterlist=mat.or.vec(length(clusters),3)
clusterlist[,1]=unlist(lengthmean(clusters)[,1])
clusterlist[,2]=unlist(av)
colnames(clusterlist)=c("size","geomean","p-value")
for (n in 1:nrow(clusterlist)){
	X=clusterlist[n,1]
	Y=clusterlist[n,2]
	pcount=0
	for (k in 1:length(outall)){
		if (length(outall[[k]])>0){
			for (l in 1:nrow(outall[[k]])){
			if (outall[[k]][l,1]==X && outall[[k]][l,2]>Y)
			pcount=pcount+1
			}}
	}
clusterlist[n,3]=pcount/numberoftrials}

	return(clusterlist)
}
