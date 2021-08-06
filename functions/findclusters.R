findclusters=function(pathorobject,threshold,names="ranks",input="object"){
if (input=="file") {dists=read.csv(file=pathorobject,row.names=1,head=T)}
if (input=="object") {dists=pathorobject}
if (names=="ranks") {row.names(dists)=1:nrow(dists);colnames(dists)=1:ncol(dists)}
clusteredobject=hclust(as.dist(dists))
groups=cutree(clusteredobject,h=threshold)
clusterindeces=count(groups)[which(count(groups)[,2]>1),1]
clusters=as.list(rep(0,length(clusterindeces)))
if(length(clusterindeces)!=0){for (n in 1:length(clusterindeces)) {m=clusterindeces[n];clusters[[n]]=names(groups)[which(groups==m)]}}
return(clusters)}

lengthmean=function(clusters){
out=mat.or.vec(length(clusters),2)
colnames(out)=c("length","mean")
if (length(clusters)>0) {for (n in 1:length(clusters)){
out[n,1]=length(as.numeric(clusters[[n]]))
out[n,2]=mean(as.numeric(clusters[[n]]))}}
return(out)}
