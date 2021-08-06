TCRdist=function(filename1,filename2,AorB="AB"){

#requires annotatenamechain.R, validateCDRs.R


###IMPORT DATA AND VALIDATE###
print(paste("***********started",Sys.time()))
'%!in%' <- function(x,y)!('%in%'(x,y))

renamecols=function(data){
if(length(which(colnames(data)=="v_gene")==1)) {colnames(data)[which(colnames(data)=="v_gene")]="bestVHit"}
if(length(which(colnames(data)=="j_gene")==1)) {colnames(data)[which(colnames(data)=="j_gene")]="bestJHit"}
if(length(which(colnames(data)=="V")==1)) {colnames(data)[which(colnames(data)=="V")]="bestVHit"}
if(length(which(colnames(data)=="J")==1)) {colnames(data)[which(colnames(data)=="J")]="bestJHit"}
if(length(which(colnames(data)=="chain")==1)) {colnames(data)[which(colnames(data)=="chain")]="class"}
if(length(which(colnames(data)=="cdr3")==1)) {colnames(data)[which(colnames(data)=="cdr3")]="aaSeqCDR3"}
if(("descrR1" %!in% colnames(data)) && ("clonotype_id" %!in% colnames(data))){data=annotatenamechain(data)}
return(data)}

data1_prevalidated=renamecols(read.csv(file=filename1,header=TRUE)); data2_prevalidated=renamecols(read.csv(file=filename2,header=TRUE))
data1_prevalidated[is.na(data1_prevalidated)]=""; data2_prevalidated[is.na(data2_prevalidated)]=""

print(paste("*********data_read",Sys.time()))

data1=validateCDRs(data1_prevalidated)
data2=validateCDRs(data2_prevalidated)

print(paste("*********validated",Sys.time()))

###LIST ALL COMBINATIONS OF ab CLONOTYPES IN FILE1 and FILE2###

findpairs=function(data){

if ("clonotype_id" %!in% colnames(data)){
split=strsplit(as.character(data$descrR1),"_")
clones=mat.or.vec(length(split),3)
for (n in 1:length(split)){clones[n,1]=unlist(split[[n]])[1]; clones[n,2]=paste(unlist(split[[n]])[2:length(unlist(split[[n]]))],collapse="_")}
clones[,3]=as.character(data$class)}

if ("clonotype_id" %in% colnames(data)) {clones=cbind(as.character(data$clonotype_id),as.character(data$consensus_id),as.character(data$class))}

X=tapply(seq_along(clones[,1]), clones[,1], identity)[unique(clones[,1])]
pairs=mat.or.vec(length(X),3)
for (k in 1:length(X)){
	family=X[k][[1]]
	if(paste(clones[family,3],collapse="")=="TRATRB") {pairs[k,]=c(names(X[k]),rev(rev(family)))
	} else if(paste(clones[family,3],collapse="")=="TRBTRA") {pairs[k,]=c(names(X[k]),rev(family))
	} else pairs[k,2]="not1A1B"}
	
pairs=pairs[which(pairs[,2]!="not1A1B"),]
colnames(pairs)=c("clonotype","TRA","TRB")
return(pairs)}

#numberofrows=(nrow(data1)-1)/2 ; X=1:numberofrows
#pairsin1=cbind(X,X+1,X+numberofrows);colnames(pairsin1)=c("clonotype","TRA","TRB")
#print(paste("skipping_findpairs",Sys.time()))

if (AorB=="AB"){pairsin1=findpairs(data1);pairsin2=findpairs(data2)}
if (AorB=="A") {pairsin1=cbind(as.character(data1$descrR1[which(data1$class=="TRA")]),which(data1$class=="TRA"),which(data1$class=="TRA")); colnames(pairsin1)=c("clonotype","TRA","TRA")
				pairsin2=cbind(as.character(data2$descrR1[which(data2$class=="TRA")]),which(data2$class=="TRA"),which(data2$class=="TRA")); colnames(pairsin2)=c("clonotype","TRA","TRA")}	
if (AorB=="B") {pairsin1=cbind(as.character(data1$descrR1[which(data1$class=="TRB")]),which(data1$class=="TRB"),which(data1$class=="TRB")); colnames(pairsin1)=c("clonotype","TRB","TRB")
		        pairsin2=cbind(as.character(data2$descrR1[which(data2$class=="TRB")]),which(data2$class=="TRB"),which(data2$class=="TRB")); colnames(pairsin2)=c("clonotype","TRB","TRB")}

print(paste("*******found_pairs",Sys.time()))

###PREPARE FOR LOOP, INCLUDING DEFINING FUNCTIONS###

N=nrow(pairsin1); M=nrow(pairsin2)
allpaircombinations=mat.or.vec(M*N,6)

for (n in 1:N) {for (m in 1:M) {allpaircombinations[(n-1)*M+m,]=append(pairsin1[n,],pairsin2[m,])}}

output=mat.or.vec(M*N,1)

print(paste(M*N,"comparisons",Sys.time()))

trim=function(string){return(substring(string,4,nchar(as.character(string))-2))}


### WORDDIST FUNCTION ###


if (AorB=="AB" || AorB=="A"){
	
CDR1a_first=data1$aaSeqCDR1[as.numeric(allpaircombinations[,2])]
CDR2a_first=data1$aaSeqCDR2[as.numeric(allpaircombinations[,2])]
CDR2.5a_first=data1$aaSeqCDR2.5[as.numeric(allpaircombinations[,2])]
CDR3a_first=trim(data1$aaSeqCDR3[as.numeric(allpaircombinations[,2])])
CDR1a_second=data2$aaSeqCDR1[as.numeric(allpaircombinations[,5])]
CDR2a_second=data2$aaSeqCDR2[as.numeric(allpaircombinations[,5])]
CDR2.5a_second=data1$aaSeqCDR2.5[as.numeric(allpaircombinations[,5])]
CDR3a_second=trim(data2$aaSeqCDR3[as.numeric(allpaircombinations[,5])])

write.table(cbind(as.character(CDR1a_first),as.character(CDR1a_second)),file="W_CDR1a.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR2a_first),as.character(CDR2a_second)),file="W_CDR2a.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR2.5a_first),as.character(CDR2.5a_second)),file="W_CDR2.5a.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR3a_first),as.character(CDR3a_second)),file="W_CDR3a.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
}

if (AorB=="AB" || AorB=="B"){

CDR1b_first=data1$aaSeqCDR1[as.numeric(allpaircombinations[,3])]
CDR2b_first=data1$aaSeqCDR2[as.numeric(allpaircombinations[,3])]
CDR2.5b_first=data1$aaSeqCDR2.5[as.numeric(allpaircombinations[,3])]
CDR3b_first=trim(data1$aaSeqCDR3[as.numeric(allpaircombinations[,3])])
CDR1b_second=data2$aaSeqCDR1[as.numeric(allpaircombinations[,6])]
CDR2b_second=data2$aaSeqCDR2[as.numeric(allpaircombinations[,6])]
CDR2.5b_second=data2$aaSeqCDR2.5[as.numeric(allpaircombinations[,6])]
CDR3b_second=trim(data2$aaSeqCDR3[as.numeric(allpaircombinations[,6])])

write.table(cbind(as.character(CDR1b_first),as.character(CDR1b_second)),file="W_CDR1b.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR2b_first),as.character(CDR2b_second)),file="W_CDR2b.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR2.5b_first),as.character(CDR2.5b_second)),file="W_CDR2.5b.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
write.table(cbind(as.character(CDR3b_first),as.character(CDR3b_second)),file="W_CDR3b.csv",row.names=FALSE,quote=FALSE,col.names=FALSE,sep=",")
}

print(paste("***written infiles",Sys.time()))

if (AorB=="AB" || AorB=="A"){
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR1a.csv -o W_CDR1a_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR2a.csv -o W_CDR2a_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR2.5a.csv -o W_CDR2.5a_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR3a.csv -o W_CDR3a_dist.csv -g 8 -m 2"), wait = TRUE)
}

if (AorB=="AB" || AorB=="B"){
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR1b.csv -o W_CDR1b_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR2b.csv -o W_CDR2b_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR2.5b.csv -o W_CDR2.5b_dist.csv -g 4 -m 2"), wait = TRUE)
system(paste0(WD,"functions/TCRdist/worddist -i W_CDR3b.csv -o W_CDR3b_dist.csv -g 8 -m 2"), wait = TRUE)
}

print(paste("computed worddists",Sys.time()))

if (AorB=="AB" || AorB=="A"){
CDR1a_dist=read.table(file="W_CDR1a_dist.csv",sep=",")
CDR2a_dist=read.table(file="W_CDR2a_dist.csv",sep=",")
CDR2.5a_dist=read.table(file="W_CDR2.5a_dist.csv",sep=",")
CDR3a_dist=read.table(file="W_CDR3a_dist.csv",sep=",")
}

if (AorB=="AB" || AorB=="B"){
CDR1b_dist=read.table(file="W_CDR1b_dist.csv",sep=",")
CDR2b_dist=read.table(file="W_CDR2b_dist.csv",sep=",")
CDR2.5b_dist=read.table(file="W_CDR2.5b_dist.csv",sep=",")
CDR3b_dist=read.table(file="W_CDR3b_dist.csv",sep=",")
}
print(paste("*****read outfiles",Sys.time()))

system("rm W_CDR*", wait = TRUE)

if(AorB=="AB"){TCRdistance = CDR1a_dist + CDR2a_dist + CDR2.5a_dist + 3*CDR3a_dist + CDR1b_dist + CDR2b_dist + CDR2.5b_dist + 3*CDR3b_dist}
if(AorB=="A"){TCRdistance = CDR1a_dist + CDR2a_dist + CDR2.5a_dist + 3*CDR3a_dist}
if(AorB=="B"){TCRdistance = CDR1b_dist + CDR2b_dist + CDR2.5b_dist + 3*CDR3b_dist}

print(paste("calculated TCRdist",Sys.time()))

outputmatrix=matrix(TCRdistance[,1],nrow=N,ncol=M,byrow=TRUE)
row.names(outputmatrix)= pairsin1[,1]
colnames(outputmatrix)= pairsin2[,1]

print(paste("**reshaped TCRdist",Sys.time()))


filename1prefix=gsub(".csv","",filename1,perl=TRUE); filename2prefix=gsub(".csv","",filename2,perl=TRUE)
write.csv(outputmatrix,paste(filename1prefix,"_",filename2prefix,"_dist",AorB,".csv",sep=""),quote=FALSE)

print(paste("*********finished!",Sys.time()))
return(outputmatrix)
}

