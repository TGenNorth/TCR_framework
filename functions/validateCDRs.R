validateCDRs=function(data){
	
###MATCH CDR1s,CDR2s,CDR2.5s###
CDRs=read.csv(file=paste0(WD,"/functions/TCRdist/CDR1s2s2.5s.csv"),head=T)

trimVname=function(invector){
	out=as.vector(invector)
	slash=grep("/",out)
	out[slash]=gsub("/.*","",out[slash])
	DV=grep("DV",out)
	out[DV]=gsub("DV.*","",out[DV])
	star=grep("\\*",out)
	out[star]=gsub("\\*.*","",out[star])
	return(out)}

CDRs$Segment=trimVname(CDRs$Segment)

data$bestVHit=trimVname(data$bestVHit)
data_CDRs=CDRs[match(data$bestVHit,CDRs[,1]),c(2,3,4)]
data$aaSeqCDR1=data_CDRs[,1]
data$aaSeqCDR2=data_CDRs[,2]
data$aaSeqCDR2.5=data_CDRs[,3]


###REMOVE UNPRODUCTIVE CDR1s,CDR2s,CDR2.5s###
censored=Reduce(union,list(which(is.na(data$aaSeqCDR1)),which(is.na(data$aaSeqCDR2)),which(is.na(data$aaSeqCDR2.5)),grep("\\*",data$aaSeqCDR3),grep("\\_",data$aaSeqCDR3),which(is.na(data$aaSeqCDR3)),which(nchar(as.character(data$aaSeqCDR3))<=5)))
selected=setdiff(1:nrow(data),censored)	
return(data[selected,])}
