loadcountsandfreqs=function(pathtofile){

#requires annotatenamechain.R, validateCDRs.R

###LOAD DATA###
counts_pre=read.csv(file=pathtofile,head=T)
if(length(grep("descrR1",colnames(counts_pre)))==0) {counts=annotatenamechain(validateCDRs(counts_pre[which(counts_pre$class!="other"),]))}
if(length(grep("descrR1",colnames(counts_pre)))==1) {counts=validateCDRs(counts_pre[which(counts_pre$class!="other"),])}

###NORMALIZE COUNTS TO FREQS###
freqs=counts
for (column in setdiff(colnames(counts),
c("nSeqCDR3","bestVHit","bestJHit","aaSeqCDR3","class","aaSeqCDR1","aaSeqCDR2","aaSeqCDR2.5","descrR1"))){
cloneFreq=0
cloneFreq[which(counts$class=="TRA")]=counts[which(counts$class=="TRA"),column]/sum(counts[which(counts$class=="TRA"),column])
cloneFreq[which(counts$class=="TRB")]=counts[which(counts$class=="TRB"),column]/sum(counts[which(counts$class=="TRB"),column])
freqs[,column]=cloneFreq}

counts<<-counts
freqs <<-freqs}
