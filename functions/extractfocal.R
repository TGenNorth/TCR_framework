extractfocal=function(counts,freqs,focalcolname){
overallclonotypelist=mat.or.vec(0,0)
data=cbind(freqs[,c("descrR1","class","nSeqCDR3","aaSeqCDR3","bestVHit","bestJHit")],freqs[,focalcolname],counts[,focalcolname]);colnames(data)[7]="cloneFreq";colnames(data)[8]="cloneCount"
return(data)}
