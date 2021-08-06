annotatenamechain=function(decomp){
As=intersect(grep("TRAV",decomp$bestVHit),grep("TRAJ",decomp$bestJHit))
Bs=intersect(grep("TRBV",decomp$bestVHit),grep("TRBJ",decomp$bestJHit))
Cs=setdiff(1:nrow(decomp),union(As,Bs))
class=mat.or.vec(nrow(decomp),1)
if (length(As)>0) class[As]="TRA"
if (length(Bs)>0) class[Bs]="TRB"
if (length(Cs)>0) class[Cs]="other"
descrR1=1:nrow(decomp)
return(cbind(descrR1,class,decomp))}
