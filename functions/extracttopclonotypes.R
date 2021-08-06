extracttopclonotypes=function(data,chain,topn,random=F){
output=list(0)
matchingchain=which(data$class==paste0("TR",chain))
if (random==F) {subset=matchingchain[order(-data[matchingchain,"cloneFreq"])][1:topn]}
if (random!=F) {subset=sample(matchingchain)[1:topn]}
return(data[subset,])}
