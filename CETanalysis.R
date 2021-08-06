WD="/Users/jaltin/Dropbox/Work-TGen/Projects/TCRPepSeq/scripts/CETcode/"			#set the directory of the CETcode folder

###SOURCE FUNCTIONS###
library(plyr)
setwd(paste0(WD,"functions/"))
files.sources=list.files()[grep("\\.R",list.files())]
invisible(sapply(files.sources,source))

setwd(paste0(WD,"workingFiles/"))
pathtocountsfile=paste0(WD,"countsTables/FluM1.csv")								#specify the input datafile
clusterlistfilename="clusterlistfile.csv"											#name the new file that will contain the clusterlist

###IMPORT COUNTS###
loadcountsandfreqs(pathtocountsfile)

###GENERATE CLUSTERLIST###
focallist=c("expandedonly_5","CD40Lflu1_5")											#identify the samples (input file columns) of interest
overallclusterlist=findCETs(focallist,clusterlistfilename=clusterlistfilename)

###GENERATE SHORTLIST TABLE###
shortlisttable=makecombinedshortlisttable(clusterlistfilename=clusterlistfilename,
focalsamples=focallist,
pathtorawcountsfile=pathtocountsfile,
maxthreshold=50,
alpha=0.01,
chain="A")