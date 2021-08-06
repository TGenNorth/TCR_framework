The CETanalysis.R script uses the loadcountsandfreqs function to read a TCR counts file (formatted like the examples shown in the countsTables folder), and then uses the findCETs function to generate an overall list of CETs for the sample columns specified in the focallist (this list is output to workingFiles/clusterlistfilename). It then uses the makecombinedshortlisttable function to generate a shortlist of significant sets for the specified chain up to the given alpha value (this list is saved as the shortlisttable object).

Intermediate working files will appear in the workingFiles directory.

Running CETanalysis.R with default parameters will re-generate the clusters shown in Figure 3 of the Ceglia et al manuscript.

CETanalysis.R has the following dependencies (located in the functions folder):

FILE						TYPE
annotatenamechain.R				R function
CETsigtest.R					R function
extractfocal.R					R function
extracttopclonotypes.R				R function
findCETs.R					R function
findclusters.R					R function
getseqlistav.R					R function
loadcountsandfreqs.R				R function
makecombinedshortlisttable.R			R function
sigtestacrossthresholds.R			R function
TCRdist.R					R function
validateCDRs.R					R function

randomtrials/100_randomtrials_chain=A.Rdat	R data file
randomtrials/100_randomtrials_chain=B.Rdat	R data file
TCRdist/CDR1s2s2.5s.csv				csv file
TCRdist/worddist				executable file compiled from C code