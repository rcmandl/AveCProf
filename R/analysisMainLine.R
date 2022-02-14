#
# Final Analysis Main Line
#
# Note that this script is a snippet from the original script and serves as an example and has not been tested (work in progress...)
#
# René Mandl
#

WORKINGDIRECTORY <- 'path to your working directory'
RFUNCTIONSDIRECTORY<-'path to the R-functions'
SUBJECT<-'subject name for particular subject'
STUDYDIR<-'study directory'

setwd(WORKINGDIRECTORY)

# Load the optimizeProfile function
source(paste0(RFUNCTIONSDIRECTORY,'/optimizeProfiles.R'))


##
## functions
##

computeHistBoundaries<-function( inData, stdFraction ){
  # return the lower and upper value based on de stdFraction
    theMean<-mean(inData)
    theVar<-var(inData,na.rm=TRUE)
    
    minVal <- theMean - ( stdFraction * (sqrt(theVar)) )
    maxVal<-  theMean + ( stdFraction * (sqrt(theVar)) )
  list(minVal, maxVal)
}

selectProfiles<-function(inData, minCurv, maxCurv, minLength, maxLength) {
  profiles<-inData[[1]]
  theLength<-inData[[2]]
  theCurvature<-inData[[3]]
  theResult<-profiles[ ( theLength > minLength ) & ( theLength < maxLength ) & ( theCurvature > minCurv )  & ( theCurvature < maxCurv ),]
  theResult
}


#####################################################################
### MAIN LINE
#####################################################################

WM_boundary<-100
PIAL_boundary<-200
offset<-30
# curvature threshold in terms of STD of histogram
FLAT_CURV_threshold<- 1.0

# thickness treshold in therms of STD of histogram
FLAT_THICK_threshold<- 0.5

##################################################
### for right hemisphere
##################################################

theDir<-paste0(STUDYDIR,'/',SUBJECT,'/right_hemi/data')
setwd(theDir)


# now for each ROI (e.g. BA1) perform the following:

# BA1
theFile<-dir()[grep("BA1",dir())]
load(theFile)
R_BA1<-list(theResult,thickness,curvature)


len<-computeHistBoundaries(R_BA1[[2]], FLAT_THICK_threshold )
curv<-computeHistBoundaries(R_BA1[[3]], FLAT_CURV_threshold )
RBA1<-selectProfiles(R_BA1,curv[[1]], curv[[2]],len[[1]],len[[2]])
RBA1<-RBA1[,(WM_boundary-offset):(PIAL_boundary+offset)]

##
## Optimize profiles
##

res<-optimizeProfiles(RBA1,typeOfData="curvature")
RBA1Optimized<-res[[1]]



##################################################
### for right hemisphere
##################################################

theDir<-paste0(STUDYDIR,'/',SUBJECT,'/left_hemi/data')
setwd(theDir)


# BA1
theFile<-dir()[grep("BA1",dir())]
load(theFile)
L_BA1<-list(theResult,thickness,curvature)

len<-computeHistBoundaries(L_BA1[[2]], FLAT_THICK_threshold )
curv<-computeHistBoundaries(L_BA1[[3]], FLAT_CURV_threshold )
LBA1<-selectProfiles(L_BA1,curv[[1]], curv[[2]],len[[1]],len[[2]])
LBA1<-LBA1[,(WM_boundary-offset):(PIAL_boundary+offset)]


res<-optimizeProfiles(LBA1,typeOfData="curvature")
LBA1Optimized<-res[[1]]

# Now collect all variabels with 'Optimzed' in their name and store them to disk in one file.
fileName<-sprintf("%s_optimized_profiles",SUBJECT)
save(list=ls(pattern="Optimized"),file=fileName)


######################################################
######################################################
######################################################
