#
# computeHistBoundaries
#
# René Mandl Version 1.0
#

computeHistBoundaries<-function( inData, stdFraction ){
  # return the lower and upper value based on de stdFraction
    theMean<-mean(inData)
    theVar<-var(inData,na.rm=TRUE)
    
    minVal <- theMean - ( stdFraction * (sqrt(theVar)) )
    maxVal<-  theMean + ( stdFraction * (sqrt(theVar)) )
  list(minVal, maxVal)
}
