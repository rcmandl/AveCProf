#
# computeAlignedProfiles
#
# This function computes the aligned profiles for one region.Simply averaging these results in the required Average profile. 
# There are a number of optional parameters (e.g. wmBoundary). 
# If not specified the information is collected from the info set in set_AveCProf_env.sh
#
# the result file contains a variable named optimizedProfiles [list structure] where te first element (optimizedProfiles[[1]])
# contains the aligned profiles. (the other variables contain info on the chosen reference profile and the applied warping parameters)
#
# Rene Mandl Version 0.8
#

computeAlignedProfiles<-function(inFileName, outFileName, wmBoundary=NULL, pialBoundary=NULL, offset=NULL, flatCurveThreshold=NULL, flatThickThreshold=NULL ) {
  
  OPTIMIZEPROFILES=Sys.getenv("AVECPROF_R_OPTIMIZEPROFILES")
  COMPUTEHISTBOUNDARIES=Sys.getenv("AVECPROF_R_COMPUTEHISTBOUNDARIES")
  
  # Load the optimizeProfile function
  source(OPTIMIZEPROFILES)
  # load the computeHistBoundaries function
  source(COMPUTEHISTBOUNDARIES)
  
  if (is.null(wmBoundary)) {
    wmBoundary=Sys.getenv("AVECPROF_WM_BOUNDARY")
  }
  
  if (is.null(pialBoundary)) {
    pialBoundary=Sys.getenv("AVECPROF_PIAL_BOUNDARY")
  }
  
  if (is.null(offset)) {
    offset=Sys.getenv("AAVECPROF_OFFSET")
  }
  
  if (is.null(flatCurveThreshold)) {
    flatCurveThreshold=Sys.getenv("AVECPROF_FLAT_CURVE_THRESHOLD")
  }
  
  if (is.null(flatThickThreshold)) {
    flatThickThreshold=Sys.getenv("AVECPROF_FLAT_THICK_THRESHOLD")
  }
  
   
  # small helper function  
  selectProfiles<-function(inData, minCurv, maxCurv, minLength, maxLength) {
    profiles<-inData[[1]]
    theLength<-inData[[2]]
    theCurvature<-inData[[3]]
    theResult<-profiles[ ( theLength > minLength ) & ( theLength < maxLength ) & ( theCurvature > minCurv )  & ( theCurvature < maxCurv ),]
    theResult
  }
  
  #
  # MAINLINE
  #
  
  ## the input file (result from prepRawContrastFiles.R and with the .RDAT extension) contains 3 variabels. theResult (the set of profiles), one column thickness and one column curvature.
  load(inFileName)
  cortexRegionData<-list(theResult,thickness,curvature)
  
  len<-computeHistBoundaries(cortexRegionData[[2]], FLAT_THICK_threshold )
  curv<-computeHistBoundaries(cortexRegionData[[3]], FLAT_CURV_threshold )
  selectedProfiles<-selectProfiles(cortexRegionData,curv[[1]],curv[[2]],len[[1]],len[[2]])
  
  croppedSelectedProfiles<-selectedProfiles[,(WM_boundary-offset):(PIAL_boundary+offset)]
  optimizedProfiles<-optimizeProfiles(croppedSelectedProfiles,typeOfData="curvature")
  
  save(optimizedProfiles,file=outFileName)

}
