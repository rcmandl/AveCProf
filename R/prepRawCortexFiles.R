############################
# prepRawCortexFiles.R
#
# This function reads in a .layer file generated using samplecortex, removes incomplete profiles and standardizes variance.
# René Mandl V 0.8
#
# NOT TESTED YET
#
############################

#################
### doPrepCortex
################

# This function reads in the profile data for specific contrast for a single layer from a single subject.

doPrepCortex<-function(subject, studyDir, hemisphere, contrast, area) {
	require(kernlab)
	
	CURRDIR=getwd()
	NEWDIR=paste(studyDir,"/",subject,"/layers/samples",sep="")
	setwd(NEWDIR)
       
	part<-sprintf("_.*%s.layer", contrast)
	pat=paste("^",hemisphere,".",area,part,sep="")

	Name=dir(pattern=pat)

	#
	# Because we are working with an extension of the cortex in both directions
    # it may well be that we are sometimes sample outside the volume.
	# these values are marked as a '*'
	# First we have to remove them and replace them by NA

	M=as.matrix(read.table(Name, na.strings=c("*")))
		   		   
	 	   
	#	
	# you can skip the first 8 columns and the last one	(e.g. start-stop coordinates nr of points)
	# 
	   
	nCols=ncol(M)
	nRows=nrow(M)
	   
	thickness=sqrt( (M[,2] - M[,5])^2 + (M[,3]- M[,6]) ^2 + (M[,4] -M[,7]) ^2 )
	curvature=M[,8]

	M=M[,9:(nCols-1)]
	   
    # now replace all '0' by NA 
	
	M[M==0]=NA
	   
	# now remove all rows for which there is at least one contrast with all NAs found
	   
	rowM=rowMeans(M,na.rm=TRUE)
	# remove rows from matrices
	M=M[!is.nan(rowM),]
	# remove same rows from thickness and curvature   
	thickness=thickness[!is.nan(rowM)]
	curvature=curvature[!is.nan(rowM)]
       
	# now recompute nCols and nRows
	nCols=ncol(M)
	nRows=nrow(M)
	   
	rowM=rowMeans(M,na.rm=TRUE)
	rowMM=matrix(rowM,nrow=nRows,ncol=nCols)
	M=M-rowMM
	  
	# now standardise variance
	rowS=apply(M,1,sd)
	rowSM=matrix(rowS,nrow=nRows,ncol=nCols)
	M=M/rowSM
	# And to ensure that no new NA's are introduced, remove them again (if any)
	rowM=rowMeans(M,na.rm=TRUE)
	# remove rows from matrices
	M=M[!is.nan(rowM),]
	# remove same rows from thickness and curvature   
	thickness=thickness[!is.nan(rowM)]
	curvature=curvature[!is.nan(rowM)]
	   
	theResult=M	   	     	   
	    
  outName=paste(studyDir,"/layermatrices/","boundaries_",subject,"_",hemisphere,"_",area,"_",contrast,".RDAT",sep="")
 
  save(theResult, thickness, curvature, file=outName)
}



###############################################
## example how to call function
###############################################

# subject = "SUBJECT001"
# studyDir = "/HOME?JOHNDO/CORTEXLAYERS
# contrast="T1_HR_deconvoluted"
# hemisphere="lh"
# area="V1"

# doPrepCortex<-function(subject, theDir, hemisphere, contrast, area)
