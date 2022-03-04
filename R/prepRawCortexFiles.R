#
# prepRawCortexFiles.R
#
# This function reads in a layer file (.layer extension) generated using samplecortex, removes incomplete profiles and standardizes variance.
# The result is an .RDAT file
#
# René Mandl V 0.8
#

doPrepCortex<-function(subject, studyDir, hemisphere, contrast, area) {
	require(kernlab)
	
	CURRDIR=getwd()
	NEWDIR=paste(studyDir,"/",subject,"/layers/samples",sep="")
	setwd(NEWDIR)
       
	part<-sprintf("_.*%s.layer", contrast)
	pat=paste("^",hemisphere,".",area,part,sep="")

	Name=dir(pattern=pat)

	# Because we are working with an extension of the cortex in both directions
  # it may well be that we are sometimes sample outside the volume.
	# these values are marked as a '*'
	# First we have to remove them and replace them by NA

	M_full=as.matrix(read.table(Name, na.strings=c("*")))
		   		   
	#	
	# you can skip the first 7 columns and the last 2 columns	as they contain nrOfPoints (1), start-stop XYZ-coordinates (6), 
	# thickness (1) and curvature(1).
	# 
	   
	nCols=ncol(M_full)
	nRows=nrow(M_full)
	   
	thickness=M_full[,(nCols - 1)]
	curvature=M_full[,nCols]

	M=M_full[,8:(nCols-2)]
	   
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
