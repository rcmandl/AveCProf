##
## optimizeProfiles
##
## Version 1.3
##
##
## In version 1.2 I added te possibility to provide a second set of profiles. Optimisation is still computed
## on the first set but the optimisation is also applied on the second set (if provided)
## This allows to compute the optimisation on the deconvoluted set and apply it on the original profiles
## Version 1.3 also includes th assymetric least squares (ASL) option to remove the baseline (from the ptw package) as well as the spline implementation
## For now we use the spline implementaton. because the ASL may require more tweaking to obtain acceptable results.
##
## Note that NA's in the profiles are not replaced by zeroes. 
##
##  Rene Mandl


Curvature<-function(theLine,delta) {
	# compute the curvature. This is done by computing the distance between
	# the  point P on the function theLine and the line segment defined by (p-delta,p+delta)
	
	dist2Line<-function(px,py,x1,y1,x2,y2) {
		# compute distance between line segment and point
		len=function(x1,y1,x2,y2) sqrt((x2-x1)^2+(y2-y1)^2)
		segmentLength=len(x1,y1,x2,y2)
		u=0
		v=0
		result=NULL
	
		s=( ((( px - x1 ) * ( x2-x1 )) + (( py-y1 ) * ( y2 - y1 ))) ) / ( segmentLength * segmentLength ) 
		if ( s < 0.000000001 || s > 1.0) {
			a=len(px,py,x1,y1)
			b=len(px,py,x2,y2)
			if (a>b) result=b
			else result=a
		} else {
			a = x1 + s * ( x2 -x1 )
			b = y1 + s * ( y2 -y1 )
			result = len(px,py,a,b)
		}
   	result
	}

    # main line Curvature
	l<-length(as.vector(theLine))
	
	result<-theLine*0
	for ( p in c((1+delta):(l-delta)) ) {
		result[p]<-dist2Line(p,theLine[p],p-delta,theLine[p-delta],p+delta,theLine[p+delta])
		# now test if curvature is negative
		if ( ( ( theLine[p - delta] + theLine[p + delta] ) / 2 ) > theLine[p] ) { result[p] = result[p] * -1 }
	}
	result
}



optimizeProfiles<-function(theProfiles,typeOfData='normal',secondSetProfiles=NULL, bestRefProfile=NULL) {



getSticks<-function(profile,delta=9) {
# require(ptw)		
#
# main line getSticks
#
profile<-as.numeric(profile)
result<-profile * 0
# compute curvature, the second argument describes the size of + and - delta
curvature<-Curvature(profile,delta)
# get spline representation to get rid of small local minima and maxima
curveSpline<-smooth.spline(as.numeric(curvature),df=15)
# get points representation of this smoothed curvature
curveSplinePoints<-predict(curveSpline,1:length(profile))$y
# now create stick representation of profile

    for (j in 2:(length(curveSplinePoints)-1)) {
            if ( abs(curveSplinePoints[j-1]) < abs(curveSplinePoints[j]) & abs(curveSplinePoints[j+1]) < abs(curveSplinePoints[j]) ) {
                   result[j]=profile[j]
            }
    }
    
result
}


	###################################
	### MAIN LINE OPTIMIZE PROFILES ###
	###################################


	require(ptw)

    delta<-9

    if (typeOfData == 'sticks') {
    	theProfiles<-as.matrix(theProfiles)
    	bestRepresentation<-bestref(theProfiles)
    	if ( is.null(bestRefProfile) ) {
    	  theModel<-theProfiles[bestRepresentation$best.ref,]
    	} else { theModel<-bestRefProfile }
    	modelSticks<-getSticks(theModel,delta)

    	profileLength<-ncol(theProfiles)
    	nrOfProfiles<-nrow(theProfiles)
    
        warpCoefs<-matrix(nrow=nrOfProfiles,ncol=3)
        
    	#
    	# create sticks for all profiles
    	#
    
    	allSticks<-theProfiles*0
    
    	for (i in 1:nrOfProfiles) {
    	    allSticks[i,]<-getSticks(theProfiles[i,],delta)
    	    # remove info outside delta+1 as it is not reliable
    	    allSticks[i,1:(delta+1)] <- 0
    	    allSticks[i, (profileLength - delta - 1): profileLength] <- 0
    	}

    	#
   		# now perform final optimisation
    	#
    
    	maxDisplacement=20
    
    	trgtm = cbind( which(modelSticks != 0, arr.ind=TRUE), modelSticks[ which(modelSticks != 0, arr.ind=TRUE) ])
    	colnames(trgtm)=c("rt","I")

    	optimizedProfiles<-theProfiles * 0

    	for (i in 1:nrOfProfiles) {
       		src=allSticks[i,]
        	srcm = cbind( which(src != 0, arr.ind=TRUE), src[ which(src != 0, arr.ind=TRUE) ])
        	colnames(srcm)=c("rt","I")
        
        	stptwResult<-stptw(list(trgtm),list(srcm),trwdth=maxDisplacement)
        
        	new.times <- warp.time(1:profileLength,stptwResult$warp.coef)
        	#result <- approx(new.times, theProfiles[i,],1:profileLength)
        	endlocation<-length(theProfiles[i,])
        	result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
        	optimizedProfiles[i,]<-result$y
        	warpCoefs[i,(1:ncol(stptwResult$warp.coef))]<-stptwResult$warp.coef
    	}
    
    	# commented out 0 substition, now missing values stay NA
    	# optimizedProfiles[is.na(optimizedProfiles)]=0
    	return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs))
    } else if ( typeOfData == 'curvature') {
    	theProfiles<-as.matrix(theProfiles)
    	bestRepresentation<-bestref(theProfiles)
    	# bestRepresentation<-bestref(Curvature(theProfiles, delta))
    	if ( is.null(bestRefProfile) ) {
    	  theModel<-theProfiles[bestRepresentation$best.ref,]
    	} else { theModel<-bestRefProfile }
    	
    	modelCurvature<-Curvature(theModel,delta)

    	profileLength<-ncol(theProfiles)
    	nrOfProfiles<-nrow(theProfiles)
    	
    	warpCoefs<-matrix(nrow=nrOfProfiles,ncol=3)
    
    	#
    	# create curvatures for all profiles
    	#
    
    	allCurvatures<-theProfiles*0
    
    	for (i in 1:nrOfProfiles) {
    	    allCurvatures[i,]<-Curvature(theProfiles[i,],delta)
    	}

    	#
   		# now perform final optimisation
    	#
    
    	maxDisplacement=20
    
    	optimizedProfiles<-theProfiles * 0

    	if ( is.null(secondSetProfiles) ) {
    	    for (i in 1:nrOfProfiles) {
       		  src=allCurvatures[i,]
        
        	  Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
        	
        	  new.times <- warp.time(1:profileLength,c(Result$warp.coef))
        	  #result <- approx(new.times, theProfiles[i,],1:profileLength)
        	  endlocation<-length(theProfiles[i,])
        	  result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
        	  optimizedProfiles[i,]<-result$y  
          	warpCoefs[i,1:(ncol(Result$warp.coef))]<-Result$warp.coef      	
    	    }
    
    	  # commented out 0 substition, now missing values stay NA
    	  # optimizedProfiles[is.na(optimizedProfiles)]=0
    	} else {
    	
    	optimizedSecondSet=secondSetProfiles*0
    	
    	for (i in 1:nrOfProfiles) {
    	  src=allCurvatures[i,]
    	  
    	  Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
    	  
    	  new.times <- warp.time(1:profileLength,c(Result$warp.coef))
    	  #result <- approx(new.times, theProfiles[i,],1:profileLength)
    	  endlocation<-length(theProfiles[i,])
    	  result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
    	  optimizedProfiles[i,]<-result$y  
    	  warpCoefs[i,1:(ncol(Result$warp.coef))]<-Result$warp.coef
    	  # apply warping to second set as well
    	  #result <- approx(new.times, secondSetProfiles[i,],1:profileLength)
    	  endlocation<-length(theProfiles[i,])
    	  result <- approx(new.times, secondSetProfiles[i,], 1:profileLength, yleft=secondSetProfiles[i,1], yright=secondSetProfiles[i,endlocation])
    	  optimizedSecondSet[i,]<-result$y
    	}
    	
    # optimizedProfiles[is.na(optimizedProfiles)]=0
    # optimizedSecondSet[is.na(optimizedSecondSet)]=0
    	}
    	if ( is.null(secondSetProfiles) ) {
    	return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs))
    	} else {
    	  return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs,optimizedSecondSet))
    	}
    	
    } else if ( typeOfData == 'normal') {
    	theProfiles<-as.matrix(theProfiles)
    	bestRepresentation<-bestref(theProfiles)
    	if ( is.null(bestRefProfile) ) {
    	  theModel<-theProfiles[bestRepresentation$best.ref,]
    	} else { theModel<-bestRefProfile }

    	profileLength<-ncol(theProfiles)
    	nrOfProfiles<-nrow(theProfiles)
 
        warpCoefs<-matrix(nrow=nrOfProfiles,ncol=3)   
 
    	#
   		# now perform final optimisation
    	#
    
    	maxDisplacement=20
    
    	optimizedProfiles<-theProfiles * 0

    	for (i in 1:nrOfProfiles) {
       		src=theProfiles[i,]
        
        	Result<-ptw(theModel,src,init.coef=c(0,1,0),trwdth=maxDisplacement)
        	        	
        	new.times <- warp.time(1:profileLength,c(Result$warp.coef))
        	#result <- approx(new.times, theProfiles[i,],1:profileLength)
        	endlocation<-length(theProfiles[i,])
        	result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
        	optimizedProfiles[i,]<-result$y 
        	warpCoefs[i,]<-Result$warp.coef 
    	}
    
    	# commented out 0 substition, now missing values stay NA
    	# optimizedProfiles[is.na(optimizedProfiles)]=0
    	
    	return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs))
    	
    	
    } else if (typeOfData == 'ALS' ) {
      # asymmetric least squares correction of data
      
      theProfiles<-as.matrix(theProfiles)
      bestRepresentation<-bestref(theProfiles)
      if ( is.null(bestRefProfile) ) {
        theModel<-theProfiles[bestRepresentation$best.ref,]
      } else { theModel<-bestRefProfile }
      modelCurvature<-baseline.corr(theModel)
      
      profileLength<-ncol(theProfiles)
      nrOfProfiles<-nrow(theProfiles)
      
      warpCoefs<-matrix(nrow=nrOfProfiles,ncol=3)
      
      #
      # create curvatures for all profiles
      #
      
      allCurvatures<-theProfiles*0
      
      allCurvatures<-baseline.corr(theProfiles)
      
      #
      # now perform final optimisation
      #
      
      maxDisplacement=20
      
      optimizedProfiles<-theProfiles * 0
      
      if ( is.null(secondSetProfiles) ) {
        for (i in 1:nrOfProfiles) {
          src=allCurvatures[i,]
          
          Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
          
          new.times <- warp.time(1:profileLength,c(Result$warp.coef))
          #result <- approx(new.times, theProfiles[i,],1:profileLength)
          endlocation<-length(theProfiles[i,])
          result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
          optimizedProfiles[i,]<-result$y  
          warpCoefs[i,]<-Result$warp.coef      	
        }
        
        # commented out 0 substition, now missing values stay NA
        # optimizedProfiles[is.na(optimizedProfiles)]=0
      } else {
        
        optimizedSecondSet=secondSetProfiles*0
        
        for (i in 1:nrOfProfiles) {
          src=allCurvatures[i,]
          
          Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
          
          new.times <- warp.time(1:profileLength,c(Result$warp.coef))
          #result <- approx(new.times, theProfiles[i,],1:profileLength)
          endlocation<-length(theProfiles[i,])
          result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
          optimizedProfiles[i,]<-result$y  
          warpCoefs[i,]<-Result$warp.coef
          # apply warping to second set as well
          #result <- approx(new.times, secondSetProfiles[i,],1:profileLength)
          endlocation<-length(theProfiles[i,])
          result <- approx(new.times, secondSetProfiles[i,], 1:profileLength, yleft=secondSetProfiles[i,1], yright=secondSetProfiles[i,endlocation])
          optimizedSecondSet[i,]<-result$y
        }
        
       # optimizedProfiles[is.na(optimizedProfiles)]=0
       #  optimizedSecondSet[is.na(optimizedSecondSet)]=0
      }
      if ( is.null(secondSetProfiles) ) {
        return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs))
      } else {
        return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs,optimizedSecondSet))
      }
      
     }  else if (typeOfData == 'spline' ) {
        # use splines (should be the same as our method)
        
        theProfiles<-as.matrix(theProfiles)
        bestRepresentation<-bestref(theProfiles)
        if ( is.null(bestRefProfile) ) {
        theModel<-theProfiles[bestRepresentation$best.ref,]
        } else { theModel<-bestRefProfile }
        
        modelCurvature<-theModel - (smooth.spline(theModel,df=7)$y)
        
        profileLength<-ncol(theProfiles)
        nrOfProfiles<-nrow(theProfiles)
        
        warpCoefs<-matrix(nrow=nrOfProfiles,ncol=2)
        
        #
        # create curvatures for all profiles
        #
        
        allCurvatures<-theProfiles*0
        
        for (i in 1:nrOfProfiles) {
          allCurvatures[i,]<-theProfiles[i,] - (smooth.spline(theProfiles[i,],df=7)$y)
        }
        
        #
        # now perform final optimisation
        #
        
        maxDisplacement=20
        
        optimizedProfiles<-theProfiles * 0
        
        if ( is.null(secondSetProfiles) ) {
          for (i in 1:nrOfProfiles) {
            src=allCurvatures[i,]
            
            Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
            
            new.times <- warp.time(1:profileLength,c(Result$warp.coef))
            #result <- approx(new.times, theProfiles[i,],1:profileLength)
            endlocation<-length(theProfiles[i,])
            result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
            optimizedProfiles[i,]<-result$y  
            warpCoefs[i,]<-Result$warp.coef      	
          }
          
          # optimizedProfiles[is.na(optimizedProfiles)]=0
        
        
        #
        # now perform final optimisation
        #
        
        maxDisplacement=20
        
        optimizedProfiles<-theProfiles * 0
        
        if ( is.null(secondSetProfiles) ) {
          for (i in 1:nrOfProfiles) {
            src=allCurvatures[i,]
            
            Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
            
            new.times <- warp.time(1:profileLength,c(Result$warp.coef))
            result <- approx(new.times, theProfiles[i,],1:profileLength)
            
            optimizedProfiles[i,]<-result$y  
            warpCoefs[i,]<-Result$warp.coef      	
          }
          
          # optimizedProfiles[is.na(optimizedProfiles)]=0
        } else {
          
          optimizedSecondSet=secondSetProfiles*0
          
          for (i in 1:nrOfProfiles) {
            src=allCurvatures[i,]
            
            Result<-ptw(modelCurvature[delta:(profileLength -delta)],src[delta:(profileLength -delta)],trwdth=maxDisplacement,optim.crit="WCC",init.coef=c(0,1))
            
            new.times <- warp.time(1:profileLength,c(Result$warp.coef))
            #result <- approx(new.times, theProfiles[i,],1:profileLength)
            endlocation<-length(theProfiles[i,])
            result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
            optimizedProfiles[i,]<-result$y  
            warpCoefs[i,]<-Result$warp.coef
            # apply warping to second set as well
            #result <- approx(new.times, secondSetProfiles[i,],1:profileLength)
            endlocation<-length(theProfiles[i,])
            result <- approx(new.times, secondSetProfiles[i,], 1:profileLength, yleft=secondSetProfiles[i,1], yright=secondSetProfiles[i,endlocation])
            optimizedSecondSet[i,]<-result$y
          }
          
        #  optimizedProfiles[is.na(optimizedProfiles)]=0
        #  optimizedSecondSet[is.na(optimizedSecondSet)]=0
        }
        if ( is.null(secondSetProfiles) ) {
          return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs))
        } else {
          return(list(optimizedProfiles, bestRepresentation$best.ref, bestRepresentation$crit.values, warpCoefs,optimizedSecondSet))
        }
      
      } } else {
    	print(" ")
    	print(" ERROR: typeOfData should be 'normal', 'curvature','sticks', 'ALS' or 'spline'")
    	print(" ")
    }
}

##
## applyStpWarp()
##

applyStpWarp<-function (theProfiles, warpCoef) {
	
	# Allows to apply one warp to all profiles
	theProfiles<-as.matrix(theProfiles)
	profileLength<-ncol(theProfiles)
    nrOfProfiles<-nrow(theProfiles)
    	
	optimizedProfiles<-theProfiles * 0
	
	for (i in 1:nrOfProfiles) {
		new.times <- warp.time(1:profileLength,c(warpCoef))
        	#result <- approx(new.times, theProfiles[i,],1:profileLength)
        	endlocation<-length(theProfiles[i,])
        	result <- approx(new.times, theProfiles[i,], 1:profileLength, yleft=theProfiles[i,1], yright=theProfiles[i,endlocation])
        	optimizedProfiles[i,]<-result$y 
    }
    
	return(list(optimizedProfiles))
}



#
# code snippet (as an example):
#
# select the profiles for each cluster
#
# profilesCluster1<-theProfiles[clusterNumbers==1,80:220]
# profilesCluster2<-theProfiles[clusterNumbers==2,80:220]
# profilesCluster3<-theProfiles[clusterNumbers==3,80:220]
#
# optimize the profiles for each individual cluster
#
#result1<-optimizeProfiles(profilesCluster1,typeOfData="sticks")
# result2<-optimizeProfiles(profilesCluster2,typeOfData="sticks")
# result3<-optimizeProfiles(profilesCluster3],typeOfData="sticks")
#
# compute the means (or used the best representations)
#
# theMeans[1,]<-colMeans(result1[[1]])
# theMeans[2,]<-colMeans(result2[[1]])
# theMeans[3,]<-colMeans(result3[[1]])
#
# Now align the means of the clusters
#
# theMeansOptimized<-optimizeProfiles(theMeans,typeOfData="normal")
#
# finally apply the warps found for the means to each indiviual profile of the clusters
#
# finalCluster1<-applyStpWarp(result1[[1]],theMeansOptimized[[4]][1,])
# finalCluster2<-applyStpWarp(result2[[1]],theMeansOptimized[[4]][2,])
# finalCluster3<-applyStpWarp(result3[[1]],theMeansOptimized[[4]][3,])
#
#  Now finalCluster1 contains all the aligned profiles for cluster 1 etc.

