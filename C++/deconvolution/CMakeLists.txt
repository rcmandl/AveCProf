# This program uses the 2 ascii files (e.g. white and pial from freesurfer)
# representing two surfaces and samples (assuming correspondence between points of
# files) the cortex values at given thickness points from the provided minc file 
#
# René Mandl version 1.0
#

cmake_minimum_required(VERSION 3.14)

PROJECT(itkLandweberDeconvolution)
OPTION(DEBUG ON)
SET(ITK_DIR /Users/renemandl/Software/ITK/ITK)
FIND_PACKAGE(ITK)
IF (ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR
	   "Cannot build withoud ITK, please set ITK_DIR.")
ENDIF(ITK_FOUND)
ADD_EXECUTABLE(itkLandweberDeconvolution itkLandweberDeconvolution.cxx)

TARGET_LINK_LIBRARIES(itkLandweberDeconvolution ${ITK_LIBRARIES})
