#!/bin/bash

# Test if the environment variable is set so that we can find all the programs needed
if [ -z $AVECPROF ]
then
echo "You have to set the AVECPROF environment variable first"
exit 2
fi

# on which operating system are we running?
UNAME=`uname`
AVECPROF_OS="LINUX"
if [ "$UNAME"  = "Darwin" ]
then
AVECPROF_OS="DARWIN"
fi
export AVECPROF_OS


# binaries
export AVECPROF_BIN_DECONVOLUTION="${AVECPROF}/C++/deconvolution/itkDeconvolution"
export AVECPROF_BIN_SAMPLECORTEX="${AVECPROF}/C++/profilesampling/samplecortex"

# bash scripts
export AVECPROF_SH_PERFORMSURFACESAMPLING="${AVECPROF}/scripts/perform_surface_sampling_in_patches.sh"
export AVECPROF_SH_CREATESURFACEMEASUREMENTS="${AVECPROF}/scripts/create_subject_surface_measurements.sh"

# R scripts
export AVECPROF_R_OPTIMIZEPROFILES="${AVECPROF}/R/optimizeProfiles.R"
export AVECPROF_R_PREPRAWCORTEXFILES="${AVECPROF}/R/prepRawCortexFles.R"
export AVECPROF_R_ANALYSISMAINLINE="${AVECPROF}/R/analysisMainLine.R"
#
# script variable default settings
# you can modify these values if needed
#

export AVECPROF_SAMPLECORTEX_EXTENTFACTOR=3
export AVECPROF_SAMPLECORTEX_NROFSTEPS=300
export AVECPROF_WM_BOUNDARY=100
export AVECPROF_PIAL_BOUNDARY=200
export AVECPROF_OFFSET=30
export AVECPROF_FLAT_THICK_THRESHOLD=0.50
export AVECPROF_FLAT_CURVE_THRESHOLD=1.0
