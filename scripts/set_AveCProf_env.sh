#!/bin/bash

# Test if the environment variable is set so that we can find all the programs needed
if [ ! -v AVECPROF ]
then
echo "You have to set the AVECPROF environment variable first"
exit 2
fi

# binaries
export DECONVOLUTION=${AVECPROF}/C++/deconvolution/itkDeconvolution
export SAMPLECORTEX=${AVECPROF}/C++/profilesampling/samplecortex

# scripts
export PERFORMSURFACESAMPLING=${AVECPROF}/scripts/perform_surface_sampling_in_patches.sh
export CREATESURFACEMEASUREMENTS=${AVECPROF}/scripts/create_subject_surface_measurements.sh

