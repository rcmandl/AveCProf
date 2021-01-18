#!/bin/zsh
#
# This script computes the cortex samples from the measurements computed in an earlier stage with create_subject_surface_measurements
# Rene Mandl
#
#

 set -x

if [ $# -ne 3 ]
then
echo "wrong number of parameters"
echo ""
echo "usage: perform_surface_sampling_in_patches <subject's freesurfer directory> <contrast.nii> <patchsize>"
echo ""
echo "Computes ascii files of cortical layer measurements that can be analysed with R"
echo "We also include curvature information (avg curvature) which is stored in the column after the second coordinate triplet."
echo
echo "NOTE: you have to run create_subject_surface_measurements first."
echo ""
echo ""
exit 1
fi

if [ ! -v AVECPROF ]
then
echo "You have to set the AVECPROF environment variable first"
exit 2
fi

CURRENTDIR=`pwd`

if [ "${AVECPROF_OS}" = "DARWIN" ]
then
CONTRAST=$2
else
CONTRAST=`readlink -l $2`
fi


PATCHSIZE=$3

if [ ! -e ${CONTRAST} ]
then
echo "Can't find contrast file ${CONTRAST}"
exit 3
fi

cd $1
if [ ! -e layers ]
then
echo "Can't find layers directory in $1"
echo " You have to run create_subject_surface_measurements first."
echo ""
exit 4
fi

cd layers

SAMPLECORTEX=${AVECPROF}/C++/profilesampling/samplecortex
LAYER22COLUMNS=${AVECPROF}/scripts/layer_to_2_columns.awk

mkdir samples

BC=`basename $CONTRAST .gz`
BC=`basename $BC .nii`

# Test if the environmnent variables are set; otherwise use a default value

AVECPROF_SAMPLECORTEX_EXTENTFACTOR="${AVECPROF_SAMPLECORTEX_EXTENTFACTOR:-3}"
AVECPROF_SAMPLECORTEX_NROFSTEPS="${AVECPROF_SAMPLECORTEX_NROFSTEPS:-300}"

for x in `ls *.pial`
  do
  BASE=`basename $x .pial`

  # Note the negateXY option; for some reason there appears to be a discrepancy between the world-to-scan transformation in the nifti header file and the internal ITK representation
  # (does not look like  a RAS <-> LPS issue). For now it is 'solved' by adding an option to correct for this (multiplying X and Y coordinates with -1).
  # With the scalarAddOn option we include curvature information into the file (thickness can be computed by the euclidian distance between a pial point and its corresponding wm point). 


  ${SAMPLECORTEX} --linearInterpolation --negateXY --extentFactor ${AVECPROF_SAMPLECORTEX_EXTENTFACTOR} --nrOfSteps ${AVECPROF_SAMPLECORTEX_NROFSTEPS} --scalarAddOn ${BASE}.curvature --flagVolume ${BASE}_flag.nii ${BASE}.white ${BASE}.pial ${CONTRAST}  > ./samples/${BASE}_${BC}.layer
 
  # we just create flag volumes to be sure you can overlay them on the contrast file to check if the sampling is valid 
  gzip ${BASE}_flag.nii
  mv ${BASE}_flag.nii.gz ./samples 
 
  split -a 5 -l ${PATCHSIZE} ./samples/${BASE}_${BC}.layer ./samples/${BASE}_${BC}_xxxxx_

  NROFPATCHES=0
  for y in `ls ./samples/${BASE}_${BC}_xxxxx_*`
    do
    ((NROFPATCHES++))
    NEWNAME=`printf "./samples/%s_%s_%05d.layer" ${BASE} ${BC} ${NROFPATCHES}`
    mv $y ${NEWNAME}
  done
done

gzip ./samples/*.layer

cd ${CURRENTDIR}
