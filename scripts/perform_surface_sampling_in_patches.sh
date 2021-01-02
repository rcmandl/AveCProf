#!/bin/bash
#
# This script computes the cortex samples from the measurements computed in an earlier stage with create_subject_surface_measurements
# Rene Mandl
#
#

# set -x

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

CONTRAST=`readlink -f $2`
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
LAYER22COLUMNS=${AVECPROF}scripts/layer_to_2_columns.awk

mkdir samples

BC=`basename $CONTRAST .gz`
BC=`basename $BC .nii`

for x in `ls *.pial`
  do
  BASE=`basename $x .pial`

  ${SAMPLECORTEX} --linearInterpolation --extentFactor 3 --nrOfSteps 300 --scalarAddOn ${BASE}.curvature --flagVolume ${BASE}.pial ${BASE}.white ${CONTRAST}  > ./samples/${BASE}_${BC}.layer

  split -a 5 -l ${PATCHSIZE} ./samples/${BASE}_${BC}.layer ./samples/${BASE}_${BC}_xxxxx_

  #
  # Now, what is done below is only of interest if you would use the forKriger option.
  #

  NROFPATCHES=0
  for y in `ls ./samples/${BASE}_${BC}_xxxxx_*`
    do
    ((NROFPATCHES++))
    NEWNAME=`printf "./samples/%s_%s_%05d.layer" ${BASE} ${BC} ${NROFPATCHES}`
    mv $y ${NEWNAME}
  done
done

cd ${CURRENTDIR}
