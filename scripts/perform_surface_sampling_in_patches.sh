#!/bin/bash
#
# This script computes the cortex samples from the measurements computed in an earlier stage with create_subject_surface_measurements_II
# Rene Mandl
#
#

# set -x

if [ $# -ne 4 ]
then
echo "wrong number of parameters"
echo ""
echo "usage: perform_surface_sampling_in_patches <subject's freesurfer directory> <contrast.mnc> <contrast_to_anat.xfm> <patchsize>"
echo ""
echo "Computes ascii files of cortical layer measurements that can be analysed with R"
echo "The difference between this script and perform_surface_sampling is that here you generate multiple files per BA region of a given size"
echo "instead of one big file"
echo 
echo "We also include curvature information (avg curvature) which is stored in the column after the second coordinate triplet."
echo
echo "NOTE: you have to run create_subject_surface_measurements first."
echo ""
echo ""
exit 1
fi

CURRENTDIR=`pwd`

CONTRAST=`readlink -f $2`
TRANSFORMATION=`readlink -f $3`
PATCHSIZE=$4

if [ ! -e ${CONTRAST} ]
then
echo "Can't find contrast file ${CONTRAST}"
exit 2
fi

if [ ! -e ${TRANSFORMATION} ]
then
echo "Can't find transformation file ${TRANSFORMATION}" 
exit 2 
fi


cd $1
if [ ! -e layers ]
then
echo "Can't find layers directory in $1"
echo " You have to run create_subject_surface_measurements_II first."
echo ""
exit 2
fi

cd layers

MINCSAMPLECORTEX=/Users/mandl/Sources/myProgsmac/myDTImac/sDTI_pipeline_mac/UMCU_minc_tools_1.0/mincsamplecortex/mincsamplecortex
LAYER22COLUMNS=/Users/mandl/Sources/myProgsmac/myDTImac/sDTI_pipeline_mac/UMCU_minc_tools_1.0/scripts/layer_to_2_columns.awk

mkdir samples

xfmconcat ${TRANSFORMATION} orig_anat_to_target.xfm ./samples/contrast_to_surface.xfm
xfminvert ./samples/contrast_to_surface.xfm ./samples/surface_to_contrast.xfm

BC=`basename $CONTRAST .gz`
BC=`basename $BC .mnc`

for x in lh.BA1.pial lh.BA4a.pial rh.BA44.pial rh.V1.pial lh.BA2.pial lh.BA4p.pial lh.entorhinal_exvivo.pial rh.BA45.pial rh.V2.pial lh.BA3a.pial lh.BA6.pial rh.BA1.pial rh.BA4a.pial lh.BA3b.pial lh.MT.pial rh.BA2.pial rh.BA4p.pial rh.entorhinal_exvivo.pial lh.BA44.pial lh.V1.pial rh.BA3a.pial rh.BA6.pial lh.BA45.pial lh.V2.pial rh.BA3b.pial rh.MT.pial
do
BASE=`basename $x .pial`
#${MINCSAMPLECORTEX} -forKriger -scalarFile ${BASE}.curvature ${BASE}.white ${BASE}.pial ${CONTRAST} ./samples/surface_to_contrast.xfm 1 1 > ./samples/${BASE}_${BC}.layer
${MINCSAMPLECORTEX} -linearInterpolate -extentFactor 3 -scalarFile ${BASE}.curvature ${BASE}.white ${BASE}.pial ${CONTRAST} ./samples/surface_to_contrast.xfm 300 0 > ./samples/${BASE}_${BC}.layer.tmp
# with larger files (e.g. higher resolution scans) the program with start and end with textmessages regarding the progress (e.g. "Reading Volume...." and "Outputting Volume.....")
# To remove these unwanted messages from ou results we use the follwing awk line
cat ./samples/${BASE}_${BC}.layer.tmp | awk '//{if ($1 != "Reading" && $1 != "Outputting") {print $0} }' > ./samples/${BASE}_${BC}.layer
rm ./samples/${BASE}_${BC}.layer.tmp

mincmath -const 1 -mult flagVolume.mnc tmp1.mnc
rm flagVolume.mnc
mv tmp1.mnc ./samples/${BASE}_${BC}.mnc
gzip ./samples/${BASE}_${BC}.mnc

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
# NEWNAME2C=`printf "./samples/%s_%s_%05d.layer2c" ${BASE} ${BC} ${NROFPATCHES}`
# set interval to 2 cause we already splitted the file in patches so there is no need to reduce the file size any further (for now)
# awk -v INTERVAL=2 -f ${LAYER22COLUMNS} < ${NEWNAME} > ${NEWNAME2C}
done


done

cd ${CURRENTDIR}
