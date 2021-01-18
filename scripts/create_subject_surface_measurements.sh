#!/bin/zsh
#
# This program creates the surface files for the pial an dn white matter surface that serve as input for the program samplecortex.
#
# René Mandl
#
# Version 1.0
#


#set -x

NPARAM=$#

if (( NPARAM < 1 || NPARAM > 2 ))
then
echo
echo "wrong number of parameters."
echo
echo  "usage: create_subject_surface_measurements <subject's freesurfer directory> [path to label files]"
echo
echo ""   
echo "Creates surface files (pial ad white) for given label files, that serve as input for samplecortex"
echo "These surfaces are in native space (e.g. brain.mgz)"
echo
echo "If no path to label files is provided then te default label files are used (found in label directory of the subject's directory)"
echo "You can create label files from annotation files using the program mri_annotation2label (part of freesurfer)"
echo "In this way you can create labels for complete atlases. Note that you have to create these labels before you call this script."
echo
echo "Requires that both Rscript (part of R) and Freesurfer are installed." 
echo
echo
exit 1
fi


# Check if we can find the required software

if [ ! -v AVECPROF ]
then
echo "You have to set the AVECPROF environment variable first"
exit 2
fi

if [ ! -v FREESURFER_HOME ]
then
echo "We make use of various Freesurfer commands, so please install Freesurfer first (environment variable FREESURFER_HOME not set)"
exit 2
fi

# check if Rscript is installed.

if ! command -v Rscript > /dev/null 2>&1
then
echo "Cannot find Rscript, please install that program first."
exit 2
fi


CURRENTDIR=`pwd`
cd $1
mkdir layers
cd layers

if (( 2 == NPARAM ))
then
LABELPATH=$2
else
LABELPATH="../label"
fi

#
# temporary copy the required cortex info to this directory
#

cp ../surf/lh.pial .
cp ../surf/lh.white .
cp ../surf/rh.pial .
cp ../surf/rh.white .
cp ../surf/lh.avg_curv .
cp ../surf/rh.avg_curv .

# get curvature information in ascii format add to coordiates (in native coordinates)

mris_convert --to-scanner -c lh.avg_curv lh.pial lh.pial.avg_curv.asc
mris_convert --to-scanner -c lh.avg_curv lh.white lh.white.avg_curv.asc
mris_convert --to-scanner -c rh.avg_curv rh.pial rh.pial.avg_curv.asc
mris_convert --to-scanner -c rh.avg_curv rh.white rh.white.avg_curv.asc

rm lh.pial lh.white rh.pial rh.white



##
## For left hemipshere
##

#
# Now construct the small R script that is used to select the points (for both the pial and the wm surface
# that are part of the particular label.

echo "# R script to select points from inner and out surface given their index in a file." > selectPoints_l.R
echo "# Assuming that the first 2 lines of the label and surface files are removed." >> selectPoints_l.R
echo "indices <- read.table(\"tmpLabel_l\")" >> selectPoints_l.R
echo "pial <- read.table(\"lh.pial.avg_curv.asc\");" >> selectPoints_l.R
echo "wm <- read.table(\"lh.white.avg_curv.asc\");" >> selectPoints_l.R
echo "i<-indices[,1]" >> selectPoints_l.R
echo "i=i+1; # indices start at 0 while R starts at 1" >> selectPoints_l.R
echo "pialSelected<-pial[i,2:4]" >> selectPoints_l.R
echo "curvatureSelected<-pial[i,5]" >> selectPoints_l.R
echo "write.table(pialSelected, file=\"PIALOUTFILE_l\",col.names=FALSE, row.names=FALSE)"  >> selectPoints_l.R
echo "write.table(curvatureSelected, file=\"CURVATUREOUTFILE_l\",col.names=FALSE, row.names=FALSE)" >> selectPoints_l.R
echo "wmSelected<-wm[i,2:4]" >> selectPoints_l.R
echo "write.table(wmSelected, file=\"WHITEOUTFILE_l\",col.names=FALSE, row.names=FALSE)" >> selectPoints_l.R

# Now for each label file in the label directory do:

for x in `ls ${LABELPATH}/lh.*.label`
do
   # first line removes the two first lines because they contain header info,
   awk '//{if (NR>2){print $0}}' < ${x} > tmpLabel_l
   Rscript selectPoints_l.R
   BASE=`basename $x .label`
   mv PIALOUTFILE_l ${BASE}.pial
   mv WHITEOUTFILE_l ${BASE}.white
   mv CURVATUREOUTFILE_l ${BASE}.curvature
   rm tmpLabel_l
done

rm selectPoints_l.R



##
## for the right hemisphere
##

#
# Now construct the small R script that is used to select the points (for both the pial and the wm surface
# that are part of the particular label.

echo "# R script to select points from inner and out surface given their index in a file." > selectPoints_r.R
echo "# Assuming that the first 2 lines of the label and surface files are removed." >> selectPoints_r.R
echo "indices <- read.table(\"tmpLabel_r\")" >> selectPoints_r.R
echo "pial <- read.table(\"rh.pial.avg_curv.asc\");" >> selectPoints_r.R
echo "wm <- read.table(\"rh.white.avg_curv.asc\");" >> selectPoints_r.R
echo "i<-indices[,1]" >> selectPoints_r.R
echo "i=i+1; # indices start at 0 while R starts at 1" >> selectPoints_r.R
echo "pialSelected<-pial[i,2:4]" >> selectPoints_r.R
echo "curvatureSelected<-pial[i,5]" >> selectPoints_r.R
echo "write.table(pialSelected, file=\"PIALOUTFILE_r\",col.names=FALSE, row.names=FALSE)"  >> selectPoints_r.R
echo "write.table(curvatureSelected, file=\"CURVATUREOUTFILE_r\",col.names=FALSE, row.names=FALSE)" >> selectPoints_r.R
echo "wmSelected<-wm[i,2:4]" >> selectPoints_r.R
echo "write.table(wmSelected, file=\"WHITEOUTFILE_r\",col.names=FALSE, row.names=FALSE)" >> selectPoints_r.R

# Now for each label file in the label directory do:

for x in `ls ${LABELPATH}/rh.*.label`
do
   # first line removes the two first lines because they contain header info,
   awk '//{if (NR>2){print $0}}' < ${x} > tmpLabel_r
   Rscript selectPoints_r.R
   BASE=`basename $x .label`
   mv PIALOUTFILE_r ${BASE}.pial
   mv WHITEOUTFILE_r ${BASE}.white
   mv CURVATUREOUTFILE_r ${BASE}.curvature
   rm tmpLabel_r
done

rm selectPoints_r.R lh.pial.avg_curv.asc lh.white.avg_curv.asc rh.pial.avg_curv.asc rh.white.avg_curv.asc

# and return back to original dir 
cd $CURRENTDIR
