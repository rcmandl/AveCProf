#!/bin/bash
#
# This program creates the surface files for the pial an dn white matter surface that serve as input for the program samplecortex.
#
# René Mandl
#
# Version 1.0
#


# set -x
if [ $# -ne 1 ]
then
echo
echo "wrong number of parameters."
echo
echo  "usage: create_subject_surface_measurements <subject's freesurfer directory>"
echo
echo ""   
echo "Creates surface files (pial ad white) that are ready to be used with samplecortex (found in layer directory)"   
echo "These surfaces are in native space (e.g. brain.mgz)"
echo ""   
exit 1
fi

CURRENTDIR=`pwd`
cd $1
mkdir layers
cd layers

#
#
#

cp ../surf/lh.pial .
cp ../surf/lh.white .
cp ../surf/rh.pial .
cp ../surf/rh.white .
cp ../surf/lh.avg_curv .
cp ../surf/rh.avg_curv .

mris_convert --cras_add --to-scanner lh.pial lh.pial.asc
mris_convert --cras_add --to-scanner lh.white lh.white.asc
mris_convert --cras_add --to-scanner rh.pial rh.pial.asc
mris_convert --cras_add --to-scanner rh.white rh.white.asc


# now adjust the files above (remove lines etc)

awk '//{if (NR>2){print $0}}' < lh.pial.asc > lh.pial-2.asc
awk '//{if (NR>2){print $0}}' < lh.white.asc > lh.white-2.asc
awk '//{if (NR>2){print $0}}' < rh.pial.asc > rh.pial-2.asc
awk '//{if (NR>2){print $0}}' < rh.white.asc > rh.white-2.asc


# get curvature information in ascii format

mris_convert -c lh.avg_curv lh.pial lh.pial.avg_curv.asc
mris_convert -c lh.avg_curv lh.white lh.white.avg_curv.asc
mris_convert -c rh.avg_curv rh.pial rh.pial.avg_curv.asc
mris_convert -c rh.avg_curv rh.white rh.white.avg_curv.asc

rm lh.pial lh.white rh.pial rh.white


##
## For left hemipshere
##

#
# Now construct the small R script that is used to select the points (for both the pial and the wm surface
# that are part of the particular label.

##
## NOTE THAT WE NEGATE THE X AND Y COORDINATES!!
##

echo "# R script to select points from inner and out surface given their index in a file." > selectPoints_l.R
echo "# Assuming that the first 2 lines of the label and surface files are removed." >> selectPoints_l.R
echo "indices <- read.table(\"tmpLabel_l\")" >> selectPoints_l.R
echo "pial <- read.table(\"lh.pial-2.asc\");" >> selectPoints_l.R
echo "i<-indices[,1]" >> selectPoints_l.R
echo "i=i+1; # indices start at 0 while R starts at 1" >> selectPoints_l.R
echo "pialSelected<-pial[i,1:3]" >> selectPoints_l.R
echo "pialSelected[,1]<-pialSelected[,1]*-1" >> selectPoints_l.R
echo "pialSelected[,2]<-pialSelected[,2]*-1" >> selectPoints_l.R
echo "write.table(pialSelected, file=\"PIALOUTFILE_l\",col.names=FALSE, row.names=FALSE)"  >> selectPoints_l.R
echo "wmSelected<-pial[i,1:3]" >> selectPoints_l.R
echo "wmSelected[,1]<-pialSelected[,1]*-1" >> selectPoints_l.R
echo "wmSelected[,2]<-pialSelected[,2]*-1" >> selectPoints_l.R
echo "write.table(wmSelected, file=\"WHITEOUTFILE_l\",col.names=FALSE, row.names=FALSE)" >> selectPoints_l.R

# Now for each label file in the label directory do:

for x in `ls ../label/lh.*.label`
do
   # first line removes the two first lines because they contain header info,
   awk '//{if (NR>2){print $0}}' < ${x} > tmpLabel_l
   Rscript selectPoints_l.R
   BASE=`basename $x .label`
   mv PIALOUTFILE_l ${BASE}.pial
   mv WHITEOUTFILE_l ${BASE}.white
   rm tmpLabel_l
done

rm selectPoints_l.R

##
## for the right hemisphere
##

#
# Now construct the small R script that is used to select the points (for both the pial and the wm surface
# that are part of the particular label.

##
## NOTE THAT WE NEGATE THE X AND Y COORDINATES!!
##

echo "# R script to select points from inner and out surface given their index in a file." > selectPoints_r.R
echo "# Assuming that the first 2 lines of the label and surface files are removed." >> selectPoints_r.R
echo "indices <- read.table(\"tmpLabel_r\")" >> selectPoints_r.R
echo "pial <- read.table(\"rh.pial-2.asc\");" >> selectPoints_r.R
echo "i<-indices[,1]" >> selectPoints_r.R
echo "i=i+1; # indices start at 0 while R starts at 1" >> selectPoints_r.R
echo "pialSelected<-pial[i,1:3]" >> selectPoints_r.R
echo "pialSelected[,1]<-pialSelected[,1]*-1" >> selectPoints_r.R
echo "pialSelected[,2]<-pialSelected[,2]*-1" >> selectPoints_r.R
echo "write.table(pialSelected, file=\"PIALOUTFILE_r\",col.names=FALSE, row.names=FALSE)"  >> selectPoints_r.R
echo "wmSelected<-pial[i,1:3]" >> selectPoints_r.R
echo "wmSelected[,1]<-pialSelected[,1]*-1" >> selectPoints_r.R
echo "wmSelected[,2]<-pialSelected[,2]*-1" >> selectPoints_r.R
echo "write.table(wmSelected, file=\"WHITEOUTFILE_r\",col.names=FALSE, row.names=FALSE)" >> selectPoints_r.R

# Now for each label file in the label directory do:

for x in `ls ../label/rh.*.label`
do
   # first line removes the two first lines because they contain header info,
   awk '//{if (NR>2){print $0}}' < ${x} > tmpLabel_r
   Rscript selectPoints_r.R
   BASE=`basename $x .label`
   mv PIALOUTFILE_r ${BASE}.pial
   mv WHITEOUTFILE_r ${BASE}.white
   rm tmpLabel_r
done

rm selectPoints_r.R

# and return back to original dir 
cd $CURRENTDIR
