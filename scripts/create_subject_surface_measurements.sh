# set -x
if [ $# -ne 1 ]
then
echo
echo "wrong number of parameters."
echo
echo  "usage: create_subject_surface_measurements_II <subject's freesurfer directory>"
echo
echo ""   
echo "Creates surface files (pial and white) that are ready to be used with mincsamplecortex (found in layer directory)" 
echo "In addition tihs version (II) also stores the curvature information in a separate file, which can be included via mincsamplecortex."  
echo "These surfaces are in the orientation of brain.finalsurfs.mgz for which the origin is set to -128,-128,-128"
echo "Here a target brain is placed in the layers directory and is named target.mnc"   
echo "Thus, the transformation that you need for mincsamplesurface is the transformation that registers your original anatomy used by freesurfer scan towards the volume (eg DTI) that contains"
echo "the cortexdata to be sampled. (Thus if you want to resample the original anatomy scan (T1) then you provide the identity matrix"
echo "(The transformation that registers the original anatomy to this target brain is named orig_anat_to_target.xfm)" 
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

# remove first 2 lines from ASCII pial and white files

mris_convert lh.pial lh.pial.asc
mv lh.pial.asc tmp.txt
cat tmp.txt | awk '//{if (NR > 2) print $0}' > lh.pial.asc
rm tmp.txt
mris_convert lh.white lh.white.asc
mv lh.white.asc tmp.txt
cat tmp.txt | awk '//{if (NR > 2) print $0}' > lh.white.asc
rm tmp.txt
mris_convert rh.pial rh.pial.asc
mv rh.pial.asc tmp.txt
cat tmp.txt | awk '//{if (NR > 2) print $0}' > rh.pial.asc
rm tmp.txt
mris_convert rh.white rh.white.asc
mv rh.white.asc tmp.txt
cat tmp.txt | awk '//{if (NR > 2) print $0}' > rh.white.asc
rm tmp.txt

mris_convert -c lh.avg_curv lh.pial lh.pial.avg_curv.asc
mris_convert -c lh.avg_curv lh.white lh.white.avg_curv.asc
mris_convert -c rh.avg_curv rh.pial rh.pial.avg_curv.asc
mris_convert -c rh.avg_curv rh.white rh.white.avg_curv.asc

rm lh.pial lh.white rh.pial rh.white lh.avg_curv rh.avg_curv

# now rh avg_curv.asc files contain both coordinates of the surface points as well as the curvature value


#
# Now construct the small octave script that is used to select those
#

echo "# octave script to select points from inner and out surface (and curvature) given their index in a file." > selectPoints.txt
echo "# You should have removed the first 2 lines of the label and surface files." >> selectPoints.txt
echo "indices=load(\"tmpLabel\");" >> selectPoints.txt
echo "pial=load(\"lh.pial.avg_curv.asc\");" >> selectPoints.txt
echo "i=indices(:,1);" >> selectPoints.txt
echo "i=i+1; % indices start at 0, in octave at 1! " >> selectPoints.txt
echo "pialSelected=pial(i,2:4);" >> selectPoints.txt
echo "curvatureSelected=pial(i,5);" >> selectPoints.txt
echo "save(\"-ascii\",\"PIALOUTFILE\",\"pialSelected\");" >> selectPoints.txt
echo "save(\"-ascii\",\"CURVEATUREOUTFILE\",\"curvatureSelected\");" >> selectPoints.txt
echo "white=load(\"lh.white.avg_curv.asc\");" >> selectPoints.txt
echo "whiteSelected=white(i,2:4);" >> selectPoints.txt
echo "save(\"-ascii\",\"WHITEOUTFILE\",\"whiteSelected\");" >> selectPoints.txt


# copy 
#

# again, don't forget to remove first 2 lines of label tex tfiles

for x in lh.BA1.label lh.BA2.label lh.BA3a.label lh.BA3b.label lh.BA44.label lh.BA45.label lh.BA4a.label lh.BA4p.label lh.BA6.label lh.MT.label lh.V1.label lh.V2.label lh.cortex.label lh.entorhinal_exvivo.label
do
   awk '//{if (NR>2){print $0}}' < ../label/${x} > tmpLabel
   octave selectPoints.txt
   BASE=`basename $x .label`
   mv PIALOUTFILE ${BASE}.pial
   mv WHITEOUTFILE  ${BASE}.white
   mv CURVEATUREOUTFILE ${BASE}.curvature
   rm tmpLabel
done

rm selectPoints.txt


echo "# octave script to select points from inner and out surface (and curvature) given their index in a file." > selectPoints.txt
echo "# You should have removed the first 2 lines of the label and surface files." >> selectPoints.txt
echo "indices=load(\"tmpLabel\");" >> selectPoints.txt
echo "pial=load(\"rh.pial.avg_curv.asc\");" >> selectPoints.txt
echo "i=indices(:,1);" >> selectPoints.txt
echo "i=i+1; % indices start at 0, in octave at 1! " >> selectPoints.txt
echo "pialSelected=pial(i,2:4);" >> selectPoints.txt
echo "curvatureSelected=pial(i,5);" >> selectPoints.txt
echo "save(\"-ascii\",\"PIALOUTFILE\",\"pialSelected\");" >> selectPoints.txt
echo "save(\"-ascii\",\"CURVEATUREOUTFILE\",\"curvatureSelected\");" >> selectPoints.txt
echo "white=load(\"rh.white.avg_curv.asc\");" >> selectPoints.txt
echo "whiteSelected=white(i,2:4);" >> selectPoints.txt
echo "save(\"-ascii\",\"WHITEOUTFILE\",\"whiteSelected\");" >> selectPoints.txt


for x in rh.BA1.label rh.BA2.label rh.BA3a.label rh.BA3b.label rh.BA44.label rh.BA45.label rh.BA4a.label rh.BA4p.label rh.BA6.label rh.MT.label rh.V1.label rh.V2.label rh.cortex.label rh.entorhinal_exvivo.label 
do
   awk '//{if (NR>2){print $0}}' < ../label/${x} > tmpLabel
   octave selectPoints.txt
   BASE=`basename $x .label`
   mv PIALOUTFILE ${BASE}.pial
   mv WHITEOUTFILE  ${BASE}.white
   mv CURVEATUREOUTFILE ${BASE}.curvature
   rm tmpLabel
done

rm selectPoints.txt

# Now create target.mnc volume
mri_convert -i ../mri/brain.finalsurfs.mgz -o target.mnc
cp target.mnc tmp.mnc
minc_modify_header -dinsert xspace:start=-128 target.mnc
minc_modify_header -dinsert yspace:start=-128 target.mnc
minc_modify_header -dinsert zspace:start=-128 target.mnc

# now compute the transformation that transforms the original anatomy scan to the adjusted target.mnc brain.
minctracc -lsq3 tmp.mnc  target.mnc orig_anat_to_target.xfm
rm tmp.mnc
# and return back to original dir 
cd $CURRENTDIR
