#!/bin/bash

# SetupTopup
workdir=TopupFiles
 


# the helper text,if no arguments are passed:
function showhelp() {
	echo 'RT_SetupTopup.sh'
	echo '-------------'
	echo 
	echo 'Call this script before after you have acquired the AP and PA images and before '
	echo 'running the rt-fMRI experiment as follows:'
	echo
	echo '<PATH-TO-RT_SetupTopup>/RT_SetupTopup.sh <NAME OF AP IMAGE> <NAME OF PA IMAGE>'
	echo
	echo 'It will make a directory with all Topup file in it in:'
	echo `pwd`"/$workdir"
	echo 
	echo 'These files will be used by the follow-up script:'
	echo
	echo 'RT_ApplyTopup.sh'
	echo
	echo 'If you wish to know the total running time, call prepended with time, like so: '
	echo 'time <PATH-TO-RT_SetupTopup>/RT_SetupTopup.sh <NAME OF AP IMAGE> <NAME OF PA IMAGE>'
	echo 'average running time for 1GE and 1SE image of 80x80x44 with TR of ~3sec: 3minutes'
	echo 'on Dell XPS13 notebook with SDD'
	echo 'Times on other systems may vary.'
	echo 
	echo
	echo 'It is probably best is you can also use the MRI to export the GE (or SE) fieldmaps'
	echo 'directly into MRIncoming'
	echo 'So when you have obtained the AP and PA scans, you can immedeately start this script,'
	echo 'wait for ~3 minutes (instruct the subject?), and then start the experiment'
	echo
	
	
}

# check whether FSL environment variable is set - otherwise instruct user
# to do it. 
if [ -z ${FSLDIR+x} ]; then 
	echo "\$FSLDIR is unset - do you have FSL installed (and its environment variable set)?" 
else 
	echo "Found FSL! - \$FSLDIR is set to '$FSLDIR'"; 
fi


# first argument should be ap image - second argument the pa image.
if [ -z $1 ]; then
	showhelp
	exit 1
fi


if [ -d "$workdir" ]; then
	echo 'The workdir already exists! -- remove it to be sure I do not overwrite anything'
	echo 'and call this script again.'
	exit 1
else
	mkdir $workdir
fi



# copy/define variables
# use fslmaths to copy files (and gzip them at the same time) - multiply by 1.
fslmaths $1 -mul 1 $workdir/ap.nii.gz
fslmaths $2 -mul 1 $workdir/pa.nii.gz
# use these from now on...
ap=$workdir/ap.nii.gz
pa=$workdir/pa.nii.gz

# check the # of slices
echo 'checking whether # slices are divisible by 2...'
nslices=`fslval $ap dim3`
rem=`bc <<< $nslices%2`
if [ ! $rem -eq 0 ]; then
	echo '!! removing last slice, since all dims needs to be divisible by 2... '
	xsize=`fslval $ap dim1`
	ysize=`fslval $ap dim2`
	
	# use fslroi to snoop off the last slice, of # slices is not divisible... then estimate topup after
	fslroi $ap $ap 0 $xsize 0 $ysize 0 $((nslices-1))
 	fslroi $pa $pa 0 $xsize 0 $ysize 0 $((nslices-1))
	
fi


# checking whether equal amount of images in ap and pa (common mistake)
# must be equal
dim4ap=`fslval $ap dim4`
dim4pa=`fslval $pa dim4`

if [ ! $dim4ap -eq $dim4pa ];then
	echo 'dimensions of your AP and PA scans do NOT match!!!! - fix this first!'
	exit 1
else
	echo $dim4ap ' timepoints in AP scan'
	echo $dim4pa ' timepoints in PA scan'
fi
	

# requirement #1 for topup: the merged file
fslmerge -t $workdir/appa $ap $pa


# requirement #2 for topup: now make an acqparams file.
acqparams=$workdir/acqparams.txt
dwelltime=0.080 # just a dummy value close to typical values. It doesn't matter in our case.
# but you might wish to change for other application(s).
touch $acqparams
for iter in `seq 1 $dim4ap`; do
	echo "0 1 0 $dwelltime" >> $acqparams
done
for iter in `seq 1 $dim4pa`; do
	echo "0 -1 0 $dwelltime" >> $acqparams
done


# now we can actually calculate the topup!!
echo "Starting Topup... making topup results for follow-up Script..."
topup --imain=$workdir/appa --datain=$acqparams --config=$FSLDIR/etc/flirtsch/b02b0.cnf --out=$workdir/topup_results --iout=$workdir/topup_Magnitudes --fout=$workdir/topup_TopupField --dfout=$workdir/topup_WarpField --rbmout=$workdir/topup_MotionMatrix --jacout=$workdir/topup_Jacobian -v

# generate the mcref image
fslroi $workdir/ap $workdir/mcref 0 1

# now - generate the t1ref image!
echo 'Making Reference Image for follow-up Script...'
applytopup --imain=$workdir/ap.nii.gz --inindex=1 --datain=$workdir/acqparams.txt --topup=$workdir/topup_results --out=$workdir/ref.nii.gz --method=jac




# done!
