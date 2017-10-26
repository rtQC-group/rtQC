#!/bin/bash

# SetupTopup
workdir=TopupFiles


# the helper text,if no arguments are passed:
function showhelp() {
	echo 'Usage: '
	echo 'Call this script before after you have acquired the AP and PA images and before '
	echo 'running the rt-fMRI experiment as follows:'
	echo
	echo 'RT_SetupTopup <NAME OF AP IMAGE> <NAME OF PA IMAGE>'
	echo
	echo "It will make a directory named $workdir with all Topup file in it"
	echo 'This will be the basis for the follow-up script:'
	echo 'RT_ApplyTopup'
	
}




# first argument should be ap image - second argument the pa image.
if [ -z $1 ]; then
	showhelp()
	exit 1
fi


# copy/define variables
cp $1 $workdir/ap.nii.gz
cp $2 $workdir/pa.nii.gz
# use these from now on...
ap=$workdir/ap.nii.gz
pa=$workdir/pa.nii.gz

# check the # of slices
echo 'checking whether # slices are divisible by 2...'
nslices=`fslval $ap dim3`
rem=`bc <<< $nslices%2`
if [ ! $rem -eq 0 ]; then
	echo ' removing last slice, since all dims needs to be divisible by 2... '
	xsize=`fslval $ap dim1`
	ysize=`fslval $ap dim2`
	
	# use fslroi to snoop off the last slice, of # slices is not divisible... then estimate topup after
	fslroi $ap $ap 0 $xsize 0 $ysize 0 $((nslices-1))
 	fslroi $pa $pa 0 $xsize 0 $ysize 0 $((nslices-1))
	
fi


# checking whether equal amount of images in ap and pa (common mistake)
# must be equal
dim4ap = `fslval $ap dim4`
dim4pa = `fslval $pa dim4`

if [ ! $dim4ap -eq $dim4pa ];then
	echo 'dimensions of your AP and PA scans do NOT match!!!! - fix this first!'
	exit 1
fi

# requirement #1 for topup: the merged file
fslmerge -t $workdir/appa $ap $pa


# requirement #2 for topup: now make an acqparams file.
acqparams=$workdir/acqparams.txt
dwelltime = 0.080 # just a dummy value close to typical values. It doesn't matter in our case.
# but you might wish to change for other application(s).
touch $acqparams
for iter in `seq 1 $dim4ap`;do
	echo 0 1 0 $dwelltime >> $acqparams
fi
for iter in `seq 1 $dim4pa`;do
	echo 0 -1 0 $dwelltime >> $acqparams
fi


# now we can actually calculate the topup!!
topup --imain=$workdir/appa --datain=acqparams.txt --config=$FSLDIR/etc/flirtsch/b02b0.cnf --out=topup_results --iout=topup_Magnitudes --fout=topup_TopupField --dfout=se_topup_WarpField --rbmout=topup_MotionMatrix --jacout=topup_Jacobian -v




# done!
