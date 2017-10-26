

Function name: rtQC_xtc2nii

Description: This function converts data acquired from the Philips XTC real-time export module (PAR/REC files) to Nifti format.

Usage: rtNII_volume = rtQC_xtc2nii(xtc_ParRec_volume, template_volume)

Arguments, inputs and outputs:

xtc_ParRec_volume ~ Realtime volume dumped by XTC in Par/Rec format

template_volume ~ a volume acquired using the same scanning sequence that is used as reference to complete missing header attributes

rtNII_volume ~ typical nii file you can readily use in SPM or other software, corresponding exactly to what you would have gotten if you had exported the data offline into dicom and then converted it to nii using dcm2nii from MRIcro.

Special notes: not applicable

Dependencies: dcm2nii from MRIcro

Contributors: Stavros



%% ability to compare stuff --

% some history tracking in rtqc structure
% 
% 