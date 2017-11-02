function rtQC_xtc2nii(xtcDir, niiOutDir, template_hdr_file, volume, nSlices)  

%{ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Function name: rtQC_xtc2nii

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Description: This function converts data acquired using the Philips XTC
 real-time export module (PAR/REC files) to NIfTI format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Usage: rtNII_volume = rtQC_xtc2nii(xtc_ParRec_volume, template_volume)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Arguments, inputs and outputs:

xtc_ParRec_volume ~ Realtime volume dumped by XTC in Par/Rec format

template_volume ~ a volume acquired using the same scanning sequence that 
is used as reference to complete missing header attributes

rtNII_volume ~ typical nii file you can readily use in SPM or other 
software, corresponding exactly to what you would have gotten if you had 
exported the data offline into dicom and then converted it to nii using 
dcm2nii from MRIcro.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Special notes: 
To add soon: storing timing data and printing at the end of the sequence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dependencies: dcm2nii from MRIcro, loadPARREC.m and SPM.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Stavros Skouras; Barcelonabeta Brain Research Center; 17/02/2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE: GPLv3 / please cite author
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Use at your own risk!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}


    [DATA INFO]=loadPARREC([xtcDir '/Dump-' sprintf('%04d',(volume-1)) '.par']);
    for slice=1:nSlices
        PAR(:,:,slice) = fliplr(DATA(:,:,slice,1,1,1)');
    end
    infoVol = spm_vol(template_hdr_file);
    infoVol.fname=[niiOutDir '\Volume' sprintf('%03d',volume) '.nii'];
    infoVol=rmfield(infoVol,'pinfo'); spm_write_vol(infoVol,PAR);
    

return
