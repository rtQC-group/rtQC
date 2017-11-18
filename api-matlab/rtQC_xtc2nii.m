function rtQC_xtc2nii(xtcDir, niiOutDir, template_hdr_file, volume, nSlices)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) rtQC Dev-Team
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Stavros Skouras; Barcelonabeta Brain Research Center; 17/02/2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
This can be accomplished using the function rtQC_xtc_QA in advance; e.g.
template_hdr = rtQC_xtc_QA(xtcDir,dcmDir,niixDir,ExpCode,SubCode,nVolumes,nSlices)

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
%}



    [DATA INFO]=loadPARREC([xtcDir '/Dump-' sprintf('%04d',(volume-1)) '.par']);
    for slice=1:nSlices
        PAR(:,:,slice) = fliplr(DATA(:,:,slice,1,1,1)');
    end
    infoVol = spm_vol(template_hdr_file);
    infoVol.fname=[niiOutDir '\Volume' sprintf('%03d',volume) '.nii'];
    infoVol=rmfield(infoVol,'pinfo'); spm_write_vol(infoVol,PAR);
    

return
