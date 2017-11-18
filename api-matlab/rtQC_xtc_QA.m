function template_hdr = rtQC_xtc_QA(xtcDir,dcmDir,niixDir,ExpCode,SubCode,nVolumes,nSlices)


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
Description: This procedure is used for Quality Assurance of realtime data 
acquired using the Philips XTC module. It ensures appropriate header information
during the neurofeedback as well as active scaling for optimal dynamic
range.

The function relies on a voxelwise comparison of the timeseries obtained both offline
in 'Philips classic dicom' format and online using the XTC real-time 
export module. The idea is to run it on your functional localizer data 
and to then use the output to ensure proper header info during a subsequent
real-time neurofeedback acquisition. The two sequences (functional 
localizer and real-time neurofeedback) must have identical scanning
parameters and be obtained during the same scanning session. It is also 
recommended to use the last volume from the functional localizer as the 
reference volume for movement correction during your realtime analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Usage:  

template_hdr=rtQC_XTC_QA(xtcDir,dcmDir,niixDir,ExpCode,SubCode,nVolumes,nSlices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Arguments, inputs and outputs:

INPUTS

xtcDir ~ path to stored data from the functional localizer, obtained using 
the online Philips XTC real-time data module.

dcmDir ~ path to stored data from the functional localizer, exported using
the 'Philips classic dicom' offline export option.

niixDir ~ path to where you want the NIfTI version of the dicom data to be
stored.

ExpCode ~ Experiment Code; a string with the name of your experiment.

SubCode ~ Subject Code; a string with a code for your subject. 

nVolumes ~ number of volumes in your sequence.

nSlices ~ number of slices in your sequence.

OUTPUT

template_hdr ~ .img header file that you can use when converting XTC real-time
PAR/REC data to NIfTI during neurofeedback with a scanning sequence  
identical to the functional localizer from the same scanning session.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Special notes: 

a)To add soon: storing timing data and printing at the end of the sequence.

b) The 0.02 data discrepancy threshold, for accepted deviation between the 
online and offline datasets was established empirically. The dynamic range 
in this case is {0 12000} so 0.02/12000 is a negligible difference.
Next versions should include an automated method to determine this threshold. 

c) Dealing with the rare case of varying PAR/REC slice numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dependencies: dcm2nii from MRIcro, loadPARREC.m and SPM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%}



copyfile([dcmDir '/*'], ['D:\NPAD\' ExpCode '/DATA/' SubCode '/realtime_localizer/'])

if ~exist([niixDir '\rtfMRI_' ExpCode '_template_hdr_file.nii'],'file')

    % Step 2: Check the correspondence between the online (par/rec) format and offline (dicom) versions of the same sequence acquisition.
    
        % create single-volume files for whole sequence with identical values to par/rec versions
        pwDir=pwd; cd D:\NPAD\Resources\Software\MRIcroGL\mricrogl
        command = ['dcm2niix -o ' [niixDir] ' ' dcmDir] ; system(command); cd(pwDir);
        fnii4D=ls([niixDir '\*.nii']); flist = spm_file_split([niixDir '\' fnii4D(1,:)]); 

        % Check data
        for volume = 1:nVolumes;
            Vol=num2str(volume)
            ParRecVol=num2str(volume-1);
            
            % load reference data
            inpFile=flist(volume).fname; infoVol = spm_vol(inpFile); imgVol  = spm_read_vols(infoVol);
            
            % load Par/Rec data while fixing orientation
            [DATA INFO]=loadPARREC([xtcDir '/Dump-' sprintf('%04d',(volume-1)) '.par']);
            for slice=1:nSlices
                PAR(:,:,slice) = fliplr(DATA(:,:,slice,1,1,1)');
            end
            
            % compare data, allowing for minimum discrepancy & display/report any inconsistencies
            if nnz(abs(imgVol(:)-PAR(:))>0.02)>0;  % 0.02 is the Data Discrepancy Threshold 
                    display(['PROBLEM WITH SLICE ' num2str(slice) ' IN VOLUME ' num2str(volume)]); 
            end
        end

    % % Step 3: Obtain a hdr template for the specific sequence you are going to be using in the realtime paradigm.
    % template_hdr_file = get_template_hdr(nii_folder_with_exactly_same_data)
    infoVol.fname=[niixDir '\rtfMRI_' ExpCode '_template_hdr_file.nii'];
    template_hdr=[niixDir '\rtfMRI_' ExpCode '_template_hdr_file.nii'];
    infoVol=rmfield(infoVol,'pinfo'); spm_write_vol(infoVol,PAR);

end

% Step 4: Insert the online conversion function "rtQC_xtc2nii(xtc_folder, nii_folder, template_hdr_file)" into your pipeline.
