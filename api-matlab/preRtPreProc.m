function output = preRtPreProc(functional0_fn, structural_fn, spm_dir)
% Function to complete pre-real-time preprocessing of structural and
% functional data from a single subject. Steps include coregistering
% structural image to initial functional image, segmenting the coregistered
% structural image into tissue types, and reslicing the segments to the
% functional resolution image grid. Makes use of spm12 batch routines. 
% If spm12 batch parameters are not explicitly set, defaults are assumed.
%
% INPUT:
% funcional0_fn     - filename of pre-real-time functional scan
% structural_fn     - filename of T1-weighted structural scan
% structural_fn     - filename of T1-weighted structural scan
% 
% OUTPUT: 
% output            - structure with filenames and data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis


output = struct;

% STEP 1 -- Coregister structural image to first dynamic image (estimate and reslice)
disp('1 - Coregistering structural to functional image space...');
coreg_estimate = struct;
% Ref
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[functional0_fn ',1']};
% Source
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = {structural_fn};
% Other
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estwrite.other = {};
% Eoptions
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% % Roptions
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
% Run
cfg_util('run',coreg_estimate.matlabbatch);
disp('done');

% STEP 2 -- Segmentation of coregistered structural image into GM, WM, CSF, etc
% (with implicit warping to MNI space, saving forward and inverse transformations)
disp('2 - Segmenting coregistered structural image into GM, WM, CSF, etc...');
segmentation = struct;
% Channel
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.vols = {structural_fn};
% Tissue
for t = 1:6
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,' num2str(t)]};
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = t-1;
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
% Warp
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.write=[1 1];
% Run
cfg_util('run',segmentation.matlabbatch);
% Save filenames
[d, fn, ext] = fileparts(structural_fn);
output.forward_transformation = [d filesep 'y_' fn '.nii'];
output.inverse_transformation = [d filesep 'iy_' fn '.nii'];
output.gm_fn = [d filesep 'c1' fn '.nii'];
output.wm_fn = [d filesep 'c2' fn '.nii'];
output.csf_fn = [d filesep 'c3' fn '.nii'];
output.bone_fn = [d filesep 'c4' fn '.nii'];
output.soft_fn = [d filesep 'c5' fn '.nii'];
output.air_fn = [d filesep 'c6' fn '.nii'];
disp('done');

% STEP 3 -- Reslice all to functional-resolution image grid
disp('3 - Reslice all generated images to functional-resolution image grid');
reslice = struct;
% Ref
reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {[functional0_fn ',1']};
% Source
source_fns = {};
for i = 1:6
    source_fns{i} = [d filesep 'c' num2str(i) fn ext];
end
source_fns{7} = structural_fn;
reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
% Roptions
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% Run
cfg_util('run',reslice.matlabbatch);
% Save filenames
output.rstructural_fn = [d filesep 'r' fn ext];
output.rgm_fn = [d filesep 'rc1' fn '.nii'];
output.rwm_fn = [d filesep 'rc2' fn '.nii'];
output.rcsf_fn = [d filesep 'rc3' fn '.nii'];
output.rbone_fn = [d filesep 'rc4' fn '.nii'];
output.rsoft_fn = [d filesep 'rc5' fn '.nii'];
output.rair_fn = [d filesep 'rc6' fn '.nii'];
disp('done');
