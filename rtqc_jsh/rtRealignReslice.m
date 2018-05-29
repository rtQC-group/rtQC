function [rfunctional_fn, MP] = rtRealignReslice(ref_fn, functional_fn)

realign = struct;
% Data
data={};
data{1}=[ref_fn ',1'];
data{2}=[functional_fn ',1'];
realign.matlabbatch{1}.spm.spatial.realign.estwrite.data = {data'};
% Eoptions
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9; % could be lowered to increase speed
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0; % only a one pass procedure
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
realign.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
% Roptions
realign.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [1 0]; %[2 0]: All images (1..n); [1 0]: Images 2..n
realign.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
realign.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
realign.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
realign.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% Run
cfg_util('run',realign.matlabbatch);
[d, fn, ext] = fileparts(functional_fn);
rfunctional_fn = [d filesep 'r' fn ext];

