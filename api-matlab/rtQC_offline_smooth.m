function sfunctional4D_fn = rtQC_offline_smooth(functional4D_fn, Nt, fwhm)
smooth = struct;
% Data
fns={};
for i = 1:Nt
    fns{i} = [functional4D_fn ',' num2str(i) ];
end
smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
% Other
smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% Run
cfg_util('run',smooth.matlabbatch);
[d, f, e] = fileparts(functional4D_fn);
sfunctional4D_fn = [d filesep 's' f e];