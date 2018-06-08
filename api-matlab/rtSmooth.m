function sfunctional_fn = rtSmooth(functional_fn, fwhm)
% Initialize
smooth = struct;
% Data
smooth.matlabbatch{1}.spm.spatial.smooth.data={functional_fn};
% Other
smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = fwhm;
smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% Run
cfg_util('run',smooth.matlabbatch);
% Output
[d,f,e]=fileparts(functional_fn);
sfunctional_fn = [d filesep 's' f e];