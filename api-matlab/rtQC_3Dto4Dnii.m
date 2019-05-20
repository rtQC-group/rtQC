function rtQC_3Dto4Dnii(vols, TR, name)
% Initialize
spm('defaults','fmri');
spm_jobman('initcfg');
conversion = struct;
% Vols
conversion.matlabbatch{1}.spm.util.cat.vols = vols;
% Name
conversion.matlabbatch{1}.spm.util.cat.name = name;
% dtype
conversion.matlabbatch{1}.spm.util.cat.dtype = 4;
% dtype
conversion.matlabbatch{1}.spm.util.cat.RT = TR;
% Run
spm_jobman('run',conversion.matlabbatch);