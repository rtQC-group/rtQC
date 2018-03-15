% rtQC_correct_field_distortion 	(provisional feature – subject to successful implementation)
% Corrects field distortion volume by volume in real-time.
% Current developers: Johan/Stavros



% this function should ask several questions (i.e. popup gui's that allow
% you to select some files, etc) and then start 'watching' a directory for
% subsequent distortion (+motion) correction.
%
% we provide an optional data set to play around with, as well as a
% README.md file (and a guide) to show 'normal' operation.
%
% test data can be acquired at: /media/data/Dropbox/Prog/RTQA/test
%
% data that we shared with Stavros: /media/data/Dropbox/Temp/RTQA 
%
% More test data (i.e. the 7T:  /media/data/Dropbox/Temp/RTQA
%
% The data format should be in so-called 
%

% rtQC_correct_field_distortion(varargin)
%function rtQC_correct_field_distortion(varargin)


s='Welcome to real-time motion + distortion correction.\nFor the best results, please remember to use Spin-Echo Ap/PA Scans\ninstead of Gradient-Echo ones due to better WM/GM contrast.\n';
fprintf(s);


% check whether SPM is on the path
% check whether these scripts are on the path, too.




% ask where the ap fieldmap is
[f_ap, ~] = spm_select(1,'image','Select AP Fieldmap',{},pwd,'.*');
f_ap=regexprep(f_ap,',1$','');

% ask where the pa fieldmap is
[f_pa, ~] = spm_select(1,'image','Select PA Fieldmap',{},pwd,'.*');
f_pa=regexprep(f_pa,',1$','');

% ask where images will start to appear.
[p_mri, ~] = spm_select(1,'dir','In which folder will fMRI data appear?',{},pwd,'.*');

% ask if the fmri data is going to be ap or pa encoded
encoding_question = 'Will fMRI data be AP or PA encoded? [A/P] ';
ph_encoding='None';
while strcmp(ph_encoding,'None')
    s=input(encoding_question,'s');
    if strcmp(s,'')
        ph_encoding = 'None';
    elseif strcmpi(s,'a')
        ph_encoding = 'A';
        applytopup_arg = 1;
    elseif strcmpi(s,'p')
        ph_encoding = 'P';
        applytopup_arg = 2;
    end
    encoding_question = '... AP or PA encoded? [A/P] ';
end


% I should do this .. the bash scripts will do this already - and also
% convert to fsl .nii.gz format using fslmaths -mul 1
% now make IN the 'incoming' directory, the directory TopupFiles
%fprintf('... creating TopupFiles Directory in %s\n',p_mri);
topupdir = [p_mri filesep 'TopupFiles'];
%mkdir(topupdir)

fprintf('... copying ap and pa fieldmaps into this directory\n');

try
    copyfile(f_ap,[p_mri filesep '.']);
catch
    disp('file already at right place');
end
try
    copyfile(f_pa,[p_mri filesep '.']);
catch
    disp('file already at right place');
end
[~,n,e]=fileparts(f_ap);
f_ap_onlyname = [n e];
[~,n,e]=fileparts(f_pa);
f_pa_onlyname = [n e];


% figure out from where we'd need to call the scripts;
tmp = which('rtQC_correct_field_distortion');
tmp = regexprep(tmp,[filesep '[^' filesep ']*$'],'');
p_scripts = [tmp filesep 'fsl-topup-files'];
p_scripts = regexprep(p_scripts,['api-matlab' filesep],'');

fprintf('... running RT_SetupTupup.sh script on command line\n    to create necessary files for RT operation\n');
s= [p_scripts filesep 'RT_SetupTopup.sh ' f_ap ' ' f_pa];


% then do the prepare_topup script
% change-dir INTO the directory just below the MRINCOMING folder
cd(p_mri);
system(['cd ' p_mri]);

s=[p_scripts filesep 'RT_SetupTopup.sh ' p_mri filesep 'ge_ap.nii ' p_mri filesep 'ge_pa.nii'];
disp(['... running ' s]);
system(s);



% then start up the folder-monitoring (it will operate on any NEW incoming 
% .nii file)
% use the ampersand to do it in the background.
s=[p_scripts filesep 'RT_ApplyTopup.sh ' topupdir ' ' p_mri ' ' num2str(applytopup_arg) ''];
disp(['... running ' s]);
system(s,'-echo');







% after checking - think of how to incorporate this into the RTQA folder
% structure (to show what you did with this dataset?)

