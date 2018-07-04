% test_spm_realign_rt

% STEPHAN CODE
% Access data
data_dir = '/Users/jheunis/Documents/MATLAB/data/';
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
subject = '0051210';
functional4D_fn = [data_dir filesep subject filesep 'rest_1' filesep 'rest.nii'];

% FROM OPENNFT: setupProcParams.m
% SPM flags
flagsSpmRealign = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6);
flagsSpmReslice = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'mask',1,'mean',0,'which', 2);
% Get motion realignment template data and volume
infoVolTempl = spm_vol([functional4D_fn ',1']);
mainLoopData.infoVolTempl = infoVolTempl;
imgVolTempl  = spm_read_vols(infoVolTempl);
dimTemplMotCorr     = infoVolTempl.dim;
matTemplMotCorr     = infoVolTempl.mat;
% Realign preset
A0=[];x1=[];x2=[];x3=[];wt=[];deg=[];b=[];
R(1,1).mat = matTemplMotCorr;
R(1,1).dim = dimTemplMotCorr;
R(1,1).Vol = imgVolTempl;

% STEPHAN CODE - Don't skip any initial volumes
P.nrSkipVol = 0;

% STEPHAN CODE
% Loop through all volumes in 4D functional dataset
for indVol = 1:120
    
    % Access current 'real-time' volume and assign to second index of R
    currentVol = spm_vol([functional4D_fn ',' num2str(indVol)]);
    R(2,1).mat = currentVol.mat;
    R(2,1).dim = currentVol.dim;
    R(2,1).Vol = spm_read_vols(currentVol);
    
    % realign (FROM OPENNFT: preprVol.m)
    [R, A0, x1, x2, x3, wt, deg, b, nrIter] = spm_realign_rt(R, flagsSpmRealign, indVol, P.nrSkipVol + 1, A0, x1, x2, x3, wt, deg, b);
    
    % MC params (FROM OPENNFT: preprVol.m). STEPHAN NOTE: I don't understand this part, but it runs fine
    tmpMCParam = spm_imatrix(R(2,1).mat / R(1,1).mat);
    if (indVol == P.nrSkipVol + 1)
        P.offsetMCParam = tmpMCParam(1:6);
    end
    P.motCorrParam(indVol,:) = tmpMCParam(1:6)-P.offsetMCParam; % STEPHAN NOTE: I changed indVolNorm to indVol due to error, not sure if this okay or wrong?
    % P.motCorrParam(indVolNorm,:) = tmpMCParam(1:6)-P.offsetMCParam;
    
    % Reslice (FROM OPENNFT: preprVol.m)
    reslVol = spm_reslice_rt(R, flagsSpmReslice);
    
end


