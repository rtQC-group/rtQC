function rtqc_jsh()
% Function that creates a GUI that allows the user to execute and visualise
% several aspects of a real-time (RT) fMRI dataset for a single subject.
% Real-time is simulated, i.e. dynamic processing is done per image of an
% existing dataset. Functioning includes:
% 
%   - Enter dataset and directories (in this m-file, prior to running rtqc_jsh)
%   - Run pre-real-time preprocessing in preparation for RT (or check if it
%     has been run previously), which does the following:
%       - Checks if T1w segments already exist
%       - If not, coregisters T1w image to functional0 image, segments
%       the coregistered T1w image into GM, WM and CSF compartments, and
%       creates a GM+WM+CSF mask.
%   - Start realtime simulation, which does the following:
%       - Realigns dynamic image to functional0 image (see inputs below)
%       - Loads dynamic head movement parameters resulting from realignment
%       - Calculates framewise displacement
%       - Calculates DVARS
%       - Smooths the dynamic image
%       - Detrends the masked voxel time series using a cumulative GLM
%       - Calculates percentage signal change for the masked voxels, using
%       an iteratively updated mean.
%       - Calculates iteratively updated tSNR
%       - Calculates number of outlier volumes based on specified FD
%       threshold.
%       - Plots: dynamic unprocessed image, tSNR image, HMPs, FD, DVARS,
%       "The PLot".
%
% 
% INPUTS:
% structural_fn     - filename of T1-weighted structural scan
% funcional0_fn     - filename of pre-real-time functional scan, a 3D nifti
%                     file is assumed.
% funcional4D_fn    - filename of main functional scan, a 4D nifti file is
%                     assumed. 
%
% 
% DEPENDENCIES:
% preRtPreProc          - preprocesses pre-real-time data (coregistration of
%                         T1w to functional0, segmentation into tissue types,
%                         reslicing to functional grid)
% createBinarySegments  - constructs 3D binary images for GM, WM and CSF based
%                         on the relative value of GM/WM/CSF tissue probability
%                         maps per voxel.
% rtRealignReslice      - real time realignment of dynamic functional image
%                         to prespecified reference image
% rtSmooth              - real time smoothing of dynamic functional image
% 
% 
% NOTES/ASSUMPTIONS:
% - This function/GUI makes use of spm12 batch scripting
% - rtqc_jsh is experimental and not tested much.
% - rtqc_jsh can be updated to run on actual real-time fMRI data. TODO
% - Nifti's are assumed, code needs to be added if conversion is necessary
% - Little to no attention was given to code optimization to increase
%   execution speed or for modularity. TODO
% - Some real-time methods might be inaccurate or suboptimal compared to
%   offline equivalent (e.g. cumulative mean, scaling, PSC, etc). This
%   needs attention, where known improvements exist. TODO
% - A sensible rework would involve (amongst other steps):
%   1) replacing all preprocessing steps currently dependent on SPM
%   scripts/batches (and which hence take ages to execute) with equivalent
%   steps from OpenNft, if possible.
%   2) replacing cumulative processes with incremental processes, where
%   applicable (e.g. cGLM vs iGLM)
%   3) allowing for different image file format conversion
%__________________________________________________________________________
% Copyright (C) 2018 Neu3CA.org
% Written by Stephan Heunis




global preproc_data rt_data
% Initialise data directories and filenames (TO BE SPECIFIED BY USER)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% data_dir = ''; % Directory where data is saved
% spm_dir = ''; % spm12 directory
% structural_fn = []; % T1-weighted structural nifti
% functional0_fn = []; % 3D initial functional scan
% functional4D_fn = []; % 4D unprocessed fMRI dataset
% rt_data.FD_threshold = 1; % Threshold in mm of framewise displacement, the value at or above which a volume is considered a movement outlier

% Example:
data_dir = '/Users/jheunis/Documents/MATLAB/rtqc_jsh/rtqc_data/0051210'; % Directory where data is saved
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12'; % spm12 directory
structural_fn = [data_dir filesep 'anat_1' filesep 'mprage.nii']; % T1-weighted structural nifti
functional0_fn = [data_dir filesep 'rest_1' filesep 'rest.nii,1']; % 3D initial functional scan
functional4D_fn = [data_dir filesep 'rest_1' filesep 'rest.nii']; % 4D unprocessed fMRI dataset
rt_data.FD_threshold = 1; % Threshold in mm of framewise displacement, at or above which a volume is considered a movement outlier
fontsize = 12;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Initialize variables and data for calcs/plots
f4D_img = spm_read_vols(spm_vol(functional4D_fn));
[Ni, Nj, Nk, Nt] = size(f4D_img);
Ndims = [Ni; Nj; Nk];
slice_number = floor(Nk/2);
[dir, fn, ext] = fileparts(functional4D_fn);
if ~exist([dir filesep fn '_00001' ext],'file')
    f4D_img_spm = spm_file_split(functional4D_fn);
end
rt_data.F_dyn_img = f4D_img(:,:,:,1);
rt_data.tSNR_dyn_img = f4D_img(:,:,:,1);
rt_data.FD_outliers = []; % Outlier volumes
rt_data.FD_sum = 0;
rt_data.outlier_counter = 0; % Outlier counter
rt_data.FD_dynamic_average = 0;
% Nt = ; % specify, or derive from functional_fn
t = 1:Nt;
MP1 = zeros(1,Nt);
MPall = cell(1,6);
for p = 1:6
    MPall{p} =  MP1;
end
FD = zeros(Nt,1);
DVARS = zeros(1,Nt);

title_1 = 'Unprocessed volume #';
title_2 = 'tSNR image';
title_3 = 'HMPs: translation';
title_4 = 'HMPs: rotation';
title_5 = 'FD';
title_6 = 'DVARS';
title_7 = 'THEPLOT';

% Create the GUI figure
fig = figure('Visible','off', 'units','normalized','outerposition',[0 0 1 1]);

% Create UIcontrols
% Position = [left bottom width height]

txt_1 = uicontrol('Style','text',...
    'Position',[20 730 70 40],...
    'String','OPTIONS:', 'fontsize', fontsize);
pb_startPreProc = uicontrol(gcf, 'Style', 'push', 'String', 'PreProc', ...
    'Position', [20 710 70 40],...
    'CallBack', @startPreProc, ...
    'UserData', 0);
pb_startRtSim = uicontrol(gcf, 'Style', 'push', 'String', 'Real-time', ...
    'Position', [20 660 70 40],...
    'CallBack', @startRtSim, ...
    'UserData', 0);
txt_dyn = uicontrol('Style','text',...
    'Position',[100 650 70 40],...
    'String','#0', 'fontsize', fontsize);

txt_2 = uicontrol('Style','text',...
    'Position',[20 590 100 40],...
    'String','REAL-TIME INFO:', 'fontsize', fontsize);
txt_out = uicontrol('Style','text',...
    'Position',[20 560 70 40],...
    'String','Outliers (FD):', 'fontsize', fontsize);
txt_out_val = uicontrol('Style','text',...
    'Position',[100 560 70 40],...
    'String','0', 'fontsize', fontsize);
txt_fdsum = uicontrol('Style','text',...
    'Position',[20 520 70 40],...
    'String','FD sum:', 'fontsize', fontsize);
txt_fdsum_val = uicontrol('Style','text',...
    'Position',[100 520 70 40],...
    'String','0', 'fontsize', fontsize);
txt_fdave = uicontrol('Style','text',...
    'Position',[20 480 70 40],...
    'String','FD average:', 'fontsize', fontsize);
txt_fdave_val = uicontrol('Style','text',...
    'Position',[100 480 70 40],...
    'String','0', 'fontsize', fontsize);


txt_plane = uicontrol('Style','text',...
    'Position',[20 210 100 40],...
    'String','Select plane', 'fontsize', fontsize);
popup = uicontrol('Style', 'popup',...
    'String', {'dim1 (X)','dim2 (Y)','dim3 (Z)'},...
    'Value', 3,...
    'Position', [20 180 100 40],...
    'Callback', @setDim);

%     Add a text uicontrol to label the slice slider.
txt_slice = uicontrol('Style','text',...
    'Position',[20 300 100 40],...
    'String',['Change slice: #' num2str(slice_number)], 'fontsize', fontsize);
% Create slider
sld_slice = uicontrol('Style', 'slider',...
    'Min',1,'Max',Ndims(popup.Value),'Value',slice_number,...
    'Position', [20 270 100 40],...
    'SliderStep', [1/Ndims(popup.Value) 1/Ndims(popup.Value)],...
    'Callback', @changeSlice);







% Create axes
ax1 = subplot(7,4,[17 21 25]); % Unprocessed dynamic image
imagesc(squeeze(f4D_img(:,:,slice_number,1))); title(title_1);  colormap bone; colorbar;

ax2 = subplot(7,4,[18 22 26]); % tSNR image
imagesc(squeeze(f4D_img(:,:,slice_number,1))); title(title_2); colormap bone; colorbar;

ax3 = subplot(7,4,[3 4]); % HMPs - translation
hold on;
title(title_3, 'fontsize', fontsize);
p1 = plot(ax3, t, MPall{1}, 'r');
p2 = plot(ax3, t, MPall{2}, 'b');
p3 = plot(ax3, t, MPall{3}, 'g');
set(ax3, 'Color', 'k');
set(ax3,'Xticklabel',[]);
axis([0 (Nt+5) -1 1]);
ylabel('mm', 'fontsize', fontsize)
grid;
hold off;

ax4 = subplot(7,4,[7 8]); % HMPs - rotation
hold on;
title(title_4, 'fontsize', fontsize);
p4 = plot(ax4, t, MPall{4}, 'r');
p5 = plot(ax4, t, MPall{5}, 'b');
p6 = plot(ax4, t, MPall{6}, 'g');
set(ax4, 'Color', 'k');
set(ax4,'Xticklabel',[]);
axis([0 (Nt+5) -0.02 0.02]);
ylabel('rad', 'fontsize', fontsize)
grid;
hold off;

ax5 = subplot(7,4,[11 12]); % Framewise displacement
p8 = plot(ax5, t, FD, 'c');
set(ax5, 'Color', 'k');
set(ax5,'Xticklabel',[]);
title(ax5, title_5, 'fontsize', fontsize);
ylabel('mm', 'fontsize', fontsize)
axis([0 (Nt+5) 0 2]);
grid;

ax6 = subplot(7,4,[15 16]); % DVARS
p7 = plot(ax6, t, DVARS, 'r');
set(ax6, 'Color', 'k');
set(ax6,'Xticklabel',[]);
title(ax6, title_6, 'fontsize', fontsize);
axis([0 (Nt+5) 0 1000]);
ylabel('a.u.', 'fontsize', fontsize)
grid;

set( findall( fig, '-property', 'Units' ), 'Units', 'Normalized' )
fig.Visible = 'on';


% UIcontrol functions

    % startPreProc
    % This function is called when user selects the 
    function startPreProc(PushButton, EventData)
        currentState = get(PushButton, 'UserData');        
        
        % Preprocess structural and f0 images
        [d, f, e] = fileparts(structural_fn);
        if exist([d filesep 'rc1' f e], 'file')
            preproc_data = struct;
            preproc_data.forward_transformation = [d filesep 'y_' f e];
            preproc_data.inverse_transformation = [d filesep 'iy_' f e];
            preproc_data.gm_fn = [d filesep 'c1' f e];
            preproc_data.wm_fn = [d filesep 'c2' f e];
            preproc_data.csf_fn = [d filesep 'c3' f e];
            preproc_data.bone_fn = [d filesep 'c4' f e];
            preproc_data.soft_fn = [d filesep 'c5' f e];
            preproc_data.air_fn = [d filesep 'c6' f e];
            preproc_data.rstructural_fn = [d filesep 'r' f e];
            preproc_data.rgm_fn = [d filesep 'rc1' f e];
            preproc_data.rwm_fn = [d filesep 'rc2' f e];
            preproc_data.rcsf_fn = [d filesep 'rc3' f e];
            preproc_data.rbone_fn = [d filesep 'rc4' f e];
            preproc_data.rsoft_fn = [d filesep 'rc5' f e];
            preproc_data.rair_fn = [d filesep 'rc6' f e];
        else
            preproc_data = preRtPreProc(functional0_fn, structural_fn, spm_dir);
        end
        
        % Create binary GM, WM and CSF masks
        [GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.1);
        
        % Initialize variables for real-time simulation
        preproc_data.I_GM = find(GM_img_bin);
        preproc_data.I_WM = find(WM_img_bin);
        preproc_data.I_CSF = find(CSF_img_bin);
        preproc_data.mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
        preproc_data.I_mask = find(preproc_data.mask_reshaped);
        preproc_data.N_maskvox = numel(preproc_data.I_mask);
        preproc_data.N_vox = Ni*Nj*Nk;
        preproc_data.line1_pos = numel(preproc_data.I_GM);
        preproc_data.line2_pos = numel(preproc_data.I_GM) + numel(preproc_data.I_WM);
        
        disp('Preprocessing done!')
    end


    function startRtSim(PushButton, EventData)
        currentState = get(PushButton, 'UserData');
        
        % Initialize variables to be updated during each real-time
        % iteration
        rt_data.F_realigned = zeros(preproc_data.N_vox, Nt);
        rt_data.F_smoothed = zeros(preproc_data.N_vox, Nt);
        rt_data.F_dyn = zeros(preproc_data.N_vox, Nt);
        rt_data.F_dyn_img = zeros(Ni, Nj, Nk);
        rt_data.F_dyn_detrended = zeros(preproc_data.N_vox, Nt);
        rt_data.F_perc_signal_change = zeros(preproc_data.N_vox, Nt);
        rt_data.F_theplot = zeros(preproc_data.N_vox, Nt+5);
        rt_data.F_cumulative_mean = zeros(preproc_data.N_vox, Nt);
        rt_data.F_cumulative_sum = zeros(preproc_data.N_vox, Nt);
        rt_data.F_stdev = zeros(preproc_data.N_vox, Nt);
        rt_data.F_tSNR = zeros(preproc_data.N_vox, Nt);
        rt_data.tSNR_dyn_img = zeros(Ni, Nj, Nk);
        rt_data.X_MP_detrending = (1:Nt)';
        rt_data.X_MP_detrending = rt_data.X_MP_detrending - mean(rt_data.X_MP_detrending); % demean non-constant regressors
        rt_data.X_MP_detrending = [rt_data.X_MP_detrending ones(Nt,1)]; % add constant regressor
        rt_data.X_Func_detrending = rt_data.X_MP_detrending;
        rt_data.MP = zeros(Nt,6);
        rt_data.MP_detrended = zeros(Nt,6);
        rt_data.MP_mm = zeros(Nt,6);
        
        % Create axes and image handles for empty THEPLOT
        ax7 = subplot(7,4,[19 20 23 24 27 28]);
        GM_img = rt_data.F_theplot(preproc_data.I_GM, :);
        WM_img = rt_data.F_theplot(preproc_data.I_WM, :);
        CSF_img = rt_data.F_theplot(preproc_data.I_CSF, :);
        all_img = [GM_img; WM_img; CSF_img];
        rt_data.im1 = imagesc(ax7, all_img); colormap(gray); caxis([-5 5]);
        hold on; 
        rt_data.p9 = line([1 Nt],[preproc_data.line1_pos preproc_data.line1_pos],  'Color', 'b', 'LineWidth', 2 );
        rt_data.p10 = line([1 Nt],[preproc_data.line2_pos preproc_data.line2_pos],  'Color', 'r', 'LineWidth', 2 );
        hold off;
        title(ax7, title_7);
        xlabel('Dynamic #')
        ylabel('THEPLOT')        
        set(findall( fig, '-property', 'Units' ), 'Units', 'Normalized' )

        [d, f, e] = fileparts(functional4D_fn);
        T = zeros(Nt,8); % Timing vector: realignment, FD, DVARS, smoothing, detrending, PSC, draw stuff, total
        
        % Start real-time simulation run
        for i = 1:Nt
            t0_start = clock;
            
            txt_dyn.String = ['#' num2str(i)];
            
            % 0: Load dynamic functional image
            fdyn_fn = [d filesep f '_' sprintf('%05d',i) e]; % filename of dynamic functional image
            rt_data.F_dyn_img = spm_read_vols(spm_vol(fdyn_fn)); % this is the unprocessed image to be used for DVARS and THEPLOT
            rt_data.F_dyn(:,i) = rt_data.F_dyn_img(:);
            
            
            % 1:    rtFD - REAL-TIME FRAMEWISE DISPLACEMENT
            % Real-time calculation of estimated absolute displacement of
            % brain, per dynamic. This can help to identify outliers.
            % a) First realign each dynamic image data to functional0 to
            % get Head Movement/Motion Parameters (HMPs / MPs)
            t1_start = clock;
            r_fdyn_fn = rtRealignReslice(functional0_fn, fdyn_fn);
            T(i,1) = etime(clock,t1_start);
            % b) Then load MPs, calculate FD, and load figure data
            t2_start = clock;
            mp = load([d filesep 'rp_' f '.txt']); % stuff todo
            rt_data.MP(i,:) = mp(end,:);
            rt_data.MP_mm(i,:) = rt_data.MP(i,:);
            rt_data.MP_mm(i,4:6) = rt_data.MP(i,4:6)*50; % 50mm = assumed radius of subject head; from Power et al paper (2014) [1].
            if i == 1
                mp_diff = zeros(1, 6); % if first dynamic is realigned to functional0, this is technically not correct. TODO
            else
                mp_diff = diff(rt_data.MP_mm(i-1:i, :));
            end
            fd = sum(abs(mp_diff),2);        
            FD(i) = fd;
            for k = 1:6
                MPall{k}(i) =  rt_data.MP(i,k);
            end
            T(i,2) = etime(clock,t2_start);
            
            
            % 2:    rtDVARS - REAL-TIME ROOT MEAN VARIANCE OF IMAGE INTENSITY DIFFERENCE
            % Real-time calculation of the variance of difference images,
            % which is an indication of how much the image intensity
            % changes between frames and can help to identify outliers.
            % Description from Power et al paper (2014) [1].
            % Method from Chris Rorden's scripts [2].
            t3_start = clock;
            if i == 1
                f_diff = (zeros(1, preproc_data.N_vox))';
            else
                f_diff = diff((rt_data.F_dyn(:, i-1:i))');
                f_diff = f_diff';
            end
            dvars = var(f_diff); % Root Mean Variance across voxels
            DVARS(:, i) = dvars;
            T(i,3) = etime(clock,t3_start);
                 
            
            % 3:    rtTHEPLOT
            % Real-time calculation of datapoints for customized Matlab version of THEPLOT (from Power, 2017 [3]).
            % This heatmap shows percentage signal change of unprocessed
            % BOLD signal alongside quality traces like FD, DVARS and
            % physiological data, and allows visual inspection of data
            % quality issues.
            % a) First smooth dynamic functional data (see paper) using
            % specified Gaussian kernel fwhm.            
            t4_start = clock;
            s_fdyn_fn = rtSmooth(fdyn_fn, [6 6 6]);
            s_fdyn_img = spm_read_vols(spm_vol(s_fdyn_fn));
            rt_data.F_smoothed(:,i) = s_fdyn_img(:);
            T(i,4) = etime(clock,t4_start);            
            % b) Then detrend masked data using cumulative GLM. Mask uses
            % combination of GM, WM and CSF masks as derived in
            % preRtPreProc step. Write last iteration's functional data to
            % static value matrix (but growing in size on each iteration). 
            % Data detrending is necessary for THEPLOT because (reference)
            t5_start = clock;
            x_func_detrending = rt_data.X_Func_detrending(1:i, :);
            beta_func = x_func_detrending\rt_data.F_smoothed(preproc_data.I_mask,1:i)'; % func = X*beta + e ==> beta = X\func ==> func_detrended = mp - X(i)*beta(i)
            F_detrended = rt_data.F_smoothed(preproc_data.I_mask,1:i)' - x_func_detrending(:, 1)*beta_func(1, :); % remove effects of all regressors except constant
            F_detrended = F_detrended';
            rt_data.F_dyn_detrended(preproc_data.I_mask,i) = F_detrended(:,i);
            T(i,5) = etime(clock,t5_start);
            % c) Then calculate percentage signal change for THEPLOT display
            % purposes, using cumulative mean.
            t6_start = clock;
            if i == 1
                rt_data.F_cumulative_mean(preproc_data.I_mask,i) = rt_data.F_dyn_detrended(preproc_data.I_mask,i);
            else
                rt_data.F_cumulative_mean(preproc_data.I_mask,i) = ((i-1)*rt_data.F_cumulative_mean(preproc_data.I_mask,i-1) + rt_data.F_dyn_detrended(preproc_data.I_mask,i))/i;
            end
            
            f_perc_signal_change = 100*(rt_data.F_dyn_detrended(preproc_data.I_mask,i)./rt_data.F_cumulative_mean(preproc_data.I_mask,i)) - 100;
            f_perc_signal_change(isnan(f_perc_signal_change))=0;
            rt_data.F_perc_signal_change(preproc_data.I_mask, i) = f_perc_signal_change;
            rt_data.F_theplot(preproc_data.I_mask, i) = rt_data.F_perc_signal_change(preproc_data.I_mask, i);
            T(i,6) = etime(clock,t6_start);
            
            
            % 4:    STATISTICS AND COUNTERS
            rt_data.F_stdev(preproc_data.I_mask, i) = std(rt_data.F_dyn(preproc_data.I_mask, 1:i), 0, 2);
            rt_data.F_tSNR(preproc_data.I_mask, i) = rt_data.F_cumulative_mean(preproc_data.I_mask,i)./rt_data.F_stdev(preproc_data.I_mask, i);
            rt_data.tSNR_dyn_img = reshape(rt_data.F_tSNR(:, i), Ni, Nj, Nk);
            if fd >= rt_data.FD_threshold
                rt_data.outlier_counter = rt_data.outlier_counter + 1;
                rt_data.FD_outliers = [rt_data.FD_outliers; i];
            end
            rt_data.FD_sum = rt_data.FD_sum + fd;
            rt_data.FD_dynamic_average = rt_data.FD_sum/i;
            txt_out_val.String = num2str(rt_data.outlier_counter);
            txt_fdsum_val.String = num2str(rt_data.FD_sum);
            txt_fdave_val.String = num2str(rt_data.FD_dynamic_average);           
            
            t7_start = clock;
            drawBrains;
            drawMPs;
            drawDVARS;
            drawFD;
            drawTHEPLOT;
            T(i,7) = etime(clock,t7_start);
            T(i,8) = etime(clock,t0_start);
            pause(0.005);
            disp(['iteration done: ' num2str(i)])
        end
       
        assignin('base','rt_data',rt_data);
        assignin('base','preproc_data',preproc_data);
        assignin('base','T',T);
        assignin('base','DVARS',DVARS);
        assignin('base','FD',FD);
    end


    function drawMPs()
        set(p1, 'YData', MPall{1});
        set(p2, 'YData', MPall{2});
        set(p3, 'YData', MPall{3});
        set(p4, 'YData', MPall{4});
        set(p5, 'YData', MPall{5});
        set(p6, 'YData', MPall{6});
        drawnow;
    end


    function drawDVARS()
        set(p7, 'YData', DVARS);
        drawnow;
    end

    function drawFD()
        set(p8, 'YData', FD);
        drawnow;
    end

    function drawTHEPLOT()
        
        GM_img = rt_data.F_theplot(preproc_data.I_GM, :);
        WM_img = rt_data.F_theplot(preproc_data.I_WM, :);
        CSF_img = rt_data.F_theplot(preproc_data.I_CSF, :);
        all_img = [GM_img; WM_img; CSF_img];
        
        set(rt_data.im1, 'CData', all_img); drawnow
    end

    function setDim(source,event)
        
        sld_slice.Max = Ndims(popup.Value);
        slice_number = floor(Ndims(popup.Value)/2);
        sld_slice.Value = slice_number;
        sld_slice.Max = Ndims(popup.Value);
        sld_slice.SliderStep = [1/Ndims(popup.Value) 1/Ndims(popup.Value)];
        txt_slice.String = ['Change slice: #' num2str(slice_number)];
        drawBrains;
    end


    function changeDyn(source,event)
        idyn = round(source.Value);
        txt_dyn.String = ['Change dynamic: #' num2str(idyn)];
        drawBrains;
    end


    function changeSlice(source,event)
        slice_number = round(source.Value);
        slice_number = min(slice_number, Ndims(popup.Value));
        
        txt_slice.String = ['Change slice: #' num2str(slice_number)];
        drawBrains;
    end


    function drawBrains()
        if popup.Value == 1
            imagesc(ax1, squeeze(rt_data.F_dyn_img(slice_number,:,:))); title(ax1, title_1); colormap(ax1, 'bone'); colorbar(ax1);
            imagesc(ax2, squeeze(rt_data.tSNR_dyn_img(slice_number,:,:))); title(ax2, title_2); colormap(ax2, 'bone'); colorbar(ax2);
        elseif popup.Value == 2
            imagesc(ax1, squeeze(rt_data.F_dyn_img(:,slice_number,:))); title(ax1, title_1); colormap(ax1, 'bone'); colorbar(ax1);
            imagesc(ax2, squeeze(rt_data.tSNR_dyn_img(:,slice_number,:))); title(ax2, title_2); colormap(ax2, 'bone'); colorbar(ax2);
        else
            imagesc(ax1, squeeze(rt_data.F_dyn_img(:,:,slice_number))); title(ax1, title_1); colormap(ax1, 'bone'); colorbar(ax1);
            imagesc(ax2, squeeze(rt_data.tSNR_dyn_img(:,:,slice_number))); title(ax2, title_2); colormap(ax2, 'bone'); colorbar(ax2);
        end
        
    end

end


% Sources: 
% [1] - Power et al. 2014. Methods to detect, characterize, and remove
% motion artifact in resting state fMRI. Neuroimage. https://doi:10.1016/j.neuroimage.2013.08.048
% [2] - Chris Rorden's 'nii_qa_moco' script: https://github.com/rordenlab/spmScripts/blob/master/nii_qa_moco.m
% [3] - Power, 2017. A simple but useful way to assess fMRI scan qualities.
% NeuroImage. https://doi.org/10.1016/j.neuroimage.2016.08.009
