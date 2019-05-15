%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,USA.
%
% Author: Stephan Heunis, <j.s.heunis@tue.nl>, 2018
%
% Function gui_setup
%
% Description:
%
% Parameters:
%
% Returns:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create tabs
gui_data.tgroup = uitabgroup('Parent', fig,...
    'SelectionChangedFcn', @tabChangedCallback);
gui_data.tab_info = uitab('Parent', gui_data.tgroup, 'Title', 'Info');
gui_data.tab_pre = uitab('Parent', gui_data.tgroup, 'Title', 'Pre QC');
gui_data.tab_online = uitab('Parent', gui_data.tgroup, 'Title', 'Online QC');
gui_data.tab_post = uitab('Parent', gui_data.tgroup, 'Title', 'Post QC');

% Set variable values
gui_data.current_dir = mfilename('fullpath');
gui_data.folders = regexp(gui_data.current_dir, filesep, 'split');
gui_data.root_dir = gui_data.folders{1};
for i = 2:(numel(gui_data.folders)-2)
    gui_data.root_dir = [gui_data.root_dir filesep gui_data.folders{i}];
end
gui_data.data_dir = '/Users/jheunis/Desktop/sample-data/sub-opennft';
% gui_data.data_dir = '/Users/jheunis/Desktop/sample-data/sub-opennft';
if ~defaults.start_clean
    gui_data.structural_fn = [gui_data.data_dir filesep 'structScan_PSC.nii'];
    gui_data.functional4D_fn = [gui_data.data_dir filesep 'fanon-0007-00001-000001-01.nii'];
    gui_data.ROI1_fn = [gui_data.data_dir filesep 'lROI_1.nii'];
    gui_data.ROI2_fn = [gui_data.data_dir filesep 'rROI_2.nii'];
else
    gui_data.structural_fn = '';
    gui_data.functional4D_fn = '';
    gui_data.ROI1_fn = '';
    gui_data.ROI2_fn = '';
end

% -- gui_data.RT_status = {'initialized', 'running', 'stopped', 'completed'}
gui_data.preProc_status = 0; % 0 = not done; 1 = done;
gui_data.preProc_step1_status = 0; % 0 = not done; 1 = done;
gui_data.preProc_step2_status = 0; % 0 = not done; 1 = done;
gui_data.preProc_step3_status = 0; % 0 = not done; 1 = done;

gui_data.ax_volumes_title = 'Realigned volume';
gui_data.ax_FD_title = 'Framewise displacement (FD)';
gui_data.ax_globalZ_title = 'Global Z-score';
gui_data.ax_tSNR_title = 'tSNR';
gui_data.ax_thePlot_title = 'The Plot';

% --- Create UIcontrols for INFO tab --- 
gui_data.panel_main = uipanel('Parent', gui_data.tab_info,...
    'Title','',...
    'Position',[.25 .15 .5 .7],...
    'fontsize', gui_data.large_font_size);
gui_data.ax1 = axes('Parent', gui_data.panel_main,...
    'Position',[.08 .35 .84 .55]);
gui_data.txt_rtQcsub = uicontrol('Parent', gui_data.panel_main,...
    'Style','text',...
    'String','An open-source collaborative framework for quality control methods in real-time fMRI',...
    'Units', 'Normalized',...
    'Position',[0.02 0.2 0.96 0.05],...
    'fontangle', 'italic',...
    'fontsize', gui_data.medium_font_size);
gui_data.pb_help = uicontrol('Parent', gui_data.panel_main,...
    'Style', 'push',...
    'String', 'Github', ...
    'Units', 'Normalized',...
    'Position', [0.1 0.02 0.2 0.1],...
    'CallBack', @linkToGithub, ...
    'UserData', 0,...
    'fontweight', 'bold',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_docs = uicontrol('Parent', gui_data.panel_main,...
    'Style', 'push',...
    'String', 'Readthedocs', ...
    'Units', 'Normalized',...
    'Position', [0.4 0.02 0.2 0.1],...
    'CallBack', @linkToDocs, ...
    'UserData', 0,...
    'fontweight', 'bold',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_start = uicontrol('Parent', gui_data.panel_main,...
    'Style', 'push',...
    'String', 'Start', ...
    'Units', 'Normalized',...
    'Position', [0.7 0.02 0.2 0.1],...
    'CallBack', @gotoTab1, ...
    'UserData', 0,...
    'fontweight', 'bold',...
    'fontsize', gui_data.button_font_size);
brain_img = imread(which(gui_data.img_fn));
image(gui_data.ax1, brain_img)
removeAxesTicks(gui_data.ax1);



%  --- Create UIcontrols for PRE QC tab --- 
% Create panels
gui_data.panel_defaults = uipanel('Parent', gui_data.tab_pre,...
    'Title','Online QC defaults',...
    'Position',[.02 .73 .48 .25],...
    'fontweight', 'bold',...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_preproc = uipanel('Parent', gui_data.tab_pre,...
    'Title','Data and preprocessing',...
    'Position',[.02 .3 .48 .4],...
    'fontweight', 'bold',...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_roi = uipanel('Parent', gui_data.tab_pre,...
    'Title','Region(s) of interest',...
    'Position',[.02 .02 .48 .25],...
    'fontweight', 'bold',...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_other = uipanel('Parent', gui_data.tab_pre,...
    'Title','Quality checks...',...
    'Position',[.52 .6 .45 .38],...
    'fontweight', 'bold',...
    'fontsize', gui_data.standard_font_size);



% Defaults panel UIControls 
gui_data.txt_panel_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style','text',...
    'String',defaults.panel_defaults_string,...
    'Units', 'Normalized',...
    'Position',[0.02 0.42 0.96 0.56],...
    'HorizontalAlignment', 'left',...
    'fontangle', 'italic',...
    'fontsize', gui_data.axes_font_size);
gui_data.txt_spm_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style','text',...
    'String','SPM12 directory',...
    'Units', 'Normalized',...
    'Position',[0.02 0.25 0.15 0.14],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_spm_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style', 'edit',...
    'String', gui_data.spm_dir, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.25 0.73 0.14],...
    'CallBack', @editSPMdir, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_spm_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.25 0.06 0.14],...
    'CallBack', @setSPMdir, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_fd_threshold_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style','text',...
    'String','FD Threshold',...
    'Units', 'Normalized',...
    'Position',[0.02 0.04 0.15 0.14],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_fd_threshold_defaults = uicontrol('Parent', gui_data.panel_defaults,...
    'Style', 'edit',...
    'String', num2str(gui_data.FD_threshold), ...
    'Units', 'Normalized',...
    'Position', [0.18 0.05 0.1 0.14],...
    'CallBack', @editFDthreshold, ...
    'HorizontalAlignment', 'right',...
    'fontsize', gui_data.button_font_size);



% Preprocessing panel UIControls
gui_data.txt_panel_preproc = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String',defaults.panel_preproc_string,...
    'Units', 'Normalized',...
    'Position',[0.02 0.7 0.96 0.28],...
    'HorizontalAlignment', 'left',...
    'fontangle', 'italic',...
    'fontsize', gui_data.axes_font_size);
gui_data.txt_structural_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','Structural .nii',...
    'Units', 'Normalized',...
    'Position',[0.02 0.58 0.15 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_structural_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style', 'edit',...
    'String', gui_data.structural_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.58 0.73 0.10],...
    'CallBack', @editStructural, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_structural_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.58 0.06 0.10],...
    'CallBack', @setStructural, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_functional_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','Functional .nii',...
    'Units', 'Normalized',...
    'Position',[0.02 0.42 0.15 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_functional_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style', 'edit',...
    'String', gui_data.functional4D_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.42 0.73 0.10],...
    'CallBack', @editFunctional, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_functional_pre = uicontrol('Parent', gui_data.panel_preproc,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.42 0.06 0.10],...
    'CallBack', @setFunctional, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.pb_preproc = uicontrol('Parent', gui_data.panel_preproc,...
    'Style', 'push',...
    'String', 'RUN PREPROC', ...
    'Units', 'Normalized',...
    'Position', [0.03 0.1 0.2 0.2],...
    'CallBack', @runPreProc, ...
    'UserData', 0,...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.txt_preproc1 = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','A. Coregister anatomical to functional space:',...
    'Units', 'Normalized',...
    'Position',[0.32 0.26 0.5 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.txt_preproc2 = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','B. Segment coregistered anatomical image:',...
    'Units', 'Normalized',...
    'Position',[0.32 0.14 0.5 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.txt_preproc3 = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','C. Reslice all images to functional space:',...
    'Units', 'Normalized',...
    'Position',[0.32 0.02 0.5 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.txt_preproc1_status = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','-',...
    'Units', 'Normalized',...
    'Position',[0.8 0.28 0.1 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_preproc2_status = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','-',...
    'Units', 'Normalized',...
    'Position',[0.8 0.16 0.1 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_preproc3_status = uicontrol('Parent', gui_data.panel_preproc,...
    'Style','text',...
    'String','-',...
    'Units', 'Normalized',...
    'Position',[0.8 0.04 0.1 0.10],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);


% ROI panel UIControls
% [0.02 0.8 0.15 0.14] - [0.18 0.8 0.73 0.14] - [0.92 0.8 0.06 0.14]
% [0.02 0.6 0.15 0.14] - [0.18 0.6 0.73 0.14] - [0.92 0.6 0.06 0.14]
gui_data.txt_panel_roi = uicontrol('Parent', gui_data.panel_roi,...
    'Style','text',...
    'String',defaults.panel_roi_string,...
    'Units', 'Normalized',...
    'Position',[0.02 0.45 0.96 0.53],...
    'HorizontalAlignment', 'left',...
    'fontangle', 'italic',...
    'fontsize', gui_data.axes_font_size);
gui_data.txt_roi1_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style','text',...
    'String','ROI 1 (.nii)',...
    'Units', 'Normalized',...
    'Position',[0.02 0.27 0.15 0.14],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_roi1_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style', 'edit',...
    'String', gui_data.ROI1_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.28 0.73 0.14],...
    'CallBack', @editROI1, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_roi1_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.28 0.06 0.14],...
    'CallBack', @setROI1, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_roi2_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style','text',...
    'String','ROI 2 (.nii)',...
    'Units', 'Normalized',...
    'Position',[0.02 0.07 0.15 0.14],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_roi2_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style', 'edit',...
    'String', gui_data.ROI2_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.08 0.73 0.14],...
    'CallBack', @editROI2, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_roi2_pre = uicontrol('Parent', gui_data.panel_roi,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.08 0.06 0.14],...
    'CallBack', @setROI2, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);



% --- Create UIcontrols for ONLINE QC tab ---
% Create panels
gui_data.panel_settings = uipanel('Parent', gui_data.tab_online,...
    'Title','Settings',...
    'Position',[.02 .7 .36 .28],...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_controls = uipanel('Parent', gui_data.tab_online,...
    'Title','Controls',...
    'Position',[.02 .58 .36 .1],...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_volumes = uipanel('Parent', gui_data.tab_online,...
    'Title','Real-time Volume Display',...
    'Position',[.02 .02 .36 .53],...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_rtvalues = uipanel('Parent', gui_data.tab_online,...
    'Title','Real-time Metrics',...
    'Position',[.4 .8 .58 .18],...
    'fontsize', gui_data.standard_font_size);
gui_data.panel_rtplots = uipanel('Parent', gui_data.tab_online,...
    'Title','Real-time QC Display',...
    'Position',[.4 .02 .58 .76],...
    'fontsize', gui_data.standard_font_size);

% Settings panel UIControls
% [0.02 0.82 0.15 0.12] - [0.18 0.82 0.73 0.12] - [0.92 0.82 0.06 0.12]
% [0.02 0.66 0.15 0.12] - [0.18 0.66 0.73 0.12] - [0.92 0.66 0.06 0.12]
% [0.02 0.5 0.15 0.12] - [0.18 0.5 0.1 0.12]
% [0.02 0.34 0.15 0.12] - [0.18 0.34 0.73 0.12] - [0.92 0.34 0.06 0.12]
% [0.02 0.18 0.15 0.12] - [0.18 0.18 0.73 0.12] - [0.92 0.18 0.06 0.12]
gui_data.txt_structural = uicontrol('Parent', gui_data.panel_settings,...
    'Style','text',...
    'String','Structural .nii',...
    'Units', 'Normalized',...
    'Position',[0.02 0.82 0.15 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_structural = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'edit',...
    'String', gui_data.structural_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.82 0.73 0.12],...
    'CallBack', @editStructural, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_structural = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.82 0.06 0.12],...
    'CallBack', @setStructural, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_functional = uicontrol('Parent', gui_data.panel_settings,...
    'Style','text',...
    'String','Functional .nii',...
    'Units', 'Normalized',...
    'Position',[0.02 0.66 0.15 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_functional = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'edit',...
    'String', gui_data.functional4D_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.66 0.73 0.12],...
    'CallBack', @editFunctional, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_functional = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.66 0.06 0.12],...
    'CallBack', @setFunctional, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_fd_threshold = uicontrol('Parent', gui_data.panel_settings,...
    'Style','text',...
    'String','FD Threshold',...
    'Units', 'Normalized',...
    'Position',[0.02 0.5 0.15 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_fd_threshold = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'edit',...
    'String', num2str(gui_data.FD_threshold), ...
    'Units', 'Normalized',...
    'Position', [0.18 0.5 0.1 0.12],...
    'CallBack', @editFDthreshold, ...
    'HorizontalAlignment', 'right',...
    'fontsize', gui_data.button_font_size);
gui_data.txt_roi1 = uicontrol('Parent', gui_data.panel_settings,...
    'Style','text',...
    'String','ROI 1 (.nii)',...
    'Units', 'Normalized',...
    'Position',[0.02 0.34 0.15 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_roi1 = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'edit',...
    'String', gui_data.ROI1_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.34 0.73 0.12],...
    'CallBack', @editROI1, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_roi2 = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.34 0.06 0.12],...
    'CallBack', @setROI1, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.txt_roi2 = uicontrol('Parent', gui_data.panel_settings,...
    'Style','text',...
    'String','ROI 2 (.nii)',...
    'Units', 'Normalized',...
    'Position',[0.02 0.18 0.15 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.small_font_size);
gui_data.edit_roi2 = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'edit',...
    'String', gui_data.ROI2_fn, ...
    'Units', 'Normalized',...
    'Position', [0.18 0.18 0.73 0.12],...
    'CallBack', @editROI2, ...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.button_font_size);
gui_data.pb_roi2 = uicontrol('Parent', gui_data.panel_settings,...
    'Style', 'push',...
    'String', '...', ...
    'Units', 'Normalized',...
    'Position', [0.92 0.18 0.06 0.12],...
    'CallBack', @setROI2, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);


% RT controls panel UIcontrols
gui_data.pb_initialize = uicontrol('Parent', gui_data.panel_controls,...
    'Style', 'push',...
    'String', 'INITIALIZE', ...
    'Units', 'Normalized',...
    'Position', [0.05 0.3 0.25 0.45],...
    'CallBack', @initialize, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.pb_startRT = uicontrol('Parent', gui_data.panel_controls,...
    'Style', 'push',...
    'String', 'START', ...
    'Units', 'Normalized',...
    'Position', [0.35 0.3 0.25 0.45],...
    'CallBack', @startRT, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.pb_stopRT = uicontrol('Parent', gui_data.panel_controls,...
    'Style', 'push',...
    'String', 'STOP', ...
    'Units', 'Normalized',...
    'Position', [0.65 0.3 0.25 0.45],...
    'CallBack', @stopRT, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);


% RT values panel UIcontrols
gui_data.txt_Nvolume = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','Acquired volumes:  0',...
    'Units', 'Normalized',...
    'Position',[0.05 0.65 0.4 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_Nvolume_valid = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','Valid volumes:  0',...
    'Units', 'Normalized',...
    'Position',[0.05 0.35 0.4 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_fdsum = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','Total FD:  0',...
    'Units', 'Normalized',...
    'Position',[0.4 0.65 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_fdave = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','Mean FD:  0',...
    'Units', 'Normalized',...
    'Position',[0.4 0.35 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);

gui_data.txt_tsnr = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','tSNR (brain):  0',...
    'Units', 'Normalized',...
    'Position',[0.7 0.78 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_gm = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','tSNR (GM):  0',...
    'Units', 'Normalized',...
    'Position',[0.7 0.55 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_wm = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','tSNR (WM):  0',...
    'Units', 'Normalized',...
    'Position',[0.7 0.32 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_csf = uicontrol('Parent', gui_data.panel_rtvalues,...
    'Style','text',...
    'String','tSNR (CSF):  0',...
    'Units', 'Normalized',...
    'Position',[0.7 0.09 0.3 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);



% RT volume panel UIcontrols
gui_data.txt_plane = uicontrol('Parent', gui_data.panel_volumes,...
    'Style','text',...
    'String','Change plane',...
    'Units', 'Normalized',...
    'Position',[0.08 0.07 0.3 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.standard_font_size);
gui_data.popup_setDim = uicontrol('Parent', gui_data.panel_volumes,...
    'Style', 'popup',...
    'String', {'dim1 (X)','dim2 (Y)','dim3 (Z)'},...
    'Units', 'Normalized',...
    'Position', [0.06 0 0.22 0.1],...
    'Value', 3,...
    'Callback', @setDim);
gui_data.txt_slice = uicontrol('Parent', gui_data.panel_volumes,...
    'Style','text',...
    'String',['Change slice: #' num2str(gui_data.slice_number)],...
    'Units', 'Normalized',...
    'Position',[0.38 0.07 0.3 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.standard_font_size);
gui_data.sld_slice = uicontrol('Parent', gui_data.panel_volumes,...
    'Style', 'slider',...
    'Units', 'Normalized',...
    'Position', [0.36 0 0.3 0.1],...
    'Min',1,...
    'Max',gui_data.Ndims(gui_data.popup_setDim.Value),...
    'Value',gui_data.slice_number,...
    'SliderStep', [1/gui_data.Ndims(gui_data.popup_setDim.Value) 1/gui_data.Ndims(gui_data.popup_setDim.Value)],...
    'Callback', @changeSlice);
gui_data.txt_img = uicontrol('Parent', gui_data.panel_volumes,...
    'Style','text',...
    'String','Change image',...
    'Units', 'Normalized',...
    'Position',[0.73 0.07 0.3 0.12],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.standard_font_size);
gui_data.popup_setImg = uicontrol('Parent', gui_data.panel_volumes,...
    'Style', 'popup',...
    'String', {'fMRI','tSNR'},...
    'Units', 'Normalized',...
    'Position', [0.73 0 0.2 0.1],...
    'Value', 1,...
    'Callback', @changeImg);

% RT plots panel UIcontrols and axes
gui_data.ax_volumesX = axes('Parent', gui_data.panel_volumes,...
    'Visible', 'off',...
    'Position',[.14 .23 .7 .7],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size);
% gui_data.ax_volumesX.Title.String = gui_data.ax_volumes_title;
gui_data.ax_volumesY = axes('Parent', gui_data.panel_volumes,...
    'Visible', 'off',...
    'Position',[.14 .23 .7 .7],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size);
% gui_data.ax_volumesY.Title.String = gui_data.ax_volumes_title;
gui_data.ax_volumesZ = axes('Parent', gui_data.panel_volumes,...
    'Visible', 'on',...
    'Position',[.14 .23 .7 .7],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size);
% gui_data.ax_volumesZ.Title.String = gui_data.ax_volumes_title;

gui_data.ax_FD = axes('Parent', gui_data.panel_rtplots,...
    'Position',[.1 .85 .8 .12],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size ,...
    'XGrid', 'on',...
    'YGrid', 'on');
gui_data.ax_FD.Title.String = gui_data.ax_FD_title;
gui_data.ax_FD.YLabel.String = 'mm';

gui_data.ax_globalZ = axes('Parent', gui_data.panel_rtplots,...
    'Position',[.1 .69 .8 .12],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size ,...
    'XGrid', 'on',...
    'YGrid', 'on');
gui_data.ax_globalZ.Title.String = gui_data.ax_globalZ_title;
gui_data.ax_globalZ.YLabel.String = 'a.u.';

gui_data.ax_tSNR = axes('Parent', gui_data.panel_rtplots,...
    'Position',[.1 .53 .8 .12],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size ,...
    'XGrid', 'on',...
    'YGrid', 'on',...
    'Title',gui_data.ax_tSNR_title);
gui_data.ax_tSNR.Title.String = gui_data.ax_tSNR_title;
gui_data.ax_tSNR.YLabel.String = 'a.u.';

gui_data.ax_thePlot = axes('Parent', gui_data.panel_rtplots,...
    'Position',[.1 .04 .8 .45],...
    'Color', 'k',...
    'FontSize', gui_data.axes_font_size);
gui_data.ax_thePlot.Title.String = gui_data.ax_thePlot_title;
gui_data.ax_thePlot.YLabel.String = 'voxels';
gui_data.ax_thePlot.XLabel.String = 'Volume #';
removeAxesXTickLabels(gui_data.ax_FD);
removeAxesXTickLabels(gui_data.ax_volumesZ);
removeAxesXTickLabels(gui_data.ax_globalZ);
removeAxesXTickLabels(gui_data.ax_tSNR);



% --- Create UIcontrols for POST-NF QC tab ---
% Create panels
gui_data.panel_rtsummary = uipanel('Parent', gui_data.tab_post,...
    'Title','Summary of Real-time Metrics',...
    'Position',[.02 .8 .96 .18],...
    'fontsize', gui_data.standard_font_size);

% Summary panel UIControls
gui_data.txt_Nvolume_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','Acquired volumes:  -',...
    'Units', 'Normalized',...
    'Position',[0.02 0.65 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_Nvolume_valid_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','Valid volumes:  -',...
    'Units', 'Normalized',...
    'Position',[0.02 0.35 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_fdsum_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','Total FD:  -',...
    'Units', 'Normalized',...
    'Position',[0.25 0.65 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_fdave_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','Mean FD:  -',...
    'Units', 'Normalized',...
    'Position',[0.25 0.35 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','tSNR (brain):  -',...
    'Units', 'Normalized',...
    'Position',[0.5 0.78 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_gm_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','tSNR (GM):  -',...
    'Units', 'Normalized',...
    'Position',[0.5 0.55 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_wm_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','tSNR (WM):  -',...
    'Units', 'Normalized',...
    'Position',[0.5 0.32 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.txt_tsnr_csf_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style','text',...
    'String','tSNR (CSF):  -',...
    'Units', 'Normalized',...
    'Position',[0.5 0.09 0.2 0.2],...
    'HorizontalAlignment', 'left',...
    'fontsize', gui_data.vals_font_size);
gui_data.pb_show_fdoutliers_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style', 'push',...
    'String', 'Show FD outliers', ...
    'Units', 'Normalized',...
    'Position', [0.8 0.6 0.18 0.3],...
    'CallBack', @showFDoutliers, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);
gui_data.pb_run_qc_post = uicontrol('Parent', gui_data.panel_rtsummary,...
    'Style', 'push',...
    'String', 'Run Offline QC', ...
    'Units', 'Normalized',...
    'Position', [0.8 0.1 0.18 0.3],...
    'CallBack', @runOfflineQC, ...
    'UserData', 0,...
    'fontsize', gui_data.button_font_size);




% Save GUI data
guidata(fig, gui_data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions

function removeAxesTicks(ax)
set(ax,'xtick',[])
set(ax,'xticklabel',[])
set(ax,'ytick',[])
set(ax,'yticklabel',[])
set(ax,'ztick',[])
set(ax,'zticklabel',[])
end

function removeAxesXTickLabels(ax)
set(ax,'xticklabel',[])
end
