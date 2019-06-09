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
% Function rtQC_defaults
%
% Description:
%
% Parameters:
%
% Returns:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function defaults = rtQC_defaults()

defaults = struct;
% ------------------------------- %
% Defaults for app mode on startup
% ------------------------------- %
defaults.demo_mode = 1; % set to 1 if the app should be started in demo mode, 0 if not
defaults.start_clean = 0; % set to 1 if sample data filenames should not be populated
defaults.use_sample_set = 2; % 0 for no sample set; 1 for OpenNFT task data; 2 for HCP resting state data
defaults.demo_visibility = 'on';
if ~defaults.demo_mode
    defaults.start_clean = 1;
    defaults.use_sample_set = 0;
    defaults.demo_visibility = 'off';
    defaults.mode_string = 'Standard Mode';
else
    defaults.mode_string = 'Demo Mode';
end

% ------------------------------- %
% SPM12 directory location
% ------------------------------- %
defaults.spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';

% ------------------------------- %
% Defaults for data and experiment
% ------------------------------- %
% Defaults for two sample datasets are provided, for when app is run in
% demo mode. Users can specify another custom dataset, although logic for
% handling would have to be added. (e.g. if gui_data.sample_set == 3 ...)
% Sample data sets should be downloaded and locations have to be specified
% Sample set 1: data available at
defaults.sample_set{1}.data_dir = '/Users/jheunis/Desktop/sample-data/sub-opennft';
defaults.sample_set{1}.structural_fn = [defaults.sample_set{1}.data_dir filesep 'structScan_PSC.nii'];
defaults.sample_set{1}.functional4D_fn = [defaults.sample_set{1}.data_dir filesep 'fanon-0007-00001-000001-01.nii'];
defaults.sample_set{1}.ROI1_fn = [defaults.sample_set{1}.data_dir filesep 'lROI_1.nii'];
defaults.sample_set{1}.ROI2_fn = [defaults.sample_set{1}.data_dir filesep 'rROI_2.nii'];
defaults.sample_set{1}.im1_dqchecks_fn = '';
defaults.sample_set{1}.im2_dqchecks_fn = '';
% Sample set 2: data available at
defaults.sample_set{2}.data_dir = '/Users/jheunis/Desktop/sample-data/sub-hcp';
defaults.sample_set{2}.structural_fn = [defaults.sample_set{2}.data_dir filesep 'mprage.nii'];
defaults.sample_set{2}.functional4D_fn = [defaults.sample_set{2}.data_dir filesep 'rest.nii'];
defaults.sample_set{2}.ROI1_fn = '';
defaults.sample_set{2}.ROI2_fn = '';
defaults.sample_set{2}.im1_dqchecks_fn = '';
defaults.sample_set{2}.im2_dqchecks_fn = '';
% Sample set 3 (user to add own data information here):
% defaults.sample_set{3}.data_dir = '';
% defaults.sample_set{3}.structural_fn = [defaults.sample_set{3}.data_dir filesep ''];
% defaults.sample_set{3}.functional4D_fn = [defaults.sample_set{3}.data_dir filesep ''];
% defaults.sample_set{3}.ROI1_fn = '';
% defaults.sample_set{3}.ROI2_fn = '';
% defaults.sample_set{3}.im1_dqchecks_fn = '';
% defaults.sample_set{3}.im2_dqchecks_fn = '';
% If app is not run in demo mode (i.e. standard mode) start with empty
% defaults
defaults.empty_set.data_dir = '';
defaults.empty_set.structural_fn = '';
defaults.empty_set.functional4D_fn = '';
defaults.empty_set.ROI1_fn = '';
defaults.empty_set.ROI2_fn = '';
defaults.empty_set.im1_dqchecks_fn = '';
defaults.empty_set.im2_dqchecks_fn = '';
% Experimental parameters
defaults.Ndims = [50 50 50]; % initialising a dimension array, keep as is
defaults.slice_number = 20; % initial slice to be displayed for entered data
defaults.Nt_default = 155; % number of volumes in fMRI timeseries
defaults.FD_threshold = 0.25; % in mm - from Power et al , 2012
defaults.FD_radius = 50; % in mm - from Power et al , 2012
defaults.fwhm = 6; % in mm - initialise to twice the image voxel size

% ------------------------------- %
% Defaults for app design and content
% ------------------------------- %
defaults.img_fn = 'rtqc_logo2_black.png';
defaults.small_font_size = 10;
defaults.standard_font_size = 13;
defaults.medium_font_size = 16;
defaults.button_font_size = 13;
defaults.large_font_size = 24;
defaults.axes_font_size = 12;
defaults.vals_font_size = 16;
defaults.panel_defaults_string = ['Please specify settings necessary for OnlineQC, including directory locations ' ...
    'threshold values and sample data selections. The framewise displacement (FD) threshold is set at the minimum value that ' ...
    'would classify a volume as a motion outlier (typically ranging from 0.2 to 0.5 mm)'];
defaults.panel_preproc_string = ['Before Online QC can commence, structural and functional data have to ' ...
    'reside in the same subject space such that structural masks (e.g. gray matter) can be mapped in real-time. ' ...
    'Please specify pre-real-time collected structural and functional data and run the preprocessing pipeline (note: this can take several minutes).'];
defaults.panel_roi_string = ['If you are interested in tracking quality metrics within specific regions of interest (ROIs), ' ...
    'please load the image volumes (note: these are assumed to be in register with the subject''s pre-collected functional data).'];
defaults.panel_dqchecks_string = ['Real-time and standard fMRI data are sometimes acquired via different procedures (e.g. image reconstructions, ' ... 
    'format transformations etc.). This can result in important differences. This QC step facilitates comparing the output of the real-time ' ...
    'acquisition process with its offline equivalent.'];
defaults.panel_dqchecksA_string = ['Both images (real-time acquired and offline) must be in NIfTI format. The latter should be exported from standard, ' ...
    'offline storage and converted via a standard offline procedure (e.g. dcm2niix, MRIConvert, etc.). If the mutual information comparison does not show ' ...
    'a close match (i.e. a high value), appropriate corrective steps should be taken.'];
defaults.panel_dqchecksB_string = 'In case real-time data are exported in PAR/REC format, we provide a streamlined procedure to rectify header info, orientation and scaling.';

