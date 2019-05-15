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
defaults.start_clean = 0; % set to 1 if sample data filenames should not be populated
defaults.sample_set = 2; % 0 for no sample set; 1 for OpenNFT task data; 2 for HCP resting state data
defaults.use_WMSCF = 0;
defaults.spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
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
defaults.panel_qchecks_string = '';
defaults.Ndims = [50 50 50];
defaults.slice_number = 20;
defaults.Nt_default = 155; % from sample data
defaults.FD_threshold = 0.25; % in mm
defaults.FD_radius = 50; % in mm
defaults.fwhm = 6; % in mm

