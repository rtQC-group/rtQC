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

defaults.spm_dir = '/Users/jheunis/Documents/MATLAB/spm12';
defaults.img_fn = 'rtqc_logo2_black.png';
defaults.structural_fn = '/Users/jheunis/Documents/MATLAB/rtQC/rtqc_jsh/0051210/anat_1/mprage.nii';
defaults.functional4D_fn = '/Users/jheunis/Documents/MATLAB/rtQC/rtqc_jsh/0051210/rest_1/rest.nii';
defaults.small_font_size = 10;
defaults.standard_font_size = 13;
defaults.button_font_size = 13;
defaults.large_font_size = 24;
defaults.axes_font_size = 12;
defaults.vals_font_size = 16;
defaults.Ndims = [50 50 50];
defaults.slice_number = 20;

defaults.FD_threshold = 0.25; % in mm
defaults.FD_radius = 50; % in mm
defaults.fwhm = 6; % in mm

