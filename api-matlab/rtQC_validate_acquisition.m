%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
% Author: Lydia Hellrung, <lydia.hellrung@econ.uzh.ch>, 2017
%
% Function rtQC_validate_acquisition
%
% Description:
%   Performs a validation of the real-time acquisition obtained using our standard audiovisual and motor paradigm
%   (see section 4). This will be based on rtQC_validate_analysis where different data will be subjected to the same analysis.
%
% Parameters:
%
% Returns:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ok, congruency] = rtQC_validate_acquisition(path_to_nii_1, path_to_nii_2, options)


ok = true;
congruency = 0;


nifti_img_1 = load_untouch_nii(path_to_nii_1);
nifti_img_2 = load_untouch_nii(path_to_nii_2);
img_1=uint16(nifti_img_1.img*1000);
img_2=uint16(nifti_img_2.img * 1000);

% check if dimensioms are the same
size_img_1 = size(img_1);
size_img_2 = size(img_2);
if (size_img_1(1) ~= size_img_2(1) ||  size_img_1(2) ~= size_img_2(2) || size_img_1(3) ~= size_img_2(3))
    ok = false;
    return;
end


mutIn_slices = zeros(1,size_img_1(3));

for indSlice = 1:size_img_1(3)
    
    slice_img_1 = img_1(:,:,indSlice);
    slice_img_2 = img_2(:,:,indSlice);
    mutIn_img_1 = MI_GG(slice_img_1, slice_img_1);
    mutIn_img_2 = MI_GG(slice_img_2, slice_img_2);
    min_mutIn = min(mutIn_img_1, mutIn_img_2);
    
    mutIn_slices(indSlice) = MI_GG(slice_img_1, slice_img_2)/min_mutIn;
end

congruency = mean(mutIn_slices);
return;
end