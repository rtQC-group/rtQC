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
% img_1=uint16(nifti_img_1.img*1000);
% img_2=uint16(nifti_img_2.img * 1000);
img_1=nifti_img_1.img*1000;
img_2=nifti_img_2.img * 1000;

% check if dimensioms are the same
size_img_1 = size(img_1);
size_img_2 = size(img_2);
if (size_img_1(1) ~= size_img_2(1) ||  size_img_1(2) ~= size_img_2(2) || size_img_1(3) ~= size_img_2(3))
    ok = false;
    return;
end


mutIn_slices = zeros(1,size_img_1(3));

for indSlice = 1:size_img_1(3)
    disp(['Slice: ' num2str(indSlice)]);
    slice_img_1 = squeeze(img_1(:,:,indSlice));
    slice_img_2 = squeeze(img_2(:,:,indSlice));
    tic;
    [mutIn_img_1, out1] = MI_GG(slice_img_1, slice_img_1);
    toc;
    tic;
    [mutIn_img_2, out2] = MI_GG(slice_img_2, slice_img_2);
    toc;
    min_mutIn = min(mutIn_img_1, mutIn_img_2);
    tic;
    [M, out3] = MI_GG(slice_img_1, slice_img_2);
    mutIn_slices(indSlice) = M/min_mutIn;
    toc;
end

congruency = mean(mutIn_slices);
return;
end



%% stephan test scripts

fn1 = '/Users/jheunis/Desktop/sample-data/sub-opennft/rstructScan_PSC.nii';
fn2 = '/Users/jheunis/Desktop/sample-data/sub-opennft/template_func.nii';
slice_nr = 20;
im1 = spm_read_vols(spm_vol(fn1));
im2 = spm_read_vols(spm_vol(fn2));
% im1 = uint16(im1*1000);
% im2 = uint16(im2*1000);

slice_im1 = squeeze(im1(:,:, slice_nr));
slice_im2 = squeeze(im2(:,:, slice_nr));
[M, out] = MI_GG(slice_im1, slice_im2);
% 
% nifti_img_1 = load_untouch_nii(fn1);
% nifti_img_2 = load_untouch_nii(fn2);
% % img_1=uint16(nifti_img_1.img*1000);
% % img_2=uint16(nifti_img_2.img*1000);
% img_1=nifti_img_1.img*1000;
% img_2=nifti_img_2.img*1000;
% slice_img_1 = squeeze(img_1(:,:,slice_nr));
% slice_img_2 = squeeze(img_2(:,:,slice_nr));
% 
% figure;
% subplot(2,2,1); imagesc(slice_im1)
% subplot(2,2,2); imagesc(slice_im2)
% subplot(2,2,3); imagesc(slice_img_1)
% subplot(2,2,4); imagesc(slice_img_2)




