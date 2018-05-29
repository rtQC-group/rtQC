function [GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(gm_fn, wm_fn, csf_fn, threshold)
% This function constructs 3D binary images for GM, WM and CSF based on the
% relative value of GM/WM/CSF tissue probability maps per voxel. It assumes
% that the image filename parameters are for images that are in the same
% space (i.e. they match in size, voxel by voxel). If a threshold is
% specified (0<= threshold <=1), the images are first thresholded before
% the binary masks are calculated. For no threshold, specify zero.
%
% INPUT:
% gm_fn     - filename of gray matter image output from SPM segment routine
% wm_fn     - filename of white matter image output from SPM segment routine
% csf_fn    - filename of CSF image output from SPM segment routine
% threshold - value between 0 and 1 for thresholding GM/WM/CSF segment
%             images before calculation of the binary images
% 
% OUTPUT: 
% Three filenames for binary 3D images
%__________________________________________________________________________
% Copyright (C) 2018 Neu3CA.org
% Written by Stephan Heunis

GM_spm = spm_vol(gm_fn);
WM_spm = spm_vol(wm_fn);
CSF_spm = spm_vol(csf_fn);

GM_img = spm_read_vols(GM_spm);
WM_img = spm_read_vols(WM_spm);
CSF_img = spm_read_vols(CSF_spm);

if threshold ~= 0
    GM_img_thresh = GM_img;
    WM_img_thresh = WM_img;
    CSF_img_thresh = CSF_img;
    GM_img_thresh(GM_img < threshold) = 0;
    WM_img_thresh(WM_img < threshold) = 0;
    CSF_img_thresh(CSF_img < threshold) = 0;
    GM_img_bin = (GM_img_thresh >= WM_img_thresh) & (GM_img_thresh >= CSF_img_thresh) & (GM_img_thresh ~= 0);
    WM_img_bin = (WM_img_thresh > GM_img_thresh) & (WM_img_thresh >= CSF_img_thresh) & (WM_img_thresh ~= 0);
    CSF_img_bin = (CSF_img_thresh > GM_img_thresh) & (CSF_img_thresh > WM_img_thresh) & (CSF_img_thresh ~= 0);
else
    GM_img_bin = (GM_img >= WM_img) & (GM_img >= CSF_img) & (GM_img ~= 0);
    WM_img_bin = (WM_img > GM_img) & (WM_img >= CSF_img) & (WM_img ~= 0);
    CSF_img_bin = (CSF_img > GM_img) & (CSF_img > WM_img) & (CSF_img ~= 0);
end


