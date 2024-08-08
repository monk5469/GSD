%%% This file will extract Beta of an ROI from a structure SPM, and it need
%%% to added in a file, spm12/spm_regions.m.

[~,Loc] = spm_XYZreg('NearestXYZ',xY.xyz,xSPM.XYZmm);
Beta  = spm_data_read(SPM.Vbeta,'xyz',xSPM.XYZ(:,Loc));
% Find the nearest XYZMM coordination from the input cordination of an ROI. 
% Convert XYZMM to XYZ, and form the converted coordination extract the Beta.
% the two-line codes need to be dded in Line 250 of spm_regions.m


xY.Beta    = Beta;
% Set Beta in a struction xY, the code also needs to be added to spm_regions.m. 




