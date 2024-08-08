% This code demonstrate our proposed algorithm, applied to the emperical
% data from SPM software.

% Authors -- Haifeng wu (whf5469@gmail.com)
% Copyright (C) 2023 School of Electrical & Information Technology
%                     Yunnan Minzu University

%   Reference:
%       SPM12 Manual: Chapter 35 Dynamic Causal Modeling for fMRI, 2016
%-------------------------------------------------------------------------

clear all;
close all;

%%%% load structure DCM from SPM software %%%%%
data_path = fileparts(mfilename('fullpath'));
% data_path = fullfileata_path,'GLM');
load(fullfile(data_path,'DCM_mod_bwd.mat')); 

%%%% load structure SPM from SPM software %%%%%
load(fullfile(data_path,'SPM.mat'));
DCM.dt_flag = 'L';       % 'L' for large TR, 'S" for small TR

DCM.SPM = SPM; 

% specify the options for the DCM analysis
options                 = DCM.options;
options.m0_coef         = 0.8 ;
options.freq_cutoff     = 60 ;
options.p0              = 0.8;  % single p0 value (for computational efficiency)
options.iter            = 100;
options.filter_str      = 5;
options.restrictInputs  = 1;
options.type            = 'E';
options.search.flag     = 0;
options.search.p0_all   = 0.1 : 0.1 : 0.9 ;
options.search.m0_all   = 0.1 : 0.1 : 0.9 ;
options.search.freq_all = 10  :  10 : 120 ;
options.filter_flag     = 'H' ;

currentTimer = tic;
[output, options] = DCM_Estimate(DCM, options); 
toc(currentTimer)

