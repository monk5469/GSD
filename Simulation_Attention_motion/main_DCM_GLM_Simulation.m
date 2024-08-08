% This code demonstrate our proposed algorithm, applied to the simulated data.

% Authors -- Haifeng wu (whf5469@gmail.com)
% Copyright (C) 2023 School of Electrical & Information Technology
%               Yunnan Minzu University
%----------------------------------------------------------------------------

clear all;
close all;

%%%% load structure DCM about the simulation.
data_path = fileparts(mfilename('fullpath'));
load(fullfile(data_path,'DCM_Sim.mat')); 
DCM.dt_flag = 'S';       % 'L' for large TR, 'S" for small TR

% specify the options for the DCM analysis
options                 = DCM.options;
options.m0_coef         = 0.8;  %initial value
options.freq_cutoff     = 75;   % cutoff freqeuncy for a low-pass filter
options.p0              = 0.15; % single p0 value (for computational efficiency)
options.iter            = 100;  % the number of iterations
options.filter_str      = 5;    % a coefficient for a filter
options.type            = 'S';  % 'S' for simulation data, 'E' for empirical data 
options.search.flag     = 0;    % 1 for seaching the best value, 0 for not searching
options.search.p0_all   = 0.1 : 0.1 : 0.9 ; % the search range for p0
options.search.m0_all   = 0.1 : 0.1 : 0.9 ; % the search range for m0
options.search.freq_all = 20  :  10 : 140 ; % the search range for cutoff frequency
options.filter_flag     = 'L' ; % 'L' for low-pass filter, 'H' for high-pass filter 
options.filtu           = 1;

currentTimer = tic;
[output, options] = DCM_Estimate(DCM, options); 
toc(currentTimer)