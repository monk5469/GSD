% This code demonstrate our proposed algorithm, applied to the emperical
% data from Henson‘s Face data.
% Data processing tools：SPM12.
% The URL for obtaining the face data：OpenfMRI ( https://www.openfmri.org/dataset/ds000117/ ; Wakeman & Henson, 2015 ).

% Authors -- Haifeng wu (whf5469@gmail.com)
% Copyright (C) 2023 School of Electrical & Information Technology
%                     Yunnan Minzu University

clear all;
close all;

subs  = [1 2 3 5 6 8 9 10 12 14  16 17 18 19 23 24];

for s = 1 : length(subs)
    
clearvars -except  subs output_subs options_subs s

data_path         = fileparts(mfilename('fullpath'));
folder_path       = fullfile(data_path, 'Subjects');
file_dcm          = 'DCM_subject%02d.mat';
file_spm          = 'SPM%02d.mat';
file_path_dcm     = fullfile(folder_path, sprintf(file_dcm, subs(s)));
file_path_spm     = fullfile(folder_path, sprintf(file_spm, subs(s)));
% data_path = fullfileata_path,'GLM');

%%%% load structure DCM from SPM software %%%%%
load(file_path_dcm); 

%%%% load structure SPM from SPM software %%%%%
load(file_path_spm);
DCM.dt_flag    = 'L';       % 'L' for large TR, 'S" for small TR
DCM.SPM        = SPM;

% specify the options for the DCM analysis
options                 = DCM.options;
options.m0_coef         = 0.1;   
options.freq_cutoff     = 400;   
options.p0              = 0.85;  % single p0 value (for computational efficiency)  0.7
options.iter            = 100;
options.filter_str      = 5;
options.restrictInputs  = 1;
options.type            = 'E';
options.search.flag     = 0;
options.search.p0_all   = 0.1 : 0.1 : 0.9 ;
options.search.m0_all   = 0.1 : 0.1 : 0.9 ;
options.search.freq_all = 10  :  10 : 120 ;
options.filter_flag     = 'H' ;
options.connection_model = 1;  %% select different connection model Settings （1 2 3）

%%%%%%%%select different connection model Settings%%%%
[DCM] = DCM_Connection_Model(DCM,options);

currentTimer = tic;
[output, options] = DCM_Estimate(DCM, options); 
toc(currentTimer)                            

output_subs(s)  = output;
options_subs(s) = options;

end





