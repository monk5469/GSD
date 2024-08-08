%---------The  author of the original code-----------------
%%Prof Richard (Rik) Henson, MA, MSc, PhD, FBA
%%Professor of Cognitive Neuroscience, Department of Psychiatry, University of Cambridge
%%Director, Cambridge Centre for Ageing and Neuroscience (CamCAN)
%%President-Past, British Neuroscience Association
%%Fellow of the British Academy

%%MRC Cognition and Brain Sciences Unit
%%University of Cambridge
%%15 Chaucer Road
%%Cambridge, CB2 7EF
%%England
%%URL: http://www.mrc-cbu.cam.ac.uk/people/rik.henson/personal
%%Source code URL：https://github.com/SMScottLee/Face_DCM_fMRI
%----------------------------------------------------------------
%This version retains only the initial parameter setting part of DCM (only the two ROI network parts are preserved)（2023）.
%---------------------------------------------------------------
%%
addpath F:\Face_DCM_fMRI-main\imaging\henson\users\spm12_updates_r7771; 
addpath F:\Face_DCM_fMRI-main\imaging\henson\users\Rik_Updates; % for updated version of spm_regions that handles NaN voxels, and for version of spm_dcm_peb_fit that allows parfor
addpath /imaging/henson/users/sl01/FacesData/Scripts;

%% General
rwd      = 'F:\Face_DCM_fMRI-main/imaging/henson/users/sl01/FacesData/';
owd      = 'F:\Face_DCM_fMRI-main/imaging/henson/users/sl01/FacesData/';
dosubs  = [1 2 3 5 6 8 9 10 12 14  16 17 18 19 23 24]; 
stat_dir = 'stat_concat_deb';
Nses = 1;

%% DCM
ROInames = {'rEVC';'rOFA';'rFFA';'lEVC';'lOFA';'lFFA'};
VOInames = {};
for v = 1:length(ROInames)
    for ses = 1:Nses
        VOInames{v,ses} = ['VOI_' ROInames{v} sprintf('_%d',ses)];
    end
end

scaleU  = 0; % Whether to scale (Z-score) B inputs - suggestion of Peter
% Only contrasts that show sig effect (2-tailed) below
contrasts = [1 1 1  1 1 1  1 1 1;
             1 1 1  1 1 1  0 0 0; % faces
             0 1 0  0 1 0  0 1 0; % immediate repetition
             0 0 1  0 0 1  0 0 1; % delayed repetition            
             1 1 1  0 0 0  0 0 0]; % fame


Ncon = size(contrasts,1);

%% 2 ROI - right OFA & FFA
VOIs = [2:3], Nvoi = length(VOIs), dirname = sprintf('Right%d',Nvoi)
% VOIs = [5:6], Nvoi = length(VOIs), dirname = sprintf('Left%d',Nvoi)
Model(1).a = ones(Nvoi); % fully connected
Model(1).b(:,:,1) = zeros(Nvoi); % no modulation for the first input
for n = 2:Ncon
    Model(1).b(:,:,n) = ones(Nvoi); % modulations for the four remaining inputs/contrasts of interest
end
Model(1).c = zeros(Nvoi,Ncon);
Model(1).c(:,1) = 1; % inputs to OFA, FFA


%% Collate DCMs into GCM
GCM = {}; nsubses = 0;
for s=1:length(dosubs)
    
    sub = dosubs(s);
    swd = sprintf('subject_%02d',sub);
    data_path = fullfile(rwd,swd,'/bold/',stat_dir);
    out_path  = fullfile(owd,swd,'DCM');
    try mkdir(out_path), end
    out_path  = fullfile(out_path,dirname);
    try mkdir(out_path), end
    
    load(fullfile(data_path,'SPM.mat'));
    
    for ses=1:Nses
        
        for v=1:Nvoi
            load(fullfile(owd,swd,stat_dir,VOInames{VOIs(v),ses}),'xY');
            DCM_allM.xY(ses,v) = xY;
        end
        
        nsubses = nsubses+1;
 
        DCM = [];
        
        for v=1:Nvoi
            DCM.xY(v) = DCM_allM.xY(ses,v);
        end;
        
        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points
        
        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        
        DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v); % Models autocorrelation
        
        DCM.U.dt   = SPM.Sess(ses).U(1).dt;
        
        for ui = 1:size(contrasts,1)
            ind = find(contrasts(ui,:)~=0);
            
            tmp = sparse(length(SPM.Sess(ses).U(1).u)-32,1);
            for c = 1:length(ind);
                tmp = tmp + contrasts(ui,ind(c)) * SPM.Sess(ses).U(ind(c)).u(33:end);
            end
            
            if scaleU & ui>1 % Don't scale first modulation (ui=1) because in C? (Peter)
                tmp = tmp/std(tmp);
            end
            
            DCM.U.u(:,ui) = tmp;
            
            tmp = SPM.Sess(ses).U(ind(1)).name{1};
            for c = 2:length(ind)
                 if contrasts(ui,ind(c)) == 1, sym = '+'; else sym = '-'; end
                 tmp = [tmp sym SPM.Sess(ses).U(ind(c)).name{1}];
            end
            DCM.U.name{ui} = tmp;  
            
            %DCM.U.idx(ui,:) = [ind zeros(Ncon-length(ind))]; % not sure this necessary
        end
        
        DCM.delays = repmat(SPM.xY.RT,Nvoi,1)/2;
        DCM.TE     = 0.03;  % 30ms on CBU
        
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 0;
        DCM.options.centre     = 1; % centre regressors (default = 1)
        DCM.options.nograph    = 1;
        DCM.options.maxnodes   = 8;
        DCM.options.maxit      = 128; %32;
        DCM.options.hidden     = [];
        DCM.options.induced    = 0;
        
        DCM.M.options = struct(); % needed else crashes in estimation

        for m = 1:length(Model)        
            DCM.a        = Model(m).a;
            DCM.b        = Model(m).b;
            DCM.b = DCM.b - diag(diag(DCM.b));
            DCM.c        = Model(m).c;
            DCM.d        = zeros(DCM.n,DCM.n,0); % needed else crashes in estimation
            
            filename = fullfile(out_path,sprintf('DCM_mod%d_ses%d.mat',m,ses));
             save(filename,'DCM');
            
            GCM{nsubses,m} = filename;    
        end
    end
end

cd(owd)
PEBwd = fullfile(owd,'PEB',dirname)
try mkdir(PEBwd); end
cd(PEBwd)
save('GCM_defined','GCM','-v7.3');




