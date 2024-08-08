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
%The changes in this version include adjustments in participant selection and the absence of extracting the ROI of the EVC area.（2023）
%----------------------------------------------------------------

%% General
rwd     = 'F:\Face_DCM_fMRI-main/imaging/henson/users/sl01/FacesData/';
owd     = 'F:\Face_DCM_fMRI-main/imaging/henson/users/sl01/FacesData/';
dosubs  = [1 2 3 5 6 8 9 10 12 14  16 17 18 19 23 24]; % all except s11 for debrief
stat_dir = 'stat_concat_deb';
Nses = 1;
EOI_con = 2; % all 9 conditions vs baseline (ALWAYS CHECK!!!)
VOI_con = 14; % subject-specific OFA & FFA - should be faces - scr (ALWAYS CHECK!!!)

%% Select VOIs and Extract Time Series
ROInames = {'rOFA';'rFFA';'lOFA';'lFFA'};  
VOI_con = [ VOI_con VOI_con VOI_con VOI_con];
SS_VOIs = [ 1 1 1 1]; % EVC independent of subject (since no contrast can define them), OFA/FFA subject-specific
Masknames = {'rOFA_FWE';'rFFA_FWE';'lOFA_FWE';'lFFA_FWE'};
Nvoi = length(ROInames);

if any(SS_VOIs)
    %Grp_VOI_peak = nan(Nvoi,3);
    Grp_VOI_peak(1,:) = [39 -82 -10];
    Grp_VOI_peak(2,:) = [42 -46 -19];
    Grp_VOI_peak(3,:) = [-36 -85 -13];
    Grp_VOI_peak(4,:) = [-39 -49 -22]; 
end

for s = 1:length(dosubs)   
    sub = dosubs(s);
    cd(rwd)
    swd = sprintf('subject_%02d',sub);
    data_path = fullfile(rwd,swd,'bold',stat_dir);      
    out_path  = fullfile(owd,swd,stat_dir);   
    if ~exist(out_path); mkdir(out_path); end
    
    matlabbatch = {};
    matlabbatch{1}.cfg_basicio.cfg_cd.dir = cellstr(data_path);
    spm_jobman('run',matlabbatch);   
    
    % EXTRACTING TIME SERIES (one per VOI PER session):
    
    for v=1:Nvoi
        for ses=1:Nses
            clear matlabbatch
            matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(data_path,'SPM.mat'));
            matlabbatch{1}.spm.util.voi.adjust = EOI_con; % This is the contrast that specifies which conditions you will model in DCM (and so remove from the data effects of all other regressors in X)
            matlabbatch{1}.spm.util.voi.session = ses;
            matlabbatch{1}.spm.util.voi.name = fullfile(out_path,ROInames{v});            
            
            matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''};
            matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = VOI_con(v);  % Faces-Scr for OFA&FFA; 9 cons-baselinse for EVC
            %matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
            matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 1;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
            matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
            
            if SS_VOIs(v) == 1
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = Grp_VOI_peak(v,:);
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10;
                matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;
            else
                str = fullfile(owd,'GroupStats',stat_dir,[Masknames{v},'.nii,1']); % To use group-defined ROIs
                matlabbatch{1}.spm.util.voi.roi{2}.mask.image = {str};
                matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.5;
            end
            matlabbatch{1}.spm.util.voi.expression = 'i1&i2';

            spm_jobman('run',matlabbatch);
        end
    end     
end

VOInames = {};
Nscnvox = nan(length(dosubs),Nvoi,2);
for v = 1:Nvoi
    for ses = 1:Nses
        VOInames{v,ses} = ['VOI_' ROInames{v} sprintf('_%d',ses) sprintf('.mat')]; % VOInames{v,ses} = ['VOI_' ROInames{v} sprintf('_%d',ses)];
        for s = 1:length(dosubs)
            sub = dosubs(s);
            swd = sprintf('subject_%02d',sub);
            out_path = fullfile(owd,swd,stat_dir);
            VOI = load(fullfile(out_path,VOInames{v,ses}));
            Nscnvox(s,v,:) = size(VOI.xY.y);
        end
    end
end

%% Examine VOI univariate stats (THIS ASSUMES NSES=1!!!)

cols = [1:9]; %9con
betas = [];
for s = 1:length(dosubs)
    sub = dosubs(s);
    swd = sprintf('subject_%02d',sub);
    
    %% Get design matrix
    load(fullfile(rwd,swd,'bold',stat_dir,'SPM.mat'));
    for v = 1:Nvoi
        VOI = load(fullfile(owd,swd,stat_dir,VOInames{v}));
        tmp = pinv(SPM.xX.xKXs.X)*VOI.Y;
        betas(s,v,:) = tmp(cols);
    end
    
end

% IniF ImmF DelF    IniU ImmU DelU    IniS ImmS DelS
cw = [
     1/6 1/6 1/6         1/6 1/6 1/6         -1/3 -1/3 -1/3; % Perception - faces > scrambled
     1 -1 0              1 -1 0               1 -1 0; % Imm Rep
     1 0 -1              1 0 -1               1 0 -1; % Del Rep
     1/3 1/3 1/3         -1/3 -1/3 -1/3       0 0 0; % Recognition - famous > unfamous     
     1 -1/2 -1/2         1 -1/2 -1/2          0 0 0; % Imm vs Del Rep
     1/2 -1/2 0          1/2 -1/2 0           -1 1 0; % Perception X Imm Rep
     1/2 0 -1/2          1/2 0 -1/2           -1 0 1; % Perception X Del Rep
     1 -1 0              -1 1 0                0 0 0; % Recognition X Imm Rep
     1 0 -1              -1 0 1                0 0 0; % Recognition X Del Rep
    ];

T=[]; p=[]; figure(10); clf
sp = [2 4 6 1 3 5];
for v=1:Nvoi
    for c=1:size(cw,1)
         dc(:,c) = squeeze(betas(:,v,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    fprintf('%s\n',VOInames{v})
    disp([cw T' p'])
    
    figure(10),subplot(3,2,sp(v)),
    %boxplot(squeeze(dc))
    %boxplot(squeeze(betas(:,v,:)))
    bar(mean(squeeze(betas(:,v,:))))
    title(VOInames{v})
end