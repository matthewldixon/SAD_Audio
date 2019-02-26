% Creates spm 1st level design matrix and estimates job 

%Script loads in basic experiment info and specifies condition
%onsets/durations and nuisance regressors: motion parameters and framewise displacement timecourse.  

%% Basic Info
clear all

spm('defaults','fmri');
spm_jobman('initcfg');

base_dir='/Path_to_data/'; %location of folder that contains each participant's data (e.g., S01, S02)
Regressor_dir='/Path_to_regressors/'; %location of matlab file containing regressor onsets, durations, and names
SessionList={'Run1' 'Run2' 'Run3'};
cd(base_dir)
dir_contents=dir; %list directory contents
sub_list = {dir_contents(3:end).name}; %list of subjects


fs = filesep;
nslices = 24; %number of slices acquired in each volume
TR = 1.5;
hpcutoff = 128; %high pass filter cutoff frequency

%%
for n = 1:length(sub_list)
    subj=sub_list{n}
    stats_dir=[base_dir sub_list{n} fs 'ses-pre' fs 'First_level_analysis_GLM1' fs]; % directory to save output files      
    cd(stats_dir); 
    if exist('SPM.mat', 'file'); delete('SPM.mat','*.nii','*.hdr');end  % Go to output directory and delete any existing files
    
  %% Load in condition onsets and durations from excel file
   
  for k = 1:length(SessionList) 
      cd(Regressor_dir);
  
    conditionFile=sprintf('conditions_run%d.mat',k); 
    conditionPath=[Regressor_dir conditionFile]; %Full pathway to conditions file
    
    %% Load in motion parameters file and FD timecourse  
    
    Functional_Dir=[base_dir sub_list{n} fs 'ses-pre' fs SessionList{k} fs]; %directory with subject's data for the current session
    cd(Functional_Dir);
    %motion_file_name=ls('rp*.txt'); % 6 motion parameters from realignment
    motion_file_name=spm_select('list',pwd,'^rp_sub*');
    %motion_file=textread([Functional_Dir fs motion_file_name]); %read in the textfile data
    motion_file=load(motion_file_name);
    
    FD_file_name=ls(fullfile(Functional_Dir,'Framewise_Displacement_timecourse.mat')); %file containing framewise displacement timecourse
    FD_file=load('Framewise_Displacement_timecourse');
    % outliers_file.Outlier_FD=[];
    R=[motion_file FD_file.fd]; %concatenate motion parameters and FD timecourse
    
    cd(Functional_Dir);
    %Create a matlab data file that will hold the nuisance (motion)
    %regressors info)
    mult_reg_file = sprintf('multireg_run%d.mat',k);
    save(mult_reg_file, 'R');
    multiregPath = [Functional_Dir fs mult_reg_file];  % Full pathway to motion regressors files
    
    %% Create Design Matrix
    cd(Functional_Dir);
    %f=spm_select('List',Functional_Dir,'swar*'); % Select smoothed normalised images
    func=ls('swar*.nii');
    func_file=[Functional_Dir func];
    func_expanded = cellstr(spm_select('expand', func_file)); %expand 4D into 3D data set
    func_data_files= cellstr(func_expanded(5:end,:)); %grab subset of volumes (ignore first 4 volumes)
      
    jobs{1}.stats{1}.fmri_spec.sess(k).scans=func_data_files; %add fMRI files for each run to the "jobs" structure 
    func = []; func_data_files = [];% clear temporary variables for next run
    
    jobs{1}.stats{1}.fmri_spec.sess(k).multi={conditionPath}; %add condition file to jobs structure
    jobs{1}.stats{1}.fmri_spec.sess(k).multi_reg={multiregPath}; %add motion regressor file to jobs structure
    jobs{1}.stats{1}.fmri_spec.sess(k).hpf = hpcutoff; % high pass filter
    jobs{1}.stats{1}.fmri_spec.sess(k).regress = struct([]);
    jobs{1}.stats{1}.fmri_spec.sess(k).cond = struct('name', {}, 'onset', {}, 'duration', {});
  end % session

%% Design parameters

jobs{1}.stats{1}.fmri_spec.dir={stats_dir}; % directory in which to save the output files
jobs{1}.stats{1}.fmri_spec.timing.units='scans'; %onsets and durations are specified in seconds
jobs{1}.stats{1}.fmri_spec.timing.RT=TR;
jobs{1}.stats{1}.fmri_spec.timing.fmri_t=nslices; %microtime resolution: the number of time bins that each TR is divided into when building regressors. Default is 16
jobs{1}.stats{1}.fmri_spec.timing.fmri_t0=nslices/2; %microtime onset: first  time-bin  at  which the regressors are resampled to coincide 
%with data acquisition.  If t0 = 1 then the regressors will be  appropriate  for  the  first  slice. If  you  want  to  temporally  
%realign  the regressors  so  that  they match responses in the middle slice then make t0 = microtime resolution /2
jobs{1}.stats{1}.fmri_spec.bases.hrf.derivs= [0 0]; % no basis functions
jobs{1}.stats{1}.fmri_spec.volt=1; %OPTIONS: 1 = do not model interactions 2 = model interactions. model interactions (Volterra) between trials or conditions 
%(e.g.,  because a BOLD response to a stimulus of a certain condition is expected to be larger when preceded by a trial of a certain other condition) 
jobs{1}.stats{1}.fmri_spec.global= 'None';% global normalisation.  Choose �scale� if you want to scale each voxel value of each scan to the global (overall) mean of that scan 
jobs{1}.stats{1}.fmri_spec.mask={fullfile(spm('Dir'),'toolbox', 'FieldMap','brainmask.nii')};% explicit masking.  analyses will be restricted to the space of the mask
jobs{1}.stats{1}.fmri_spec.cvi= 'AR(1)';% serial correlations
jobs{1}.stats{1}.fmri_spec.fact = struct('name', {}, 'levels', {});% no factorial design

%% run model specification

cd(stats_dir);
save specify.mat jobs % save job
disp(['RUNNING model specification for subject ' sub_list{n}]);
spm_jobman('run','specify.mat'); % run the job
clear jobs

% Ensure implicit masking is switched off so that no regions are
% automatically masked out due to low signal
load SPM
SPM.xM.TH = repmat(-Inf, length(SPM.xM.TH), 1);
SPM.xM.I = 0;
for k = 1:length(SessionList)
    for u = 1:length(SPM.Sess(k).U)
        SPM.Sess(k).U(u).orth = 0;
    end
end
save SPM SPM

%% Estimate
% setup job structure for model estimation and estimate

jobs{1}.stats{1}.fmri_est.spmmat = {[stats_dir 'SPM.mat']}; %Select SPM file 
save estimate.mat jobs %Save estimate job
disp(['RUNNING model estimation for subject ' sub_list{n}])
spm_jobman('run','estimate.mat'); %Run estimate job
clear jobs

%% Contrasts

jobs{1}.stats{1}.con.spmmat={[stats_dir fs 'SPM.mat']}; %specify SPM file
jobs{1}.stats{1}.con.delete=1; %delete existing contrasts
Contrast_Names={'reappraise_minus_react' 'accept_minus_react'};
%regressor order: 1) react, 2) reappraise, 3) mindful, 4) neg story, 5) neut story 6) neut beliefs, 7) ratings, 8) asterisk, 9) blank, 
contrast_vec={[-1 1 0 0 0 0 0 0 0] [-1 0 1 0 0 0 0 0 0]}; 
      
for i=1:length(contrast_vec)
         jobs{1}.stats{1}.con.consess{i}.tcon.name=Contrast_Names{i};
           ContrastVector=[contrast_vec{i} zeros(1,7) ...  %Run 1 regressor weights and zeros for 7 motion regressors (6 + framewise displacement)
           contrast_vec{i} zeros(1,7) ...  
           contrast_vec{i} zeros(1,7) ... 
           zeros(1,3)];                           %session constants 
        
        jobs{1}.stats{1}.con.consess{i}.tcon.convec=ContrastVector; %add contrast vector to jobs structure
        jobs{1}.stats{1}.con.consess{i}.tcon.sessrep='none';
end
    
    cd(stats_dir);
    save contrasts.mat jobs
    disp(['RUNNING contrast specification for subject number  ' sub_list{n}]);
    spm_jobman('run','contrasts.mat');
    disp(['Contrasts created for ' num2str(n) ' subjects']);
    clear jobs  


end
