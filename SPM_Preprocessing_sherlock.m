function SPM_Preprocessing_sherlock(subj)
%
% This script initiates SPM to conduct several preprocessing steps for
% functional MRI data.
%A structure called "jobs" is created and saved as a .mat file; it holds all of the parameters for
%each preprocessing step and then is executed.
%
%Sequence:
%1)realigment
%2)slice-time correction (optional)
%3)coregistration
%4)segmentation
%5)Normalization
%6)Smoothing
%

%% Experiment parameters
Defaults = spm_get_defaults;
%spm('Defaults','fMRI'); % Initialise SPM defaults
spm_jobman('initcfg'); %
fs=filesep;
TR= 1.5;
NumSlices=24;
SliceAcquisition = [1:2:NumSlices 2:2:NumSlices]; %This specifies that slices were collected in an interleaved fashion starting from the bottom.
Smoothing_Kernal = [8 8 8]; %(FWHM)
numRuns=3; % number of functional runs
data_dir='/path_to_data/';
Session='ses-pre';

jobs_dir=[data_dir subj fs Session fs 'Jobs/'];



%% Unzip and expand files 
%
for k=1:numRuns
func_dir=[data_dir subj fs Session fs 'Run' int2str(k) fs];
cd(func_dir)
%func_file=gunzip([func_dir '*.nii.gz']); %unzip functional data if necessary
func=ls('*.nii');
func_file=[func_dir func]

func_expanded = cellstr(spm_select('expand', func_file)); %expand 4D into 3D data set
func_data_files{:,k} = cellstr(func_expanded(5:end,:)); %grab subset of volumes (ignore first 4 volumes)
end

anat_dir=[data_dir subj fs Session fs 'anat' fs];
cd(anat_dir)
%anat_file=gunzip([anat_dir '*.nii.gz']); %unzip anatomical data if necessary
anat=ls('*.nii'); 
anat_file=[anat_dir anat];

anat_expanded = cellstr(spm_select('expand', anat_file));


cd(jobs_dir); %Go to jobs dir, where we will save all the job files

%% REALIGNMENT - estimate & reslice (i.e, motion correction)
for k=1:numRuns
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1,k}=func_data_files{k}; %cellstr(func_files);
end
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp=4; % set interpolation to 4th degree
matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality=1;      % default 0.9
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1]; % Create a mean image, and reslice other files

%clear func_files func
save('realign.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch

%% Compute framewise displacement
%adapted from from Chris Rorden's nii_qa_moco script, which is based on Powers et al. http://www.ncbi.nlm.nih.gov/pubmed/22019881 

for k=1:numRuns
    cd([data_dir subj fs Session fs 'Run' int2str(k) fs]);
    spm_select('List', fullfile(data_dir,subj, '/', Session,'/Run1/'), 'rp*');
    motion_params_file=spm_select('List', strcat(data_dir,subj, '/', Session,'/Run',int2str(k),'/'), '/*.txt'); %select the rp*.txt file
    rp=load(motion_params_file);
    rpdx = rp;%preallocate, set
    for i = 1:6
        rpdx(2:end,i) = rp(2:end,i) - rp(1:end-1,i); %compute temporal diffs
    end
    %fprintf('Assuming motion parameters saved in format "Xmm Ymm Zmm Xrad Yrad Zrad" \n');
    for i = 4:6 %convert rotations from radians to mm
        rpdx(1:end,i) = rpdx(1:end,i) * 50; %brain assumed 50mm radius: circle circumference=2*pi*r http://en.wikipedia.org/wiki/Great-circle_distance
    end

rpdx = abs(rpdx);
fd = sum(rpdx,2);
save('Framewise_Displacement_timecourse.mat', 'fd'); %Save framewise displacement timecourse, to be used as nuissance regressor in first level GLM
end

cd(jobs_dir);
    

%% Slice Time correction
for k=1:numRuns
    r_func_file=spm_select('List', strcat(data_dir,subj, '/', Session,'/Run', int2str(k),'/'), '^rsub*'); %grab realigned data
    r_func_expanded = cellstr(spm_select('expand',strcat(data_dir,subj, '/', Session,'/Run',int2str(k),'/', r_func_file)));
    matlabbatch{1}.spm.temporal.st.scans{1, k} = r_func_expanded; %cellstr(func_files);  %put the functional files (f) in a cell array within the jobs structure that SPM works with.
end

matlabbatch{1}.spm.temporal.st.nslices=NumSlices;
matlabbatch{1}.spm.temporal.st.tr=TR;
matlabbatch{1}.spm.temporal.st.ta=TR-TR/NumSlices;
matlabbatch{1}.spm.temporal.st.so=SliceAcquisition;
matlabbatch{1}.spm.temporal.st.refslice=SliceAcquisition(floor(length(SliceAcquisition)/2));
matlabbatch{1}.spm.temporal.st.prefix='a';

save('slicetiming.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch

%% COREGISTRATION
%(align functional and structural data in space)
%Align the structural to the mean functional using normalized mutual
%information (NMI)
mean_func=spm_select('List', fullfile(data_dir,subj, '/', Session,'/Run1/'), 'mean*'); 
matlabbatch{1}.spm.spatial.coreg.estimate.ref = cellstr([data_dir subj fs Session fs 'Run1' fs mean_func]);
matlabbatch{1}.spm.spatial.coreg.estimate.source = cellstr(anat_file);

save('coreg.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch


%% NORMALIZATION - estimate and write (functional) 
conCat_files=[];
for k=1:numRuns
a_func_file=spm_select('List', strcat(data_dir,subj, '/', Session,'/Run', int2str(k),'/'), '^a*'); %grab slicetime corrected data
a_func_expanded = cellstr(spm_select('expand',strcat(data_dir,subj, '/', Session,'/Run',int2str(k),'/', a_func_file)));
conCat_files=[conCat_files;a_func_expanded]; %concatenate files from all runs
clear a_func_expanded
end


matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = anat_expanded; % image to align - use anatomical to determine parameters
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = conCat_files; %data to be normalized

matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [3 3 3]; %voxel sizes
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';


save('norm_func.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);

clear func func_files conCat_files
clear matlabbatch
        
%% NORMALIZATION - estimate and write (structural) 
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = anat_expanded; % use anatomical to determine parameters
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = anat_expanded; %data to be normalized

matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1]; %voxel sizes
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

save('norm_struct.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);
clear matlabbatch
%}
        
%% SMOOTHING
%spatial smoothing of the signal to minimize the effect of individual
%differences in anatomy and activation patterns.
conCat_files_normalized=[];

for k=1:numRuns
w_func_file=spm_select('List', strcat(data_dir,subj, '/', Session,'/Run', int2str(k),'/'), '^w*'); %grab normalized data
w_func_expanded = cellstr(spm_select('expand',strcat(data_dir,subj, '/', Session,'/Run',int2str(k),'/', w_func_file)));
conCat_files_normalized=[conCat_files_normalized;w_func_expanded]; %concatenate files from all runs
clear w_func_expanded
end

matlabbatch{1}.spm.spatial.smooth.data = conCat_files_normalized;
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

save('smooth.mat', 'matlabbatch');
spm_jobman('run', matlabbatch);

clear W_func_file W_func_expanded conCat_files_normalized
clear matlabbatch
