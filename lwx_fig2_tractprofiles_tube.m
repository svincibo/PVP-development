% set file paths
rootDir = '/Volumes/240/lwx/proj-5e849e65952fef3dcd7a1700/sub-128/';

t1 = fullfile(rootDir, 'dt-neuro-anat-t1w.tag-acpc_aligned.tag-preprocessed.id-5f075035880bc1cc69a708c4', 't1.nii.gz');
wbfg = fullfile(rootDir, 'dt-neuro-track-tck.id-5f07ac47f7e21ab7856e1f75/track.tck');
wmc = fullfile(rootDir, 'dt-neuro-wmc.tag-cleaned.id-5f07eda8f7e21a9b1a6f1afa/classification.mat');
measureFile = fullfile(rootDir, 'dt-neuro-tensor.id-5f07ac2ef7e21a766a6e1ebb/fa.nii.gz'); % or any other diffusion measure nifti

% Add path to plotting code.
addpath(genpath('/Volumes/240/lwx/bsc-track-plotting'));

% Add path to plotting code.
addpath(genpath('/Volumes/240/lwx/plot-tract-tubes'));

% Add path to spm.
addpath(genpath('~/Documents/spm12'));

% load data
anat = niftiRead(t1);
measure= niftiRead(measureFile);
wbFG = fgRead(wbfg);
load(wmc);

% extract track of interest as a fg
tractIndex = 38;
fg_classified = bsc_makeFGsFromClassification_v5(classification,wbFG);
fg = fg_classified{tractIndex};

% set color of tract
if tractIndex == 34 %left mdlfspl
    
    % set save name
    saveName = fullfile('mdlfspl');
    
    % set color of streamlines
    c = [42, 102, 0]/255; % green
    
    % set subset
    subset = 100; % suggestion = 100
    
    % set radius
    radius = 4;
    
elseif tractIndex == 32 %left mdlfang
    
    % set save name
    saveName = fullfile('mdlfang');
    
    % set color of streamlines
    c = [207 255 226]/255; % sea foam
    
    % set subset
    subset = 100; % suggestion = 100
    
    % set radius
    radius = 4;
    
    elseif tractIndex == 36 %left parc
    
    % set save name
    saveName = fullfile('parc');
    
    % set color of streamlines
    c = [64 224 208]/255; % turquoise
    
    % set subset
    subset = 100; % suggestion = 100
    
    % set radius
    radius = 4;
    
    elseif tractIndex == 38 %left tpc
    
    % set save name
    saveName = fullfile('tpc');
    
    % set color of streamlines
    c = [27 102 87]/255; % dark turquoise
    
    % set subset
    subset = 100; % suggestion = 100
    
    % set radius
    radius = 4;
    
% elseif tractIndex == 10 %lifof
%     
%     % set save name
%     saveName = fullfile('lifof');
%     
%     % set color of streamlines
%     c = [142 198 255]/255; % light blue
%     
%     % set subset
%     subset = 100; % suggestion = 100
%     
%     % set radius
%     radius = 4;
%     
% elseif tractIndex == 30 %lilf
%     
%     % set save name
%     saveName = fullfile('lilf');
%     
%     % set color of streamlines
%     c = [0 127 255]/255; % dark blue
%     
%     % set subset
%     subset = 100; % suggestion = 100
%     
%     % set radius
%     radius = 4;
%     
% elseif tractIndex == 14 %lslf1and2
%     
%     % set save name
%     saveName = fullfile('slf1and2');
%     
%     % set color of streamlines
%     c = [237 177 32]/255; % burnt yellow
%     
%     % set subset
%     subset = 1000; % suggestion = 100
%     
%     % set radius
%     radius = 15;
%     
% elseif tractIndex == 16 %lslf3
%     
%     % set save name
%     saveName = fullfile('slf3');
%     
%     % set color of streamlines
%     c = [240 221 165]/255; % light burnt yellow
%     
%     % set subset
%     subset = 1000; % suggestion = 100
%     
%     % set radius
%     radius = 9;
%     
end

% set number of nodes
numnodes = 200;

% Specify desired view: sagittal, coronal, axial.
figView = 'sagittal';

if strcmp(figView, 'sagittal')
    slices = [-1 0 0];
elseif strcmp(figView, 'coronal')
    slices = [0 -1 0];
elseif strcmp(figView, 'axial')
    slices = [0 0 -1];
end

% specify savedir
saveDir = '/Volumes/240/lwx/tractprofiletubes/';

tractTubes_v2(anat,fg,numnodes,slices,subset,c,measure,radius,figView,saveDir,saveName)
