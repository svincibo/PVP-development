function [] = main()

% set file paths
t1 = fullfile('path/to/t1.nii.gz');
wbfg = fullfile('path/to/track.tck');
wmc = fullfile('path/to/classification.mat');
measureFile = fullfile('path/to/fa.nii.gz'); % or any other diffusion measure nifti

% load data
anat = niftiRead(t1);
measure= niftiRead(measureFile);
wbFG = fgRead(wbfg);
load(wmc);

% extract track of interest as a fg
tractIndex = 40 % left meyers loop
fg_classified = bsc_makeFGsFromClassification_v5(classification,wbFG);
fg = fg_classified{tractIndex};

% set subset
subset = 100;

% set number of nodes
numnodes = 200;

% set slices
% if axial; if coronal; slices = [0 -1 0]; if sagittal; slices = [-1 0 0] 
slices = [0 0 -1];

% set figview
figView = 'axial';

% set radius
radius = 2;

% set save name
saveName = 'test';

tractTubes(anat,fg,numnodes,slices,subset,measure,radius,figView,saveName)
