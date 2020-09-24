clear all; close all; clc
format shortG

% Specify colors.
ventralcolor = [189 80 199]/255; % purple
vtpcolor = [161 95 84]; % brown
dorsalcolor = [224 83 114]/255; % salmon

blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Add path to plotting code.
addpath(genpath('/Volumes/240/lwx/bsc-track-plotting'));

% Should we save figures?
save_figures = 'yes';

% Specify desired view: sagittal, coronal, axial.
figView = 'sagittal';

if strcmp(figView, 'sagittal')
    slices = {[-1 0 0],[0 -1 0],[0 0 -1]};
elseif strcmp(figView, 'coronal')
    slices = {[-1 0 0],[0 -1 0],[0 0 -1]};
elseif strcmp(figView, 'axial')
    slices = {[-1 0 0],[0 -1 0],[0 0 -1]};
end

% Specify save directory.
saveDir = '/Volumes/240/lwx/mba_plots/';

% Specify number of streamlines to be plotted per tract.
subSample = 2000;

% Specify colors.
ifof = [142 198 255]/255; % light blue 
ilf = [0 127 255]/255; % dark blue

slf12 = [237 177 32]/255; % burnt yellow
slf3 = [240 221 165]/255; % light burnt yellow

parc = [64 224 208]/255; % turquoise
tpc =  [27 102 87]/255; % dark turquoise
mdlfspl = [42, 102, 0]/255; % green
mdlfang =  [207 255 226]/255; % sea foam

vof = [147 112 219]/255; % medium purple

fat = [240 128 128]/255; % light coral

% Set rootDir for where the data are.
rootDir = '/Volumes/Seagate/';

% Get contents of the directory where the tract measures for this subject are stored.
subfolders = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) ~= '.');

% Keep only names that are subject folders.
subfolders = subfolders(arrayfun(@(x) x.name(1), subfolders) == 's');

% Load in each tract's tractography measures for this subject.
for i = 25 %1:size(subfolders, 1)
    
    % Grab subID.
    sub(i) = str2num(subfolders(i).name(end-2:end));
    
    % Display current sub ID.
    disp(subfolders(i).name)
    
    % Get location of track.tck.
    temp = dir(fullfile(subfolders(i).folder, subfolders(i).name,  '/dt-neuro-track*/track.tck'));
    wbFG = fgRead(fullfile(temp.folder, temp.name)); clear temp;
    
    % Get location of classification.mat.
    temp = dir(fullfile(subfolders(i).folder, subfolders(i).name,  '/dt-neuro-wmc*/classification.mat'));
    classification = fullfile(temp.folder, temp.name); clear temp;
    
    % Get location of t1.
    temp = dir(fullfile(subfolders(i).folder, subfolders(i).name,  '/dt-neuro-anat*/t1.nii.gz'));
    t1 = niftiRead(fullfile(temp.folder, temp.name)); clear temp;
    
    % Set save directory and make, if necessary.
    saveDir = ['/Volumes/240/lwx/mba_plots/' num2str(sub(i))];
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
%     % Plot VH pathway.
%     subSelect = [10 30]; % left IFOF, left ILF
%     colors = num2cell(cat(1, ifof, ilf), 2);
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification,t1, figView, saveDir, subSelect, colors, subSample, slices)
%     
%     % Plot DH pathway.
%     subSelect = [14 16]; % left SLF12, left SLF3
%     colors = num2cell(cat(1, slf12, slf3), 2);
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
%     
%         
%     % pArc only.
%     subSelect = [36]; % left parc, left tpc, left mdlfspl, left mdlfang
%     colors = num2cell(cat(1, parc, tpc, mdlfspl, mdlfang), 2);
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
%     
%     
%     
%     % Plot VTP pathway.
%     subSelect = [36 38 34 32]; % left parc, left tpc, left mdlfspl, left mdlfang
%     colors = num2cell(cat(1, parc, tpc, mdlfspl, mdlfang), 2);
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
%     
%     % Plot VOF pathway.
%     subSelect = [58]; % left vof
%     colors = {vof};
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
%     
%     % Plot FAT pathway.
%     subSelect = [28]; % left fat
%     colors = {fat};
%     bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
% 
   % All pathways.
    subSelect = [58 10 30 36 38 34 32 14 16 28]; % left vof, left IFOF, left ILF, left parc, left tpc, left mdlfspl, left mdlfang, left SLF12, left SLF3, left FAT
    colors = num2cell(cat(1, vof, ifof, ilf, parc, tpc, mdlfspl, mdlfang, slf12, slf3, fat), 2);
    bsc_plotClassifiedStreamsAdaptive_v2(wbFG, classification, t1, figView, saveDir, subSelect, colors, subSample, slices)
    
end

% % Old color options. 
% ifof = [230 230 250]/255; % lavender
% ilf = [147 112 219]/255; % medium purple
% 
% slf12 = [139 35 35]/255; % brown4
% slf3 = [240 128 128]/255; % light coral
% 
% parc = [188 143 143]/255; % rosy brown
% tpc =  [139 69 19]/255; % saddle brown
% mdlfspl = [92 64 51]/255; % dark brown
% mdlfang =  [205 170 125]/255; % burlywood
% 
% vof = [0 127 255]/255; % slate blue
% 
% fat = [64 224 208]/255; % turquoise



