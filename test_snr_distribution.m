clear all; close all; clc
format shortG
ls

details = 'sub-103-multishell-posttopuppostdenoise';

% Set working directories.
rootDir = '/Volumes/240/lwx';

% Get the directory where the dwi image is stored.
dwi_data_path = dir(fullfile(rootDir, 'testing', details, 'dwi.nii.gz'));

% Load nii file.
dwi = niftiRead(fullfile(dwi_data_path.folder, dwi_data_path.name));

% Read in the bval files to identify b = 0 volumes.
bval = dlmread(fullfile(rootDir, 'testing', details, 'dwi.bvals'));

% Get the directory where the mask image is stored.
mask_data_path = dir(fullfile(rootDir, 'testing', details, 'mask_mask.nii.gz'));

% Load nii file, dimensions 1, 2, and 3.
mask = niftiRead(fullfile(mask_data_path.folder, mask_data_path.name));

% Find the non b = 0 volumes, dimension 4.
dwi_idx = find(bval ~= 0);

% Calculate the standard deviation of the background noise.
std_noise = std(double(dwi.data(1:10, 1:10, 1:10, dwi_idx)), [], 'all');

% Calculate the mean of the signal inside at each voxel inside the brain across non-b0 volumes.
signal = double(dwi.data) .* double(mask.data);
mean_signal = mean(signal(:, :, :, dwi_idx), 4);

% Compute SNR in only the brain voxels and only the non-b0 volumes.
snr = signal./std_noise; % mean(y(inside mask))/std(y(outside mask))

% Assign to dwi.data for easy output.
dwi.data = snr;

% Write out new nifti file with values that are the probability that that voxel will have ICVF = 1.
niftiWrite(dwi, fullfile(rootDir, 'testing', details, [details '_b_snr.nii.gz']))







