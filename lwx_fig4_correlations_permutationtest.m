clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

multcompcorrection = 'no';

group = 'children'; % adults, children

dorsalcolor = [236 176 32]/255; %burnt yellow
ventralcolor= [14 114 184]/255; % blue

%% READ IN DATA AND ORGANIZE.

% Children.
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));
d = d(d.group_age3 ~= 3, :);

tpc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC_z')))));
pArc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc_z')))));
mdlfspl = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl_z')))));
mdlfang = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang_z')))));
vof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF_z')))));
aslant = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant_z')))));
slf12 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2_z')))));
slf3 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3_z')))));
ilf = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightILF_z')))));
ifof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF_z')))));

% Compute correlations among tracts.
[rho_child, pmat_child] = corr(cat(2, aslant, slf12, slf3, mdlfang, mdlfspl, tpc, pArc, ilf, ifof, vof), 'rows', 'pairwise');
    
clear d tpc pArc mdlfspl mdlfang vof aslant slf12 slf3 ilf ifof

% Adults
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));
d = d(d.group_age3 == 3, :);
    
% SELECT the measurements of the tracts that I care about.
tpc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC'))))); tpc = (tpc-nanmean(tpc))/nanstd(tpc);
pArc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc'))))); pArc = (pArc-nanmean(pArc))/nanstd(pArc);
mdlfspl = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl'))))); mdlfspl = (mdlfspl-nanmean(mdlfspl))/nanstd(mdlfspl);
mdlfang = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang'))))); mdlfang = (mdlfang-nanmean(mdlfang))/nanstd(mdlfang);
vof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF'))))); vof = (vof-nanmean(vof))/nanstd(vof);
aslant = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant'))))); aslant = (aslant-nanmean(aslant))/nanstd(aslant);
slf12 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2'))))); slf12 = (slf12-nanmean(slf12))/nanstd(slf12);
slf3 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3'))))); slf3 = (slf3-nanmean(slf3))/nanstd(slf3);
ilf = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightILF'))))); ilf = (ilf-nanmean(ilf))/nanstd(ilf);
ifof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF'))))); ifof = (ifof-nanmean(ifof))/nanstd(ifof);

% Compute correlations among tracts.
[rho_adult, pmat_adult] = corr(cat(2, aslant, slf12, slf3, mdlfang, mdlfspl, tpc, pArc, ilf, ifof, vof), 'rows', 'pairwise');

% z_aslant, z_slf12, z_slf3, z_mdlfang, z_mdlfspl, z_tpc, z_pArc, z_ilf, z_ifof, z_vof
G = [0 1 1 3 3 3 3 2 2 0];

% Get real distribution of vv-vd difference: Ventral-vertical correlation greater than dorsal-vertical correlation.
vv_adult = rho_adult(find(G==2), find(G==3)); vv_adult = vv_adult(:);
vd_adult = rho_adult(find(G==1), find(G==3)); vd_adult = vd_adult(:);
vv_child = rho_child(find(G==2), find(G==3)); vv_child = vv_child(:);
vd_child = rho_child(find(G==1), find(G==3)); vd_child = vd_child(:);

%% Permutation testing for interaction between correlations pair (dorsal-vertical and dorsal-ventral correlations) and age group (children and adults).

% Get real distribution.
interaction = vv_child - vd_child - (vv_adult - vd_adult); 
mu_interaction = nanmean(interaction(:));
sigma_interaction = nanstd(interaction(:));

% Create null distribution of vv-vd difference, where the correlation category (i.e., vertical-ventral, vertical-dorsal) and 
% age category (i.e., children, adult) labelling are broken.
null_dis = [vv_child(:); vd_child(:); vv_adult(:); vd_adult(:)];
for i = 1:10000
    
    % Randomly select a permutation with replacement.
    this_vv_child = randsample(null_dis, size(vv_child, 1), true);
    
    % Randomly select a permutation with replacement.
    this_vd_child = randsample(null_dis, size(vd_child, 1), true);
    
    % Randomly select a permutation with replacement.
    this_vv_adult = randsample(null_dis, size(vv_adult, 1), true);
    
    % Randomly select a permutation with replacement.
    this_vd_adult = randsample(null_dis, size(vd_adult, 1), true);
    
    % Get the correlation at that location.
    interaction_null(i) = mean(this_vv_child) - mean(this_vd_child) - (mean(this_vv_adult) - mean(this_vd_adult));
    
    clear this_vv this_vd
    
end
mu_interaction_null=nanmean(interaction_null);
sigma_interaction_null=nanstd(interaction_null);

% Test for significance.
z = (mu_interaction - mu_interaction_null)./sigma_interaction_null;
p = 1-normcdf(abs(z), 0, 1);

disp(['z = ' num2str(z) ', p = ' num2str(p)]);

disp(['mean vv child = ' num2str(mean(vv_child, 'all')) ', mean vd child = ' num2str(mean(vd_child, 'all'))]);
disp(['mean vv adult = ' num2str(mean(vv_adult, 'all')) ', mean vd adult = ' num2str(mean(vd_adult, 'all'))]);

