% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% for later plotting: http://people.duke.edu/~jmp33/matlab/plotting_intro.html

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% beh_measure = 'age'; %age, lit, vm, fm
measure = {'volume'};
roi = 'lobe'; %'tractendpoints';
streamline_min = 100;

% Read in lut for glasser atlas.
if strcmp(roi, 'lobe')
    lobe = readtable([rootDir 'supportFiles/LUT_glasser.csv']);
elseif strcmp(roi, 'tractendpoints')
    lobe = readtable([rootDir 'supportFiles/LUT_tractendpoints.csv']);
end

occipital = lobe.Name(lobe.Lobe == 1);
ventral = lobe.Name(lobe.Lobe == 2);
parietal = lobe.Name(lobe.Lobe == 3);
frontal = lobe.Name(lobe.Lobe == 4);

%% WHITE MATTER MEASURES
for w = 1:length(measure)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load(fullfile(rootDir, 'supportFiles', ['LWX_CM_data_' measure{w} '_singleshell.mat']))
    d = data_tbl;
    clear data_tbl
    
    % Convert into array and header for ease.
    %     data_all_in = table2array(data_tbl);
    %     data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Select tois from d.
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices tracts of interest: col.
        t_idx(k) = strcmp(d.Properties.VariableNames{k}, 'subID') ...
            || ismember(d.Properties.VariableNames{k}, occipital) || ismember(d.Properties.VariableNames{k}, ventral) ...
            || ismember(d.Properties.VariableNames{k}, parietal) ... 
            || ismember(d.Properties.VariableNames{k}, frontal) || strcmp(d.Properties.VariableNames{k}, 'icv');
        
    end
    roi = d(:, t_idx);
    
    % Return to table format.
    %     roi = array2table(t, 'VariableNames', d.Properties.VariableNames(t_idx));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(table2array(roi(find(d.gp_age ~= 3), :)));
    m_wtwg2 = nanmean(table2array(roi(find(d.gp_age == 3), :)));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(table2array(roi(find(d.gp_age ~= 3), :)));
    std_wtwg2 = nanstd(table2array(roi(find(d.gp_age == 3), :)));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2;
    r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, size(find(d.gp_age~=3))), repmat(r_max_wtwg2, size(find(d.gp_age==3))));
    r_min = cat(1, repmat(r_min_wtwg1, size(find(d.gp_age~=3))), repmat(r_min_wtwg2, size(find(d.gp_age==3))));
    
    % Display.
    disp([measure{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(table2array(roi) > r_max))
        disp(['Replaced ' num2str(numel(find(table2array(roi) > r_max))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.'])
        [idx, idy] = find(table2array(roi) > r_max);
        if ~isempty(idx)
            roi(idx, idy) = {NaN};
        end
    else
        disp('No data points were above 3 standard deviations of the within-tract, within-group mean.')
    end
    %Min
    if ~isempty(find(table2array(roi) < r_min))
        disp(['Replaced ' num2str(numel(find(table2array(roi) < r_min))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.'])
        [idx, idy] = find(table2array(roi) < r_min);
        if ~isempty(idx)
            roi(idx, idy) = {NaN};
        end
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.')
    end
    
    % Create lobe averages.
    for k = 1:length(roi.Properties.VariableNames)
        
        o_idx(k) = ismember(roi.Properties.VariableNames{k}, occipital);
        v_idx(k) = ismember(roi.Properties.VariableNames{k}, ventral);
        p_idx(k) = ismember(roi.Properties.VariableNames{k}, parietal);
        f_idx(k) = ismember(roi.Properties.VariableNames{k}, frontal);
        
    end
    
    % Get mean in each cortical region.
    odata = nanmean(table2array(roi(:, o_idx)), 2);
    vdata = nanmean(table2array(roi(:, v_idx)), 2);
    pdata = nanmean(table2array(roi(:, p_idx)), 2);
    fdata = nanmean(table2array(roi(:, f_idx)), 2);
    
%     % Creat lobe snr averages.
%     for k = 1:length(d.Properties.VariableNames)
%         
%         o_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(occipital, '_snr'));
%         v_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(ventral, '_snr'));
%         p_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(parietal, '_snr'));
%         f_snr_idx(k) = ismember(d.Properties.VariableNames{k}, strcat(frontal, '_snr'));
%         
%     end
%     
%     % Get mean in each cortical region.
%     osnrdata = nanmean(table2array(d(:, o_snr_idx)), 2);
%     vsnrdata = nanmean(table2array(d(:, v_snr_idx)), 2);
%     psnrdata = nanmean(table2array(d(:, p_snr_idx)), 2);
%     fsnrdata = nanmean(table2array(d(:, f_snr_idx)), 2);
%     
%             osnr = nanmean(table2array(d(:, o_snr_idx)), 2);
%         vsnr = nanmean(table2array(d(:, v_snr_idx)), 2);
%         psnr = nanmean(table2array(d(:, p_snr_idx)), 2);
%         fsnr = nanmean(table2array(d(:, f_snr_idx)), 2);
       
       
        % Output csv file for ANOVA in SPSS. (Matlab doesn't handle Mixed Model
        % ANOVAs well when the between-group variable is correlated with subID
        % (e.g., when between-group variable is something like age groups).
        t_out = array2table(cat(2, d.subID, d.gp_age, d.gp_lit, d.gp_vm, ...
            d.gp_fm, d.cov_age, d.cov_sex, d.c_lit, d.c_vm, d.c_fm, d.BeeryVMI, d.BeeryVP, d.BeeryMC, d.gPegs_dom_forward, d.gPegs_nondom_forward, ...
            d.WJIV_LetterWordIdentification, d.WJIV_Spelling, d.WJIV_WordAttack, d.WJIV_SpellingOfSounds, table2array(roi(:, 2:end)), odata, vdata, pdata, fdata), ...
            'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', 'c_lit', 'c_vm', 'c_fm', 'BeeryVMI', 'BeeryVP', 'BeeryMC', 'gPegs_dom_forward', 'gPegs_nondom_forward', ...
            'WJIV_LetterWordIdentification', 'WJIV_Spelling', 'WJIV_WordAttack', 'WJIV_SpellingOfSounds', roi.Properties.VariableNames{2:end}, 'occipital', 'ventral', 'parietal', 'frontal'});

    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['LWX_CM_data_forSPSS_' measure{w} '_singleshell.csv']));
    fid = fopen(fullfile(rootDir, 'supportFiles', ['LWX_CM_data_forSPSS_' measure{w} '_singleshell.csv']));
    fclose(fid)
    
end


