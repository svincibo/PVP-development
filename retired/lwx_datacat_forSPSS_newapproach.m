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
wm_measure_here = {'fa'};
streamline_min = 100;

% Load wm data from lwx_qa_tractstats.m.
load([rootDir 'supportFiles/LWX_data_streamlinecount.mat']);
streamline_idx = table2array(streamlinecounts) >= streamline_min;

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure_here{w} '_singleshell.mat']))
    d = data_tbl;
    clear data_tbl
    
    % Convert into array and header for ease.
    %     data_all_in = table2array(data_tbl);
    %     data_all_in_header = data_tbl.Properties.VariableNames;
    
    % Select tois from d.
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices tracts of interest: col.
        t_idx(k) = strcmp(d.Properties.VariableNames{k}, 'subID') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2') || strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftSLF3') || strcmp(d.Properties.VariableNames{k}, 'rightSLF3') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftAslant') || strcmp(d.Properties.VariableNames{k}, 'rightAslant') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftILF') || strcmp(d.Properties.VariableNames{k}, 'rightILF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftIFOF') || strcmp(d.Properties.VariableNames{k}, 'rightIFOF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftTPC') || strcmp(d.Properties.VariableNames{k}, 'rightTPC') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftpArc') || strcmp(d.Properties.VariableNames{k}, 'rightpArc') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftVOF') || strcmp(d.Properties.VariableNames{k}, 'rightVOF') ...
            || strcmp(d.Properties.VariableNames{k}, 'leftMDLFang') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFang');
        
    end
    t_temp = d(:, t_idx);
    
    % Change cells with low streamline counts to NaN.
    for k = 1:length(t_temp.Properties.VariableNames)
        
        % Reorganize streamline_idx into lostream_idx so that the columns are in the same order as toi_idx.
        col_idx = find(strcmp(streamlinecounts.Properties.VariableNames, t_temp.Properties.VariableNames{k}));
        
        if ~isempty(col_idx)
            
            t(:, k) = streamline_idx(:, col_idx).*table2array(t_temp(:, k));
            
        end
        
    end
    t(t == 0) = NaN;
    
    % Return to table format.
    toi = array2table(t, 'VariableNames', d.Properties.VariableNames(t_idx));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(table2array(toi(find(d.gp_age == 1), :)));
    m_wtwg2 = nanmean(table2array(toi(find(d.gp_age == 2), :)));
    m_wtwg3 = nanmean(table2array(toi(find(d.gp_age == 3), :)));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(table2array(toi(find(d.gp_age == 1), :)));
    std_wtwg2 = nanstd(table2array(toi(find(d.gp_age == 2), :)));
    std_wtwg3 = nanstd(table2array(toi(find(d.gp_age == 3), :)));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2;
    r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3;
    r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, size(find(d.gp_age==1))), repmat(r_max_wtwg2, size(find(d.gp_age==2))), ...
        repmat(r_max_wtwg3, size(find(d.gp_age==3))));
    r_min = cat(1, repmat(r_min_wtwg1, size(find(d.gp_age==1))), repmat(r_min_wtwg2, size(find(d.gp_age==2))), ...
        repmat(r_min_wtwg3, size(find(d.gp_age==3))));
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(table2array(toi) > r_max))
        disp(['Replaced ' num2str(numel(find(table2array(toi) > r_max))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.'])
        idx = find(table2array(toi) > r_max);
        if ~isempty(idx)
            toi(idx) = NaN;
        end
    else
        disp('No data points were above 3 standard deviations of the within-tract, within-group mean.')
    end
    %Min
    if ~isempty(find(table2array(toi) < r_min))
        disp(['Replaced ' num2str(numel(find(table2array(toi) < r_min))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.'])
        idx = find(table2array(toi) > r_min);
        if ~isempty(idx)
            toi(idx) = NaN;
        end
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.')
    end
    
    %     % Make hemisphere averaged column for each tract.
    %     toi_havg = cat(2, nanmean(toi(:, [1 11]), 2), nanmean(toi(:, [2 12]), 2), nanmean(toi(:, [3 13]), 2), nanmean(toi(:, [4 14]), 2), nanmean(toi(:, [5 15]), 2), ...
    %         nanmean(toi(:, [6 16]), 2), nanmean(toi(:, [7 17]), 2), nanmean(toi(:, [8 18]), 2), nanmean(toi(:, [9 19]), 2), nanmean(toi(:, [10 20]), 2));
    %     toi_havg_header = erase(data_all_in_header(1:10), 'left');
    
    % Make tract group
    
    
    % Output csv file for ANOVA in SPSS. (Matlab doesn't handle Mixed Model
    % ANOVAs well when the between-group variable is correlated with subID
    % (e.g., when between-group variable is something like age groups).
    t_out = array2table(cat(2, d.subID, d.gp_age, d.gp_lit, d.gp_vm, ...
        d.gp_fm, d.cov_age, d.cov_sex, table2array(toi(:, 2:end))), ...
        'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', toi.Properties.VariableNames{2:end}});
    
    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell.csv']));
    fid = fopen(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell.csv']));
    fclose(fid)
    
    %     % Output z-scored file for SPSS.
    toi_z = (nanmean(table2array(toi), 1) - table2array(toi))./nanstd(table2array(toi), [], 1);
%     toi_havg_z = (nanmean(toi_havg, 1) - toi_havg)./nanstd(toi_havg, [], 1);
    
    %     temp = nanmean(toi(:, hv == 1), 2);
    %     toi_meanhd_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    %
    %     temp = nanmean(toi(:, hv == 2), 2);
    %     toi_meanhv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    %
    %     temp = nanmean(toi(:, hv == 3), 2);
    %     toi_meanv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    t_out_z = array2table(cat(2, d.subID, d.gp_age, d.gp_lit, d.gp_vm, ...
        d.gp_fm, d.cov_age, d.cov_sex, toi_z(:, 2:end)), ...
        'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', toi.Properties.VariableNames{2:end}});
    writetable(t_out_z, fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell_z.csv']));
    fid = fopen(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell_z.csv']));
    fclose(fid)
    
end


