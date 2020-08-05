clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% beh_measure = 'age'; %age, lit, vm, fm
wm_measure_here = {'fa'}; 

%% WHITE MATTER MEASURES
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    load(fullfile(rootDir, 'supportFiles', ['LWX_data_' wm_measure_here{w} '_singleshell.mat']))
    
    wm = data_tbl;

    % Should outliers be removed? If so, which subIDs?
    remove_outliers = 'yes';
    if strcmp(remove_outliers, 'yes')
        
        % Identify outliers to be removed - conservative removal.
        %         outlier = [108 126 212 214 318];
        % 108, snr is below 2 SD of group mean
        % 126, dwi image has major distortions, visual inspection
        % 212, physical anomaly that precludes tracking of vertical tracks, visual inspection
        % 214, major motion artifacts, visual inspection
        % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
        
        % Identify outliers to be removed - liberal removal.
        outlier = [108 116 119 125 126 206 212 214 303 317 318];
        % 116, FD > 2
        % 119, FD > 2
        % 125, FD > 2
        % 206, FD > 2
        % 303, SNR < 15
        % 317, FD > 2
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(wm.subID, outlier);
        
        % Remove outliers.
        wm = wm(~idx_outlier, :);
        
    end
    
    % Convert into array and header for ease.
    wm_header = wm.Properties.VariableNames;
    wm_data = table2array(wm);
    
    % Get index matrices for hypothesis-driven grouping of WM tracts. 
    for k = 1:length(wm_header)
        
        % Indices of horizontal tracts in dorsal stream.
        dorsal_idx(k) = strcmp(wm_header{k}, 'leftSLF1And2') || strcmp(wm_header{k}, 'rightSLF1And2') ...
            || strcmp(wm_header{k}, 'leftSLF3') || strcmp(wm_header{k}, 'rightSLF3');
        
        % Indices of horizontal tracts in ventral stream.
        ventral_idx(k) = strcmp(wm_header{k}, 'leftILF') || strcmp(wm_header{k}, 'rightILF') ...
            || strcmp(wm_header{k}, 'leftIFOF') || strcmp(wm_header{k}, 'rightIFOF');
        
        
        % Indices of posterior vertical tracts.
        vertical_idx(k) = strcmp(wm_header{k}, 'leftTPC') || strcmp(wm_header{k}, 'rightTPC') ...
            || strcmp(wm_header{k}, 'leftpArc') || strcmp(wm_header{k}, 'rightpArc') ...
            || strcmp(wm_header{k}, 'leftMDLFspl') || strcmp(wm_header{k}, 'rightMDLFspl') ...
            || strcmp(wm_header{k}, 'leftMDLFang') || strcmp(wm_header{k}, 'rightMDLFang');
        
        % Indices of vof.
        vof_idx(k) = strcmp(wm_header{k}, 'leftVOF') || strcmp(wm_header{k}, 'rightVOF');
        
        % Indices of aslant
        aslant_idx(k) = strcmp(wm_header{k}, 'leftAslant') || strcmp(wm_header{k}, 'rightAslant');
        
        % Set the grouping variable for each tract group and the tracts that are not of interest (=0).
        if dorsal_idx(k) == 1 %tract of interest
            
            gp_tract(k) = 1;
            
        elseif ventral_idx(k) == 1
            
            gp_tract(k) = 2;
            
        elseif vertical_idx(k) == 1
            
            gp_tract(k) = 3;
            
        elseif vof_idx(k) == 1
            
            gp_tract(k) = 4;
            
        elseif aslant_idx(k) == 1
            
            gp_tract(k) = 5;
            
        else
            
            gp_tract(k) = 0;
            
        end
        
    end
    
    % Select the measurements of the tracts that I care about and convert all zeros to NaN.
    toi = wm_data(:, find(gp_tract ~= 0)); toi(toi==0) = NaN;
    
    % Update the tract list so that we can grab the tract name.
    wm_header = wm_header(gp_tract~=0);
    
    % Update the tract indexing to correspond with toi dimensions.
    gp_tract = gp_tract(find(gp_tract~=0));
    
    % Remove outliers: z-scores that are more than 3 standard deviations from the within-tract, within-group mean.
    % Find within-tract, within-group mean.
    m_wtwg1 = nanmean(toi(find(wm.gp_age == 1), :));
    m_wtwg2 = nanmean(toi(find(wm.gp_age == 2), :));
    m_wtwg3 = nanmean(toi(find(wm.gp_age == 3), :));
    
    % Find within-tract, within-group mean.
    std_wtwg1 = nanstd(toi(find(wm.gp_age == 1), :));
    std_wtwg2 = nanstd(toi(find(wm.gp_age == 2), :));
    std_wtwg3 = nanstd(toi(find(wm.gp_age == 3), :));
    
    % Get range of acceptable data.
    r_max_wtwg1 = m_wtwg1 + 3*std_wtwg1;
    r_min_wtwg1 = m_wtwg1 - 3*std_wtwg1;
    r_max_wtwg2 = m_wtwg2 + 3*std_wtwg2;
    r_min_wtwg2 = m_wtwg2 - 3*std_wtwg2;
    r_max_wtwg3 = m_wtwg3 + 3*std_wtwg3;
    r_min_wtwg3 = m_wtwg3 - 3*std_wtwg3;
    
    % Organize max and min into matrix to make indexing easier.
    r_max = cat(1, repmat(r_max_wtwg1, [size(toi(find(wm.gp_age==1)))]), repmat(r_max_wtwg2, size(toi(find(wm.gp_age==2)))), ...
        repmat(r_max_wtwg3, size(toi(find(wm.gp_age==3)))));
    r_min = cat(1, repmat(r_min_wtwg1, [size(toi(find(wm.gp_age==1)))]), repmat(r_min_wtwg2, size(toi(find(wm.gp_age==2)))), ...
        repmat(r_min_wtwg3, size(toi(find(wm.gp_age==3)))));
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Replace outliers with NaN.
    % Max
    if ~isempty(find(toi > r_max))
                disp(['Replaced ' num2str(numel(find(toi > r_max))) ' data points that were above 3 standard deviations of the within-tract, within-group mean with NaN.'])
                toi(find(toi > r_max)) = NaN;
    else
        disp('No data points were above 3 standard deviations of the within-tract, within-group mean.')
    end
    %Min
    if ~isempty(find(toi < r_min))
                disp(['Replaced ' num2str(numel(find(toi < r_min))) ' data points that were below 3 standard deviations of the within-tract, within-group mean with NaN.'])
                toi(find(toi < r_min)) = NaN;
    else
        disp('No data points were below 3 standard deviations of the within-tract, within-group mean.')
    end
    
    % Make hemisphere averaged column for each tract.
    toi_havg = cat(2, nanmean(toi(:, [1 11]), 2), nanmean(toi(:, [2 12]), 2), nanmean(toi(:, [3 13]), 2), nanmean(toi(:, [4 14]), 2), nanmean(toi(:, [5 15]), 2), ...
        nanmean(toi(:, [6 16]), 2), nanmean(toi(:, [7 17]), 2), nanmean(toi(:, [8 18]), 2), nanmean(toi(:, [9 19]), 2), nanmean(toi(:, [10 20]), 2));  
    toi_havg_header = erase(wm_header(1:10), 'left');
    
    % Output csv file for ANOVA in SPSS so that it is available, if needed. (Matlab doesn't handle Mixed Model ANOVAs well when the between-group variable is correlated with subID
    % (e.g., when between-group variable is something like age groups).
    t_out = array2table(cat(2, wm.subID, wm.gp_age, wm.gp_lit, wm.gp_vm, ...
        wm.gp_fm, wm.cov_age, wm.cov_sex, toi(:, 1:end), toi_havg(:, 1:end)), ...
        'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', wm_header{:}, toi_havg_header{:}});
    
    % Write.
    writetable(t_out, fullfile(rootDir, 'supportFiles', ['lwx_' wm_measure_here{w} '_singleshell_datareorg.csv']));
    
%     % Output z-scored file for SPSS.
    toi_z = (nanmean(toi, 1) - toi)./nanstd(toi, [], 1);
    toi_havg_z = (nanmean(toi_havg, 1) - toi_havg)./nanstd(toi_havg, [], 1);

%     temp = nanmean(toi(:, hv == 1), 2);
%     toi_meanhd_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
%     
%     temp = nanmean(toi(:, hv == 2), 2);
%     toi_meanhv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
%     
%     temp = nanmean(toi(:, hv == 3), 2);
%     toi_meanv_z = (nanmean(temp, 1) - temp)./nanstd(temp, [], 1); clear temp
    
    t_out_z = array2table(cat(2, wm.subID, wm.gp_age, wm.gp_lit, wm.gp_vm, ...
        wm.gp_fm, wm.cov_age, wm.cov_sex, toi_z(:, 1:end), toi_havg_z(:, 1:end)), ...
        'VariableNames', {'subID', 'group_age', 'group_lit', 'group_vm', 'group_fm', 'cov_age', 'cov_sex', wm_header{:}, toi_havg_header{:}});
    writetable(t_out_z, fullfile(rootDir, 'supportFiles', ['lwx_' wm_measure_here{w} '_singleshell_datareorg_z.csv']));
    
end


