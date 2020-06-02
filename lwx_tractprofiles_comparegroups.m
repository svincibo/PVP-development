% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

w_measures = {'fa'}; %, 'md'};

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

save_figures = 'yes';

% Load wm data from lwx_datacat.m.
load([rootDir 'supportFiles/LWX_data_fa_singleshell.mat']);

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - conservative removal.
%     outlier = [108 126 212 214 318];%
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
    
else
    
    outlier = [];
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

fcount = 0;
for w = 1:length(w_measures)
    
    wm_measure = w_measures{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    % Load in each tract's tractography measures for this subject.
    sub_count = 0;
    for i = 1:size(grp_contents, 1)
        
        % Display current sub ID.
        disp(grp_contents(i).name)
        
        % Update subject counter for when not all subjects are used/needed.
        sub_count = sub_count + 1;
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractprofiles = dir(fullfile(grp_contents(i).folder, grp_contents(i).name, 'dt-neuro-tractprofile*', 'profiles', '*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractprofiles = sub_contents_tractprofiles(arrayfun(@(x) x.name(1), sub_contents_tractprofiles) ~= '.');
        
        for j = 1:size(sub_contents_tractprofiles)
            
            % Preallocate based on number of subjects(size(grp_contents)) and number of tracts (size(sub_contents...)).
            if i == 1 && j == 1
                
                tract = {}; m = NaN(size(grp_contents, 1), size(sub_contents_tractprofiles, 1));
                
            end
            
            % Read in data for this subject and this tract.
            data_temp = readtable(fullfile(sub_contents_tractprofiles(j).folder, sub_contents_tractprofiles(j).name));
            
            % Get middle 80%.
            start = size(data_temp, 1)*.1;
            stop = size(data_temp, 1)*.9;
            
            % Read in mean WM measure.
            if strcmp(wm_measure, 'ad')
                
                m_wm(:, j, sub_count) = data_temp.ad_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.ad_2(start:stop);
                
            elseif strcmp(wm_measure, 'fa')
                
                m_wm(:, j, sub_count) = data_temp.fa_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.fa_2(start:stop);
                
            elseif strcmp(wm_measure, 'md')
                
                m_wm(:, j, sub_count) = data_temp.md_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.md_2(start:stop);
                
            elseif strcmp(wm_measure, 'rd')
                
                m_wm(:, j, sub_count) = data_temp.rd_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.rd_2(start:stop);
                
            elseif strcmp(wm_measure, 'icvf')
                
                m_wm(:, j, sub_count) = data_temp.icvf_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.icvf_2(start:stop);
                
            elseif strcmp(wm_measure, 'isovf')
                
                m_wm(:, j, sub_count) = data_temp.isovf_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.isovf_2(start:stop);
                
            elseif strcmp(wm_measure, 'od')
                
                m_wm(:, j, sub_count) = data_temp.od_1(start:stop);
                sd_wm(:, j, sub_count) = data_temp.od_2(start:stop);
                
            end
            
            % Grab tract name for grouping variable.
            tract(:, j, sub_count) = repmat({sub_contents_tractprofiles(j).name(1:end-13)}, 161, 1);
            
            % Grab subID.
            sub(:, j, sub_count) = repmat(str2num(grp_contents(i).name(end-2:end)), 161, 1);
            
            clear data_temp
            
        end % sub_contents
        
    end % group_contents
    
    % Find empty cells and fill with 'empty'.
    t = find(cellfun(@isempty,tract));
    tract(t) = {'empty'};
    
    % Get group indices for tract names.
    G = findgroups(tract(:));
    for i = 1:length(sub)
        if sub(i) < 200
            group(i) = 1;
        elseif sub(i) < 300
            group(i) = 2;
        else
            group(i) = 3;
        end
    end
    group = group';
    
    % Get a list of unique tract names.
    list_tract = unique(tract);
    
    % Get a list of unique sub IDs.
    subID = unique(sub);
    
    % Plot tract profiles for each tract.
    for k = 1:size(list_tract, 1)
        
        % Only plot for tracts of interest and non-empty tracts.
        if strcmp(list_tract{k},'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
                || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
                || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
                || strcmp(list_tract{k}, 'leftArc') || strcmp(list_tract{k}, 'rightArc') ...
                || strcmp(list_tract{k}, 'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3') ...
                || strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
                || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
                || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
                || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
                || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
                || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang') ...
                && ~strcmp(list_tract{k}, 'empty')
            
            % Find entries that are for this tract.
            t_idx = strcmp(tract, list_tract{k});
            
            % Open a new figure for this tract.
            fcount = fcount + 1;
            figure(fcount)
            
            yc_count = 0; oc_count = 0; a_count = 0;
            for s = 1:length(subID)
                
                % Only include subjects who are not outliers.
                if subID(s)~=0 && ~ismember(subID(s), outlier)
                    
                    % Find entries that are for this subject.
                    s_idx = sub == subID(s);
                    
                    % Subset the thing so that we only plot for this tract and for this subject.
                    t_temp = m_wm(find(t_idx == 1 & s_idx == 1));
                    
                    if ~isempty(t_temp)
                        
                        % Code the plot for subject and keep data for inspection (yc, oc, a).
                        if subID(s) < 200
                            
                            yc_count = yc_count + 1;
                            
                            % Young child.
                            plot(t_temp, 'LineStyle', '-', 'Color', [yc_color .2])
                            
                            % Collect.
                            yc(:, yc_count) = t_temp;
                            
                        elseif subID(s) < 300
                            
                            oc_count = oc_count + 1;
                            
                            % Older child.
                            plot(t_temp, 'LineStyle', '-', 'Color', [oc_color .2])
                            
                            % Collect.
                            oc(:, oc_count) = t_temp;
                            
                        else
                            
                            a_count = a_count + 1;
                            
                            % Adult.
                            plot(t_temp, 'LineStyle', '-', 'Color', [a_color .2])
                            
                            % Collect.
                            a(:, a_count) = t_temp;
                            
                        end
                        hold on;
                        
                    end % if subID{s} ~= 0
                    
                end %if exist
                
            end %sub
            
            % Plot means and standard deviations.
            plot(mean(yc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', yc_color(1:3))
            hi = mean(yc, 2) + std(yc, 0, 2); lo = mean(yc, 2) - std(yc, 0, 2); x = (1:size(mean(yc, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], yc_color(1:3));
            set(hp1, 'facecolor', yc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(mean(oc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', oc_color(1:3))
            hi = mean(oc, 2) + std(oc, 0, 2); lo = mean(oc, 2) - std(oc, 0, 2); x = (1:size(mean(oc, 2),1))';
            hp2 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], oc_color(1:3));
            set(hp2, 'facecolor', oc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(mean(a, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', a_color(1:3))
            hi = mean(a, 2) + std(a, 0, 2); lo = mean(a, 2) - std(a, 0, 2); x = (1:size(mean(a, 2),1))';
            hp3 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], a_color(1:3));
            set(hp3, 'facecolor', a_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            ylabel(wm_measure);
            xlabel(list_tract{k});
            xlim([0 160])
            box off;
            
            print(fullfile(rootDir, 'plots', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-depsc')
            
            hold off;
            
        end % if toi
        
    end %tract
    
end % for w




