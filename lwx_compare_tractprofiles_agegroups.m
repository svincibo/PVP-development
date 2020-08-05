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
streamline_min = 100;
alphastat = 0.01;
binsize = 5;

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray
hold on;
linewidth = 1.5;
linestyle = 'none';
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;
ylimlo = 0.20; ylimhi = 0.70;

save_figures = 'yes';

% Load wm data from lwx_qa_tractstats.m.
load([rootDir 'supportFiles/LWX_data_streamlinecount.mat']);

% Load streamline count data from lwx_datacat.m.
load([rootDir 'supportFiles/LWX_data_fa_singleshell.mat']);

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - conservative removal.
    %     outlier = [108 126 214 318];%
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 214, major motion artifacts, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Identify outliers to be removed - liberal removal.
    outlier = [108 116 119 125 126 206 214 303 317 318];
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
    
    % Preallocate array for the z scores for yc vs oc.
    z_keep = NaN(size(list_tract, 1), (size(m_wm, 1)-1)/binsize);
    z_keep_tractname = cell(size(z_keep));
    
    % Plot tract profiles for each tract.
    for k = 1:size(list_tract, 1)
        
        % Only plot for tracts of interest and non-empty tracts.
        if strcmp(list_tract{k},'leftSLF1And2') || strcmp(list_tract{k}, 'rightSLF1And2') ...
                || strcmp(list_tract{k},'leftSLF3') || strcmp(list_tract{k}, 'rightSLF3') ...
                || strcmp(list_tract{k}, 'leftILF') || strcmp(list_tract{k}, 'rightILF') ...
                || strcmp(list_tract{k}, 'leftIFOF') || strcmp(list_tract{k}, 'rightIFOF') ...
                || strcmp(list_tract{k}, 'leftAslant') || strcmp(list_tract{k}, 'rightAslant') ...
                || strcmp(list_tract{k}, 'leftTPC') || strcmp(list_tract{k}, 'rightTPC') ...
                || strcmp(list_tract{k}, 'leftpArc') || strcmp(list_tract{k}, 'rightpArc') ...
                || strcmp(list_tract{k}, 'leftMDLFspl') || strcmp(list_tract{k}, 'rightMDLFspl') ...
                || strcmp(list_tract{k}, 'leftVOF') || strcmp(list_tract{k}, 'rightVOF') ...
                || strcmp(list_tract{k}, 'leftMDLFang') || strcmp(list_tract{k}, 'rightMDLFang') ...
                && ~strcmp(list_tract{k}, 'empty')
            
            disp(list_tract{k})
            
            % Find entries that are for this tract.
            t_idx = strcmp(tract, list_tract{k});
            
            % Open a new figure for this tract.
            fcount = fcount + 1;
            figure(fcount)
            hold on;
            
            yc_count = 0; oc_count = 0; a_count = 0;
            for s = 1:length(subID)
                
                % Identify subjects who have less than 100 streamlines for this track.
                col_idx = find(strcmp(streamlinecounts.Properties.VariableNames, list_tract{k}));
                subID_stream = streamlinecounts.subID(find(table2array(streamlinecounts(:, col_idx)) < streamline_min));
                
                % Only include subjects who are not outliers and who have at least 100 streamlines for this tract.
                if subID(s)~=0 && ~ismember(subID(s), outlier) && ~ismember(subID(s), subID_stream)
                    
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
                        
%                         hold on;
                        
                        clear t_temp;
                        
                    end % if subID{s} ~= 0
                    
%                     % Only include subjects who are not outliers and who
%                     % have at least 100 streamlines for this tract. Set
%                     % cells with less than 100 streamlines to NaN.
%                 elseif subID(s)~=0 && ~ismember(subID(s), outlier) && ismember(subID(s), subID_stream)
                    
                end %if exist
                
            end %sub
            
            multcomp_correction = length(1:binsize:size(yc, 1));
            %% Perform bootstrap test between child groups.
            b_count = 0;
            for b = 1:binsize:size(yc, 1)-binsize
                
                b_count = b_count + 1;
                
                % Subset data for this testing bin.
                youngerchildren = yc(b:b+binsize, :);
                olderchildren = oc(b:b+binsize, :);
                
                % Bootstrap to get error bars because unequal sample sizes.
                for r = 1:10000
                    
                    this_c = randsample(youngerchildren(:), size(youngerchildren, 2), true);
                    this_a = randsample(olderchildren(:), size(olderchildren, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_real(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % creat the null distribution
                null_dis = cat(2, youngerchildren, olderchildren);
                
                for r = 1:10000
                    
                    this_c = randsample(null_dis(:), size(youngerchildren, 2), true);
                    this_a = randsample(null_dis(:), size(olderchildren, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_null(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % Confidence interval.
                ci = prctile(diff_real, [100*alphastat/2, 100*(1-alphastat/2)]);
                
                % p-value
                mu_diff_real = nanmean(diff_real);
                mu_diff_null = nanmean(diff_null);
                sigma_diff_null = nanstd(diff_null);
                
                z = (mu_diff_real - mu_diff_null)./sigma_diff_null;
                p = 1-normcdf(abs(z), 0, 1);
                
%                 z_keep(k, b_count) = z;
%                 z_keep_tractname{k, b_count} = list_tract{k};
%                 
                if p <= 0.05
                    
%                     patch([b b+binsize-.25 b+binsize-.25 b], [0.20 0.20 0.70, 0.70], [0.9290 0.6940 0.1250], 'FaceColor', [0.9290 0.6940 0.1250], ...
%                         'EdgeColor', [0.9290 0.6940 0.1250], 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2) %burnt yellow
                    
                    if p <= .001/multcomp_correction
                        
                        scatter((b+b+binsize)/2-4, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        scatter((b+b+binsize)/2, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        scatter((b+b+binsize)/2+4, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        
                    elseif p <= 0.01/multcomp_correction
                        
                        scatter((b+b+binsize)/2-2, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        scatter((b+b+binsize)/2+2, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        
                    else
                        
                        scatter((b+b+binsize)/2, 0.60, 'Marker', '*', 'MarkerEdgeColor', [0.9290 0.6940 0.1250], 'SizeData', 100)
                        
                    end
                    
                end % end if p<=0.05
                
            end % for binsize
                     
            %% Perform bootstrap test between older children and adult groups.
            b_count = 0;
            for b = 1:binsize:size(yc, 1)-binsize
                
                b_count = b_count + 1;
                
                % Subset data for this testing bin.
                olderchildren = oc(b:b+binsize, :);
                adults = a(b:b+binsize, :);

                % Bootstrap to get error bars because unequal sample sizes.
                for r = 1:10000
                    
                    this_c = randsample(olderchildren(:), size(olderchildren, 2), true);
                    this_a = randsample(adults(:), size(adults, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_real(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % creat the null distribution
                null_dis = cat(2, olderchildren, adults);
                
                for r = 1:10000
                    
                    this_c = randsample(null_dis(:), size(olderchildren, 2), true);
                    this_a = randsample(null_dis(:), size(adults, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_null(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % Confidence interval.
                ci = prctile(diff_real, [100*alphastat/2, 100*(1-alphastat/2)]);
                
                % p-value
                mu_diff_real = nanmean(diff_real);
                mu_diff_null = nanmean(diff_null);
                sigma_diff_null = nanstd(diff_null);
                
                z = (mu_diff_real - mu_diff_null)./sigma_diff_null;
                p = 1-normcdf(abs(z), 0, 1);
                
%                 z_keep(k, b_count) = z;
%                 z_keep_tractname{k, b_count} = list_tract{k};
                
                if p <= 0.05
                    
%                     patch([b b+binsize-.25 b+binsize-.25 b], [0.20 0.20 0.70, 0.70], [0 1 1], 'FaceColor', [0 1 1], ...
%                         'EdgeColor', [0 1 1], 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2) %cyan
                    
                    if p <= .001/multcomp_correction
                        
                        scatter((b+b+binsize)/2-4, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2+4, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                        
                    elseif p <= 0.01/multcomp_correction
                        
                        scatter((b+b+binsize)/2-2, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2+2, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                        
                    else
                        
                        scatter((b+b+binsize)/2, 0.625, 'Marker', '*', 'MarkerEdgeColor', [0,206,209]/255, 'SizeData', 100)
                                        
%                         z_keep(k, b_count) = z;
%                         z_keep_tractname{k, b_count} = list_tract{k};
                
                    end
                    
                end % end if p<=0.05
                
            end % for binsize
                 
            %% Perform bootstrap test between *all* children and adult groups.
            b_count = 0;
            for b = 1:binsize:size(yc, 1)-binsize
                
                b_count = b_count + 1;

                % Subset data for this testing bin.
                children = cat(2, yc(b:b+binsize, :), oc(b:b+binsize, :));
                adults = a(b:b+binsize, :);

                % Bootstrap to get error bars because unequal sample sizes.
                for r = 1:10000
                    
                    this_c = randsample(children(:), size(children, 2), true);
                    this_a = randsample(adults(:), size(adults, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_real(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % creat the null distribution
                null_dis = cat(2, children, adults);
                
                for r = 1:10000
                    
                    this_c = randsample(null_dis(:), size(children, 2), true);
                    this_a = randsample(null_dis(:), size(adults, 2), true);
                    
                    this_mean_c = nanmean(this_c);
                    this_mean_a = nanmean(this_a);
                    
                    diff_null(r) = this_mean_a - this_mean_c;
                    
                    clear this_c this_a this_mean_c this_mean_a
                    
                end
                
                % Confidence interval.
                ci = prctile(diff_real, [100*alphastat/2, 100*(1-alphastat/2)]);
                
                % p-value
                mu_diff_real = nanmean(diff_real);
                mu_diff_null = nanmean(diff_null);
                sigma_diff_null = nanstd(diff_null);
                
                z = (mu_diff_real - mu_diff_null)./sigma_diff_null;
                p = 1-normcdf(abs(z), 0, 1);
                
%                 z_keep(k, b_count) = z;
%                 z_keep_tractname{k, b_count} = list_tract{k};
                
                if p <= 0.05
                    
%                     patch([b b+binsize-.25 b+binsize-.25 b], [0.20 0.20 0.70, 0.70], [0 1 1], 'FaceColor', [0 1 1], ...
%                         'EdgeColor', [0 1 1], 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2) %cyan
                    
                    if p <= .001/multcomp_correction
                        
                        scatter((b+b+binsize)/2-4, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2+4, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                        
                    elseif p <= 0.01/multcomp_correction
                        
                        scatter((b+b+binsize)/2-2, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                        scatter((b+b+binsize)/2+2, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                        
                    else
                        
                        scatter((b+b+binsize)/2, 0.65, 'Marker', '*', 'MarkerEdgeColor', [60, 179, 113]/255, 'SizeData', 100)
                                        
                        z_keep(k, b_count) = z;
                        z_keep_tractname{k, b_count} = list_tract{k};
                
                    end
                    
                end % end if p<=0.05
                
            end % for binsize
            
            % Plot means and 95% confidence intervals (calculated from
            % standard error: 1.96*SE). Use nanmean/nanstd because one oc subject is missing TPC.
            plot(nanmean(yc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', yc_color(1:3))
            hi = nanmean(yc, 2) + 1.96*nanstd(yc, 0, 2)/sqrt(size(~isnan(yc), 2)); lo = nanmean(yc, 2) - 1.96*nanstd(yc, 0, 2)/sqrt(size(~isnan(yc), 2)); x = (1:size(nanmean(yc, 2),1))';
            hp1 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], yc_color(1:3));
            set(hp1, 'facecolor', yc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(nanmean(oc, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', oc_color(1:3))
            hi = nanmean(oc, 2) + 1.96*nanstd(oc, 0, 2)/sqrt(size(~isnan(oc), 2)); lo = nanmean(oc, 2) - 1.96*nanstd(oc, 0, 2)/sqrt(size(~isnan(oc), 2)); x = (1:size(nanmean(oc, 2),1))';
            hp2 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], oc_color(1:3));
            set(hp2, 'facecolor', oc_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            plot(nanmean(a, 2), 'LineWidth', 3, 'LineStyle', '-', 'Color', a_color(1:3))
            hi = nanmean(a, 2) + 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)); lo = nanmean(a, 2) - 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)); x = (1:size(nanmean(a, 2),1))';
            hp3 = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], a_color(1:3));
            set(hp3, 'facecolor', a_color(1:3), 'edgecolor', 'none', 'facealpha', .2);
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 160];
            xax.TickValues = [0 80 160];
            xax.TickLabels = {'20', '100', '180'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylimlo ylimhi];
            yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylimlo, '%1.2f'), num2str((ylimlo+ylimhi)/2, '%1.2f'), num2str(ylimhi, '%1.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            %     a.TitleFontWeight = 'normal';
            box off
            
            % legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
            % legend box off
            
            title(list_tract{k})
            g.XLabel.String = 'Location along tract';
            g.XLabel.FontSize = fontsize;
            g.XLabel.FontAngle = fontangle;
            
            g.YLabel.String = 'Fractional Anisotropy (FA)';
            g.YLabel.FontSize = fontsize;
            pbaspect([1 1 1])
            
            print(fullfile(rootDir, 'plots-singleshell', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_tractprofiles_' wm_measure '_' list_tract{k}]), '-depsc')
            
            hold off;
            
            % Open a new figure for the mean plot.
            fcount = fcount + 1;
            figure(fcount)
            
            gr = [1 2 3];
            gscatter([0.5 1.5 2.5], [nanmean(yc, 'all') nanmean(oc, 'all') nanmean(a, 'all')], gr, cat(1, yc_color, oc_color, a_color))
            hold on;
            plot([0.5 0.5], [nanmean(yc, 2) + 1.96*nanstd(yc, 0, 2)/sqrt(size(~isnan(yc), 2)) nanmean(yc, 2) - 1.96*nanstd(yc, 0, 2)/sqrt(size(~isnan(yc), 2))], 'LineStyle', '-', 'Color', yc_color(1:3))
            plot([1.5 1.5], [nanmean(oc, 2) + 1.96*nanstd(oc, 0, 2)/sqrt(size(~isnan(oc), 2)) nanmean(oc, 2) - 1.96*nanstd(oc, 0, 2)/sqrt(size(~isnan(oc), 2))], 'LineStyle', '-', 'Color', oc_color(1:3))
            plot([2.5 2.5], [nanmean(a, 2) + 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2)) nanmean(a, 2) - 1.96*nanstd(a, 0, 2)/sqrt(size(~isnan(a), 2))], 'LineStyle', '-', 'Color', a_color(1:3))
            
            legend off
            
            % xaxis
            xax = get(gca, 'xaxis');
            xax.Limits = [0 3];
            xax.TickValues = [0.5 1.5 2.5];
            xax.TickLabels = {'Young Children', 'Older Children', 'Adults'};
            xax.TickDirection = 'out';
            xax.TickLength = [xticklength xticklength];
            xax.FontName = fontname;
            xax.FontSize = fontsize;
            xax.FontAngle = fontangle;
            
            % yaxis
            yax = get(gca,'yaxis');
            yax.Limits = [ylimlo ylimhi];
            yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
            yax.TickDirection = 'out';
            yax.TickLabels = {num2str(ylimlo, '%1.2f'), num2str((ylimlo+ylimhi)/2, '%1.2f'), num2str(ylimhi, '%1.2f')};
            yax.FontName = fontname;
            yax.FontSize = fontsize;
            
            % general
            g = gca;
            %     a.TitleFontWeight = 'normal';
            box off
            
            % legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
            % legend box off
            
            title(list_tract{k})
            
            g.YLabel.String = 'Fractional Anisotropy (FA)';
            g.YLabel.FontSize = fontsize;
            pbaspect([1 1 1])
            
            print(fullfile(rootDir, 'plots-singleshell', ['plot_meancomparison_singleshell_' wm_measure '_' list_tract{k}]), '-dpng')
            print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_meancomparison_singleshell_' wm_measure '_' list_tract{k}]), '-depsc')
            
            hold off;
            
            clear a oc yc 
            
        end % if toi
        
    end %tract
    
end % for w

% Get only non-empty cells.
figure
hold on;
idx = find(~isnan(z_keep));
z_new = z_keep(idx);
z_new_tractname = z_keep_tractname(idx);
tn = unique(z_new_tractname);

for r = 1:length(tn)
   
    idx_temp = find(strcmp(z_new_tractname, tn(r)));
    
    scatter(z_new(idx_temp), r*ones(size(z_new(idx_temp))), '*k', 'SizeData', 100);
    
    clear idx_temp
    
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [1.5 3];
xax.TickValues = [1.5 2.25 3];
xax.TickLabels = {'1.5', '2.25', '3.00'};
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
ylimlo = 0; ylimhi = 20;
yax.Limits = [ylimlo ylimhi];
yax.TickValues = 1:19;
yax.TickDirection = 'out';
yax.TickLabels = tn;
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
g = gca;
%     a.TitleFontWeight = 'normal';
box off

% legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'southeast');
% legend box off

g.XLabel.String = 'Effect Size (dprime)';
g.XLabel.FontSize = fontsize;

pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', ['plot_zhist_singleshell_' wm_measure]), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_zhist_singleshell_' wm_measure]), '-depsc')

hold off;
