clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

measure = 'volume'; % 'fa', 'volume', 'gmd', 'md'

hemisphere = 'both'; %left, right , both

save_figures = 'yes';
alphastat = 0.66; % to return 1 SD, for 95% CI use .05

color_adults = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
color_children = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed.
    %     outlier = [108 126 318];%
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    
    
else
    
    outlier = [];
    
end

% Set up plot and measure-specific details.
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4];
xlim_lo = min(xtickvalues)-0.5; xlim_hi = max(xtickvalues)+0.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

if strcmp(measure, 'fa')
    ylimlo = 0.12; ylimhi = 0.25;
elseif strcmp(measure, 'volume')
    ylimlo = 0.001; ylimhi = .005;
%             ylimlo = 0.0000006; ylimhi = .0000025;
elseif strcmp(measure, 'snr')
    ylimlo = 5; ylimhi = 32;
    elseif strcmp(measure, 'gmd')
    ylimlo = 0.8; ylimhi = .95;
end

%% WHITE MATTER MEASURES

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
d = readtable(fullfile(rootDir, 'supportFiles', ['LWX_CM_data_forSPSS_' measure '_singleshell.csv']));

% Get index for outliers to be removed.
idx_keep = find(~ismember(d.subID, outlier));

% Remove outliers.
d = d(idx_keep, :);

% Get easy index for age group.
group = d.group_age;

occipital_child = d.occipital(group ~=3);
ventral_child = d.ventral(group ~= 3);
parietal_child = d.parietal(group ~= 3);
frontal_child = d.frontal(group ~= 3);

occipital_adult = d.occipital(group == 3);
ventral_adult = d.ventral(group == 3);
parietal_adult = d.parietal(group == 3);
frontal_adult = d.frontal(group == 3);
    
%% Plot means including individual data.

figure(1)
hold on;

coloralpha = .1;

% Means (do this first for legend and then second to keep it on top layer).
child_mean = [nanmean(occipital_child) nanmean(ventral_child) mean(parietal_child) mean(frontal_child)];
adult_mean = [nanmean(occipital_adult) nanmean(ventral_adult) mean(parietal_adult) mean(frontal_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

% Individual data points for children.
scatter(repmat(1, size(occipital_child)), occipital_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(ventral_child)), ventral_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(parietal_child)), parietal_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(frontal_child)), frontal_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Individual data points for adults.
scatter(repmat(1, size(occipital_adult)), occipital_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(2, size(ventral_adult)), ventral_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(3, size(parietal_adult)), parietal_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)
scatter(repmat(4, size(frontal_adult)), frontal_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Means (second time to put it on top layer).
child_mean = [nanmean(occipital_child) mean(ventral_child) mean(parietal_child) mean(frontal_child)];
adult_mean = [nanmean(occipital_adult) mean(ventral_adult) mean(parietal_adult) mean(frontal_adult)];

child_sd = [std(occipital_child, 'omitnan') std(ventral_child) std(parietal_child) std(frontal_child)];
adult_sd = [std(occipital_adult, 'omitnan') std(ventral_adult) std(parietal_adult) std(frontal_adult)];

% child_ci = 1.96*[std(occipital_child)/sqrt(vof_n_child) std(ventral_child)/sqrt(hv_n_child) std(parietal_child)/sqrt(pv_n_child) std(frontal_child)/sqrt(hd_n_child) std(fat_child)/sqrt(hd_n_child)];
% adult_ci = 1.96*[std(occipital_adult)/sqrt(vof_n_adult) std(ventral_adult)/sqrt(hv_n_adult) std(parietal_adult)/sqrt(pv_n_adult) std(frontal_adult)/sqrt(hd_n_adult) std(fat_adult)/sqrt(hd_n_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

errorbar(xval, adult_mean, adult_sd, 'Color', color_adults, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
errorbar(xval, child_mean, child_sd, 'Color', color_children, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'Occipital', 'Ventral', 'Parietal', 'Frontal'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
if strcmp(measure, 'volume')
    yax.TickLabels = {num2str(ylimlo, '%2.1e'), num2str((ylimlo+ylimhi)/2, '%2.1e'), num2str(ylimhi, '%2.1e')};
else
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
end
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

% legend({'Children', 'Adults'}, 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', fontsize);
% legend('boxoff');

if strcmp(measure, 'fa')
    a.YLabel.String = 'Fractional Anisotropy (FA)';
elseif strcmp(measure, 'volume')
    a.YLabel.String = 'Gray Matter Volume (GMV) Proportion';
    elseif strcmp(measure, 'gmd')
    a.YLabel.String = 'Gray Matter Density (GMD)';
elseif strcmp(measure, 'snr')
    a.YLabel.String = 'tSNR';
end

a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_cm_ind_' measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_cm_ind_' measure '_' hemisphere]), '-depsc')
    
end

hold off;

% % Plot bootstrapped difference between children and adults.
% figure(3)
% hold on;
%
% % Bootstrap to get diff and error bars because unequal sample sizes.
% [vof_diff, vof_ci] = bootstrap_diff_unequalsamplesizes(occipital_child, occipital_adult, alphastat);
% [hv_diff, hv_ci] = bootstrap_diff_unequalsamplesizes(ventral_child, ventral_adult, alphastat);
% [pv_diff, pv_ci] = bootstrap_diff_unequalsamplesizes(parietal_child, parietal_adult, alphastat);
% [hd_diff, hd_ci] = bootstrap_diff_unequalsamplesizes(frontal_child, frontal_adult, alphastat);
% [fat_diff, fat_ci] = bootstrap_diff_unequalsamplesizes(fat_child, fat_adult, alphastat);
%
% diff = [vof_diff hv_diff pv_diff hd_diff fat_diff];
% sd = [vof_ci hv_ci pv_ci hd_ci fat_ci]; % returns 1 SD when alphastat = .66
%
% xval = linspace(1, length(diff), length(diff));
%
% scatter(xval, diff, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
% errorbar(xval, diff, sd, 'Color', 'k', 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
%
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [xlim_lo xlim_hi];
% xax.TickValues = xtickvalues;
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
% xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
%
% % yaxis
% ylimlo = -0.02; ylimhi = 0.10;
% yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
% yax.TickValues = [ylimlo 0 (ylimlo+ylimhi)/2 ylimhi];
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = {num2str(ylimlo, '%2.2f'), '0', num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
%
% plot([xlim_lo xlim_hi], [0 0], 'k:')
%
% % general
% a = gca;
% %     a.TitleFontWeight = 'normal';
% box off
%
% a.YLabel.String = 'Modularity Index (MI)';
% a.YLabel.FontSize = fontsize;
% a.YLabel.FontAngle = fontangle;
% pbaspect([1 1 1])
%
% %     pos=get(gca,'Position');
% %     pos1=pos-[0 .02 0 0];
% %     set(gca,'Position', pos1);
%
% % Write.
% if strcmp(save_figures, 'yes')
%
%     print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-dpng')
%     print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-depsc')
%
% end
%
% hold off;
%
% %% =========================================================================================== %%
%
% function [diff_out, ci_out] = bootstrap_diff_unequalsamplesizes(group1, group2, alpha)
%
% % get sample sizes of each group
% n1 = length(group1);
% n2 = length(group2);
%
% % get distrution of differences
% for r = 1:10000
%
%     diff_real(r) = group2(randi(n2)) - group1(randi(n1));
%
% end
%
% diff_out = nanmean(diff_real);
% ci_temp = prctile(diff_real, [100*alpha/2, 100*(1-alpha/2)]); % treats NaNs as missing values and removes them
% ci_out = (ci_temp(2) - ci_temp(1));
%
% end
%




