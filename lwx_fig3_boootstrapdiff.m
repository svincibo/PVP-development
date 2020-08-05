clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

wm_measure = 'fa';

hemisphere = 'both'; %left, right , both

save_figures = 'yes';
alphastat = 0.05;

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
darkgray = [105 105 105]/255;
lightgray = [211 211 211]/255;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 100;
xtickvalues = [1 2 3 4 5];
xlim_lo = 0.5; xlim_hi = 5.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = 0.40; ylimhi = 0.55;

%% WHITE MATTER MEASURES

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
d = readtable(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure '_singleshell.csv']));

% Get index for outliers to be removed.
idx_keep = find(~ismember(d.subID, outlier));

% Remove outliers.
d = d(idx_keep, :);

% Get easy index for age group.
group = d.group_age;

% Plot means.
figure(1)
hold on;

% VOF
if strcmp(hemisphere, 'both')
    vof_child = nanmean(cat(2, d.leftVOF(group ~= 3), d.rightVOF(group ~= 3)), 2);
    vof_adult = nanmean(cat(2, d.leftVOF(group == 3), d.rightVOF(group == 3)), 2);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
elseif strcmp(hemisphere, 'left')
    vof_child = d.leftVOF(group ~= 3);
    vof_adult = d.leftVOF(group == 3);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
elseif strcmp(hemisphere, 'right')
    vof_child = d.rightVOF(group ~= 3);
    vof_adult = d.rightVOF(group == 3);
    vof_n_child = length(find(~isnan(vof_child))); vof_n_adult = length(find(~isnan(vof_adult)));
end

% Horizontal, Ventral
if strcmp(hemisphere, 'both')
    hv_child = nanmean(cat(2, d.leftILF(group ~= 3), d.rightILF(group ~= 3), ...
        d.leftIFOF(group ~= 3), d.rightIFOF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 3), d.rightILF(group == 3), ...
        d.leftIFOF(group == 3), d.rightIFOF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(hemisphere, 'left')
    hv_child = nanmean(cat(2, d.leftILF(group ~= 3), d.leftIFOF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 3), d.leftIFOF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(hemisphere, 'right')
    hv_child = nanmean(cat(2, d.rightILF(group ~= 3), d.rightIFOF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.rightILF(group == 3), d.rightIFOF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
end

% Vertical, Posterior
if strcmp(hemisphere, 'both')
    pv_child = nanmean(cat(2, d.leftMDLFang(group ~= 3), d.rightMDLFang(group ~= 3), ...
        d.leftMDLFspl(group ~= 3), d.rightMDLFspl(group ~= 3), ...
        d.leftTPC(group ~= 3), d.rightTPC(group ~= 3), ...
        d.leftpArc(group ~= 3), d.rightpArc(group ~= 3)), 2);
    pv_adult = nanmean(cat(2, d.leftMDLFang(group == 3), d.rightMDLFang(group == 3), ...
        d.leftMDLFspl(group == 3), d.rightMDLFspl(group == 3), ...
        d.leftTPC(group == 3), d.rightTPC(group == 3), ...
        d.leftpArc(group == 3), d.rightpArc(group == 3)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
elseif strcmp(hemisphere, 'left')
    pv_child = nanmean(cat(2, d.leftMDLFang(group ~= 3), d.leftMDLFspl(group ~= 3), d.leftTPC(group ~= 3), d.rightTPC(group ~= 3), d.leftpArc(group ~= 3)), 2);
    pv_adult = nanmean(cat(2, d.leftMDLFang(group == 3), d.leftMDLFspl(group == 3), d.leftTPC(group == 3), d.leftpArc(group == 3)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
elseif strcmp(hemisphere, 'right')
    pv_child = nanmean(cat(2, d.rightMDLFang(group ~= 3), d.rightMDLFspl(group ~= 3), d.rightTPC(group ~= 3), d.rightTPC(group ~= 3), d.rightpArc(group ~= 3)), 2);
    pv_adult = nanmean(cat(2, d.rightMDLFang(group == 3), d.rightMDLFspl(group == 3), d.rightTPC(group == 3), d.rightpArc(group == 3)), 2);
    pv_n_child = length(find(~isnan(pv_child))); pv_n_adult = length(find(~isnan(pv_adult)));
end

% Horizontal, Dorsal
if strcmp(hemisphere, 'both')
    hd_child = nanmean(cat(2, d.leftSLF1And2(group ~= 3), d.rightSLF1And2(group ~= 3), ...
        d.leftSLF3(group ~= 3), d.rightSLF3(group ~= 3)), 2);
    hd_adult = nanmean(cat(2, d.leftSLF1And2(group == 3), d.rightSLF1And2(group == 3), ...
        d.leftSLF3(group == 3), d.rightSLF3(group == 3)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
elseif strcmp(hemisphere, 'left')
    hd_child = nanmean(cat(2, d.leftSLF1And2(group ~= 3), d.leftSLF3(group ~= 3)), 2);
    hd_adult = nanmean(cat(2, d.leftSLF1And2(group == 3), d.leftSLF3(group == 3)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
elseif strcmp(hemisphere, 'right')
    hd_child = nanmean(cat(2, d.rightSLF1And2(group ~= 3), d.rightSLF3(group ~= 3)), 2);
    hd_adult = nanmean(cat(2, d.rightSLF1And2(group == 3), d.rightSLF3(group == 3)), 2);
    hd_n_child = length(find(~isnan(hd_child))); hd_n_adult = length(find(~isnan(hd_adult)));
end

% FAT
if strcmp(hemisphere, 'both')
    fat_child = nanmean(cat(2, d.leftAslant(group ~= 3), d.rightAslant(group ~= 3)), 2);
    fat_adult = nanmean(cat(2, d.leftAslant(group == 3), d.rightAslant(group == 3)), 2);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
elseif strcmp(hemisphere, 'left')
    fat_child = d.leftAslant(group ~= 3);
    fat_adult = d.leftAslant(group == 3);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
elseif strcmp(hemisphere, 'right')
    fat_child = d.rightAslant(group ~= 3);
    fat_adult = d.rightAslant(group == 3);
    fat_n_child = length(find(~isnan(fat_child))); fat_n_adult = length(find(~isnan(fat_adult)));
end

child_mean = [mean(vof_child) mean(hv_child) mean(pv_child) mean(hd_child) mean(fat_child)];
adult_mean = [mean(vof_adult) mean(hv_adult) mean(pv_adult) mean(hd_adult) mean(fat_adult)];

child_ci = 1.96*[std(vof_child)/sqrt(vof_n_child) std(hv_child)/sqrt(hv_n_child) std(pv_child)/sqrt(pv_n_child) std(hd_child)/sqrt(hd_n_child) std(fat_child)/sqrt(hd_n_child)];
adult_ci = 1.96*[std(vof_adult)/sqrt(vof_n_adult) std(hv_adult)/sqrt(hv_n_adult) std(pv_adult)/sqrt(pv_n_adult) std(hd_adult)/sqrt(hd_n_adult) std(fat_adult)/sqrt(hd_n_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

scatter(xval, child_mean, 'Marker', 'o', 'SizeData', markersize, 'MarkerFaceColor', lightgray, 'MarkerEdgeColor', lightgray)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', markersize, 'MarkerFaceColor', darkgray, 'MarkerEdgeColor', darkgray)

errorbar(xval, adult_mean, adult_ci, 'Color', darkgray, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
errorbar(xval, child_mean, child_ci, 'Color', lightgray, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
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
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adults'}, 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', fontsize);
legend('boxoff');

a.YLabel.String = 'Fractional Anisotropy (FA)';
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_' wm_measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_' wm_measure '_' hemisphere]), '-depsc')
    
end

hold off;


% Plot bootstrapped difference between children and adults.
figure(2)
hold on;

% Bootstrap to get diff and error bars because unequal sample sizes.
[vof_diff, vof_ci] = bootstrap_diff_unequalsamplesizes(vof_child, vof_adult, alphastat);
[hv_diff, hv_ci] = bootstrap_diff_unequalsamplesizes(hv_child, hv_adult, alphastat);
[pv_diff, pv_ci] = bootstrap_diff_unequalsamplesizes(pv_child, pv_adult, alphastat);
[hd_diff, hd_ci] = bootstrap_diff_unequalsamplesizes(hd_child, hd_adult, alphastat);
[fat_diff, fat_ci] = bootstrap_diff_unequalsamplesizes(fat_child, fat_adult, alphastat);

diff = [vof_diff hv_diff pv_diff hd_diff fat_diff];
ci = [vof_ci hv_ci pv_ci hd_ci fat_ci];

xval = linspace(1, length(diff), length(diff));

scatter(xval, diff, 'Marker', 'o', 'SizeData', markersize, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
errorbar(xval, diff, ci, 'Color', darkgray, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xlabels = {'VOF', 'Ventral Horizontal', 'Posterior Vertical', 'Dorsal Horizontal', 'FAT'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
ylimlo = -0.1; ylimhi = 0.15;
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

plot([xlim_lo xlim_hi], [0 0], 'k:')

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'Difference in FA: Adult - Children';
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-dpng')
    print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_diff_' wm_measure '_' hemisphere]), '-depsc')
    
end

hold off;

%% =========================================================================================== %%

function [diff_out, ci_out] = bootstrap_diff_unequalsamplesizes(group1, group2, alpha)

% get sample sizes of each group
n1 = length(group1);
n2 = length(group2);

% get distrution of differences
for r = 1:10000

    diff_real(r) = group2(randi(n2)) - group1(randi(n1));
        
end

diff_out = nanmean(diff_real);
ci_temp = prctile(diff_real, [100*alpha/2, 100*(1-alpha/2)]);
ci_out = (ci_temp(2) - ci_temp(1))/2;

end





