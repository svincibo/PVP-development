clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

wm_measure = 'fa';

tract = 'mdlfspl'; % ilf, ifof, ilfifof, slf1and2, slf3, slf1and2slf3

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
markersize = 150;
xtickvalues = [1];
xlim_lo = 0.5; xlim_hi = 1.5;
fontname = 'Arial';
fontsize = 20;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;

ylimlo = 0.35; ylimhi = 0.60;

%% WHITE MATTER MEASURES

% Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
d = readtable(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure '_singleshell.csv']));

% Get index for outliers to be removed.
idx_keep = find(~ismember(d.subID, outlier));

% Remove outliers.
d = d(idx_keep, :);

% Get easy index for age group.
group = d.group_age;

% VOF
if strcmp(tract, 'ilf')
    hv_child = nanmean(cat(2, d.leftILF(group ~= 3), d.rightILF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 3), d.rightILF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'ifof')
    hv_child = nanmean(cat(2, d.leftIFOF(group ~= 3), d.rightIFOF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftIFOF(group == 3), d.rightIFOF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'ilfifof')
    hv_child = nanmean(cat(2, d.leftILF(group ~= 3), d.rightILF(group ~= 3), ...
        d.leftIFOF(group ~= 3), d.rightIFOF(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftILF(group == 3), d.rightILF(group == 3), ...
        d.leftIFOF(group == 3), d.rightIFOF(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'slf1and2')
    hv_child = nanmean(cat(2, d.leftSLF1And2(group ~= 3), d.rightSLF1And2(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftSLF1And2(group == 3), d.rightSLF1And2(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'slf3')
    hv_child = nanmean(cat(2, d.leftSLF3(group ~= 3), d.rightSLF3(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftSLF3(group == 3), d.rightSLF3(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'slf1and2slf3')
    hv_child = nanmean(cat(2, d.leftSLF1And2(group ~= 3), d.rightSLF1And2(group ~= 3), ...
        d.leftSLF3(group ~= 3), d.rightSLF3(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftSLF1And2(group == 3), d.rightSLF1And2(group == 3), ...
        d.leftSLF3(group == 3), d.rightSLF3(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'pArctpcmdlfangmdlfspl')
    hv_child = nanmean(cat(2, d.leftpArc(group ~= 3), d.rightpArc(group ~= 3), ...
        d.leftTPC(group ~= 3), d.rightTPC(group ~= 3), ...
        d.leftMDLFang(group ~= 3), d.rightMDLFang(group ~= 3), ...
        d.leftMDLFspl(group ~= 3), d.rightMDLFspl(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftpArc(group == 3), d.rightpArc(group == 3), ...
        d.leftTPC(group == 3), d.rightTPC(group == 3), ...
        d.leftMDLFang(group == 3), d.rightMDLFang(group == 3), ...
        d.leftMDLFspl(group == 3), d.rightMDLFspl(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'pArc')
    hv_child = nanmean(cat(2, d.leftpArc(group ~= 3), d.rightpArc(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftpArc(group == 3), d.rightpArc(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'tpc')
    hv_child = nanmean(cat(2, d.leftTPC(group ~= 3), d.rightTPC(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftTPC(group == 3), d.rightTPC(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'mdlfang')
    hv_child = nanmean(cat(2, d.leftMDLFang(group ~= 3), d.rightMDLFang(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftMDLFang(group == 3), d.rightMDLFang(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
elseif strcmp(tract, 'mdlfspl')
    hv_child = nanmean(cat(2, d.leftMDLFspl(group ~= 3), d.rightMDLFspl(group ~= 3)), 2);
    hv_adult = nanmean(cat(2, d.leftMDLFspl(group == 3), d.rightMDLFspl(group == 3)), 2);
    hv_n_child = length(find(~isnan(hv_child))); hv_n_adult = length(find(~isnan(hv_adult)));
end

% % Plot means.
% figure(1)
% hold on;
%
% child_mean = [nanmean(hv_child)];
% adult_mean = [nanmean(hv_adult)];
%
% child_sd = [nanstd(hv_child)];
% adult_sd = [nanstd(hv_adult)];
%
% xval = linspace(1, length(child_mean), length(child_mean));
%
% scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
% scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)
%
% errorbar(xval, adult_mean, adult_sd, 'Color', color_adults, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
% errorbar(xval, child_mean, child_sd, 'Color', color_children, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
%
% % xaxis
% xax = get(gca, 'xaxis');
% xax.Limits = [xlim_lo xlim_hi];
% xax.TickValues = xtickvalues;
% xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
% xlabels = {tract};
% xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
% xax.TickLabels = xlabels;
% xax.FontName = fontname;
% xax.FontSize = fontsize;
% xax.FontAngle = fontangle;
%
% % yaxis
% yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
% yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
% yax.FontName = fontname;
% yax.FontSize = fontsize;
%
% % general
% a = gca;
% %     a.TitleFontWeight = 'normal';
% box off
%
% legend({'Children', 'Adults'}, 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', fontsize);
% legend('boxoff');
%
% a.YLabel.String = 'Fractional Anisotropy (FA)';
% a.YLabel.FontSize = fontsize;
% pbaspect([1 1 1])
%
% %     pos=get(gca,'Position');
% %     pos1=pos-[0 .02 0 0];
% %     set(gca,'Position', pos1);
%
% % % Write.
% % if strcmp(save_figures, 'yes')
% %
% %     print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_' wm_measure '_' tract]), '-dpng')
% %     print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_' wm_measure '_' tract]), '-depsc')
% %
% % end
% %
% % hold off;

%% Plot same as before, but including individual data.

figure(2)
hold on;

coloralpha = .1;

child_mean = [nanmean(hv_child)];
adult_mean = [nanmean(hv_adult)];

child_sd = [nanstd(hv_child)];
adult_sd = [nanstd(hv_adult)];

xval = linspace(1, length(child_mean), length(child_mean));

% Means (do this first for legend and then second to keep it on top layer).
scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

% Individual data points for children.
scatter(repmat(1, size(hv_child)), hv_child, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Individual data points for adults.
scatter(repmat(1, size(hv_adult)), hv_adult, 'Marker', 'o', 'SizeData', markersize/2, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults, 'MarkerFaceAlpha', coloralpha, 'MarkerEdgeAlpha', coloralpha)

% Plot means again so that it is on top layer.
scatter(xval, child_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_children, 'MarkerEdgeColor', color_children)
scatter(xval, adult_mean, 'Marker', 'o', 'SizeData', 2*markersize, 'MarkerFaceColor', color_adults, 'MarkerEdgeColor', color_adults)

errorbar(xval, adult_mean, adult_sd, 'Color', color_adults, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
errorbar(xval, child_mean, child_sd, 'Color', color_children, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [xlim_lo xlim_hi];
xax.TickValues = xtickvalues;
xax.TickDirection = 'out';
xax.TickLength = [0 0];
xlabels = {tract};
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
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

%     pos=get(gca,'Position');
%     pos1=pos-[0 .02 0 0];
%     set(gca,'Position', pos1);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots-singleshell', ['plot_fig3_ind_' wm_measure '_' tract]), '-dpng')
    print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fig3_ind_' wm_measure '_' tract]), '-depsc')
    
end

hold off;

