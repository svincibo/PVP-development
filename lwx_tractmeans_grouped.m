clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

wm_measure_here = {'fa'};

h = {'both'}; %left, right , both

save_figures = 'yes';

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

%% WHITE MATTER MEASURES
f = 0;
for w = 1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    data_tbl = readtable(fullfile(rootDir, 'supportFiles', ['LWX_data_forSPSS_' wm_measure_here{w} '_singleshell.csv']));
    
    % Get indices of subjects whose white matter values are all NaN.
    idx_notnan = ~any(isnan(table2array(data_tbl(:, 8:end))), 2);
    
    % Remove subjects whose white matter values are all NaN from the table and from the array.
    data_tbl = data_tbl(idx_notnan, :);
    
    % Get index for outliers to be removed.
    idx_keep = find(~ismember(data_tbl.subID, outlier));
    
    % Remove outliers.
    data_tbl = data_tbl(idx_keep, :);
    
    % Display.
    disp([wm_measure_here{w}]);
    
    % Get easy index for age group.
    group = data_tbl.group_age;
    
    % Set up plot and measure-specific details.
    capsize = 0;
    marker = 'o';
    linewidth = 1.5;
    linestyle = 'none';
    markersize = 10;
    xtickvalues = [1 2];
    xlim_lo = 0.5; xlim_hi = 2.5;
    fontname = 'Arial';
    fontsizex = 16; fontsizey = 16;
    fontsize = 16;
    fontangle = 'italic';
    yticklength = 0;
    xticklength = 0.05;
    if strcmp(wm_measure_here{w}, 'odi')
        ylimlo = 0.10; ylimhi = 0.25;
        yname = 'Orientation Dispersion Index (ODI)';
    elseif strcmp(wm_measure_here{w}, 'ndi')
        ylimlo = 0.45; ylimhi = 0.70;
        yname = 'Neurite Density Index (NDI)';
    elseif strcmp(wm_measure_here{w}, 'isovf')
        ylimlo = 0; ylimhi = 1;
        yname = 'Isotropic Volume Fraction (ISOVF)';
    elseif strcmp(wm_measure_here{w}, 'fa')
        ylimlo = 0.40; ylimhi = 0.60;
        yname = 'Fractional Anisotropy (FA)';
    end
    
    % Plot means.
    f = f+1;
    figure(f)
    
    % Horizontal, Ventral
    if strcmp(h, 'both')
        toplot_1 = cat(1, data_tbl.leftILF(group ~= 3), data_tbl.rightILF(group ~= 3), ...
            data_tbl.leftIFOF(group ~= 3), data_tbl.rightIFOF(group ~= 3));
        toplot_3 = cat(1, data_tbl.leftILF(group == 3), data_tbl.rightILF(group == 3), ...
            data_tbl.leftIFOF(group == 3), data_tbl.rightIFOF(group == 3));
        n1 = length(toplot_1)/2; n3 = length(toplot_3)/2;
    elseif strcmp(h, 'left')
        toplot_1 = data_tbl.leftILF(group == 1); toplot_2 = data_tbl.leftILF(group == 2); toplot_3 = data_tbl.leftILF(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    elseif strcmp(h, 'right')
        toplot_1 = data_tbl.rightILF(group == 1); toplot_2 = data_tbl.rightILF(group == 2); toplot_3 = data_tbl.rightILF(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    end
    markeredgecolor = [0 0.19 0.56]; markerfacecolor = [0 0.19 0.56]; % dark blue
    bar(xtickvalues-.3, [mean(toplot_1) mean(toplot_3)], 'BarWidth', .2, 'FaceColor', markerfacecolor, 'EdgeColor', markeredgecolor)
    hold on;
    errorbars = errorbar(xtickvalues-.3, [mean(toplot_1) mean(toplot_3)], ...
        1.96*[std(toplot_1, 0, 1)/sqrt(n1), std(toplot_3, 0, 1)/sqrt(n3)]);
    set(errorbars, 'Color', markeredgecolor, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    d1 = mean(toplot_3)-mean(toplot_1);
    
    % Horizontal, Dorsal
    if strcmp(h, 'both')
        toplot_1 = cat(1, data_tbl.leftSLF1And2(group ~= 3), data_tbl.rightSLF1And2(group ~= 3), ...
            data_tbl.leftSLF3(group ~= 3), data_tbl.rightSLF3(group ~= 3));
        toplot_3 = cat(1, data_tbl.leftSLF1And2(group == 3), data_tbl.rightSLF1And2(group == 3), ...
            data_tbl.leftSLF3(group == 3), data_tbl.rightSLF3(group == 3));
        n1 = length(toplot_1)/2; n3 = length(toplot_3)/2;
    elseif strcmp(h, 'left')
        toplot_1 = data_tbl.leftSLF1And2(group == 1); toplot_2 = data_tbl.leftSLF1And2(group == 2); toplot_3 = data_tbl.leftSLF1And2(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    elseif strcmp(h, 'right')
        toplot_1 = data_tbl.rightSLF1And2(group == 1); toplot_2 = data_tbl.rightSLF1And2(group == 2); toplot_3 = data_tbl.rightSLF1And2(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    end
    markeredgecolor = [1.00, 0.49, 0.00]; markerfacecolor = [1.00, 0.49, 0.00]; % dark amber
    bar(xtickvalues-.1, [mean(toplot_1) mean(toplot_3)], 'BarWidth', .2, 'FaceColor', markerfacecolor, 'EdgeColor', markeredgecolor)
    hold on;
    errorbars = errorbar(xtickvalues-.1, [mean(toplot_1) mean(toplot_3)], ...
        1.96*[std(toplot_1, 0, 1)/sqrt(n1), std(toplot_3, 0, 1)/sqrt(n3)]);
    set(errorbars, 'Color', markeredgecolor, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    d2 = mean(toplot_3)-mean(toplot_1);
    
    % Vertical, Anterior
    if strcmp(h, 'both')
        toplot_1 = cat(1, data_tbl.leftAslant(group ~= 3), data_tbl.rightAslant(group ~= 3));
        toplot_3 = cat(1, data_tbl.leftAslant(group == 3), data_tbl.rightAslant(group == 3));
        n1 = length(toplot_1)/2; n3 = length(toplot_3)/2;
    elseif strcmp(h, 'left')
        toplot_1 = data_tbl.leftAslant(group == 1); toplot_2 = data_tbl.leftAslant(group == 2); toplot_3 = data_tbl.leftAslant(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    elseif strcmp(h, 'right')
        toplot_1 = data_tbl.rightAslant(group == 1); toplot_2 = data_tbl.rightAslant(group == 2); toplot_3 = data_tbl.rightAslant(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    end
    markeredgecolor = [0.23 0.48 0.34]; markerfacecolor = [0.23 0.48 0.34]; % sea green
    bar(xtickvalues+.1, [mean(toplot_1) mean(toplot_3)], 'BarWidth', .2, 'FaceColor', markerfacecolor, 'EdgeColor', markeredgecolor)
    errorbars = errorbar(xtickvalues+.1, [mean(toplot_1) mean(toplot_3)], ...
        1.96*[std(toplot_1, 0, 1)/sqrt(n1), std(toplot_3, 0, 1)/sqrt(n3)]);
    set(errorbars, 'Color', markeredgecolor, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    d3 = mean(toplot_3)-mean(toplot_1);
    
    % Vertical, Posterior
    if strcmp(h, 'both')
        toplot_1 = cat(1, data_tbl.leftMDLFang(group ~= 3), data_tbl.rightMDLFang(group ~= 3), ...
            data_tbl.leftMDLFspl(group ~= 3), data_tbl.rightMDLFspl(group ~= 3), ...
            data_tbl.leftTPC(group ~= 3), data_tbl.rightTPC(group ~= 3), ...
            data_tbl.leftpArc(group ~= 3), data_tbl.rightpArc(group ~= 3)); %, ...
%             data_tbl.leftVOF(group ~= 3), data_tbl.rightVOF(group ~= 3));
        toplot_3 = cat(1, data_tbl.leftMDLFang(group == 3), data_tbl.rightMDLFang(group == 3), ...
            data_tbl.leftMDLFspl(group == 3), data_tbl.rightMDLFspl(group == 3), ...
            data_tbl.leftTPC(group == 3), data_tbl.rightTPC(group == 3), ...
            data_tbl.leftpArc(group == 3), data_tbl.rightpArc(group == 3));%, ...
%             data_tbl.leftVOF(group == 3), data_tbl.rightVOF(group == 3));
        n1 = length(toplot_1)/2; n3 = length(toplot_3)/2;
    elseif strcmp(h, 'left')
        toplot_1 = data_tbl.leftMDLFang(group == 1); toplot_2 = data_tbl.leftMDLFang(group == 2); toplot_3 = data_tbl.leftMDLFang(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    elseif strcmp(h, 'right')
        toplot_1 = data_tbl.rightMDLFang(group == 1); toplot_2 = data_tbl.rightMDLFang(group == 2); toplot_3 = data_tbl.rightMDLFang(group == 3);
        n1 = length(toplot_1); n2 = length(toplot_2); n3 = length(toplot_3);
    end
    markeredgecolor = [0.67 0.15 0.31]; markerfacecolor = [0.67 0.15 0.31]; % dark pink
    bar(xtickvalues+.3, [mean(toplot_1) mean(toplot_3)], 'BarWidth', .2, 'FaceColor', markerfacecolor, 'EdgeColor', markeredgecolor)
    errorbars = errorbar(xtickvalues+.3, [mean(toplot_1) mean(toplot_3)], ...
        1.96*[std(toplot_1, 0, 1)/sqrt(n1), std(toplot_3, 0, 1)/sqrt(n3)]);
    set(errorbars, 'Color', markeredgecolor, 'LineWidth', linewidth, 'LineStyle', linestyle, 'CapSize', capsize);
    d4 = mean(toplot_3)-mean(toplot_1);
    
    % xaxis
    xax = get(gca, 'xaxis');
    xax.Limits = [xlim_lo xlim_hi];
    xax.TickValues = xtickvalues;
    xax.TickDirection = 'out';
    xax.TickLength = [xticklength xticklength];
    xlabels = {'Children', 'Adults'};
    xlabels = cellfun(@(x) strrep(x, ',', '\newline'), xlabels, 'UniformOutput', false);
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
    
    legend({'', 'Horizontal, Ventral', '', 'Horizontal, Dorsal', '', ...
        'Vertical, Anterior (FAT)', '', 'Vertical, Posterior'}, 'Location', 'northwest', 'Orientation', 'vertical');
    legend('boxoff');
    
    a.YLabel.String = yname;
    a.YLabel.FontSize = fontsize;
    pbaspect([1 1 1])
    
    %     pos=get(gca,'Position');
    %     pos1=pos-[0 .02 0 0];
    %     set(gca,'Position', pos1);
    
    % Write.
    if strcmp(save_figures, 'yes')
        
        print(fullfile(rootDir, 'plots-singleshell', ['plot_singleshell_anova_' wm_measure_here{w} '_yoa_' h{1}]), '-dpng')
        print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_singleshell_anova_' wm_measure_here{w} '_yoa_' h{1}]), '-depsc')
        
    end
    
    hold off;
    
%     f=f+1;
%     figure(f)
%     bar([1 2 3 4], [d1 d2 d3 d4])
%     ylim([-.01 .1]);
    
end


