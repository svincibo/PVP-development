% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/N/dc2/projects/lifebid/development/LWX_developmentOfVerticalWM/';

alphastat = 0.01;

hemisphere = 'both';% left, right, both
% if strcmp(hemisphere, 'right')
%     opt = 10;
% elseif strcmp(hemisphere, 'left')
%     opt = 0;
% else
%     opt = 20;
% end
wm_measure_here = {'fa', 'ad', 'md', 'rd', 'od', 'icvf', 'isovf'}; 
% lo = [.49 1.37 .81 .52 .21 .65 .11];
% mid = [.53 1.43 .87 .59 .23 .75 .13];
% hi = [.57 1.49 .93 .66 .25 .85 .15];
lo = [-.01 0 0 .65 -0.06 0 -.10];
mid = [0.00 .03 .06 .75 0 0.12 0];
hi = [0.10 .06 .12 .85 0.06 0.24 0.10];

verticalcolor = [0 128 128]/255; %teal
aslantcolor = [184 132 120]/255;
mdlfangcolor = [244 173 81]/255;
mdlfsplcolor = [43 102 112]/255;
vofcolor = [136 89 65]/255;
tpccolor = [116 207 227]/255;
pArccolor = [128 68 145]/255;

capsize = 0;
marker = 'o';
linewidth = 0.5;
linestyle = '-';
markersize = 100;
fontname = 'Arial';
fontsizex = 16; fontsizey = 16;
fontangle = 'italic';
fontcolor = [0 0 0];
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
ytickvalues = 1:13;
save_figures = 'yes';
ylim_lo = min(ytickvalues)-.5; ylim_hi = max(ytickvalues)+.5;
alphablend = .5;

%% WHITE MATTER MEASURES
for w = 1%1:length(wm_measure_here)
    
    % Read in data (from LWX_devOfVerticalWM_v3_loadData.m).
    data_tbl = readtable(fullfile(rootDir, 'supportFiles', ['LWX_devOfVerticalWM_forSPSS_' wm_measure_here{w} '_sosdenoised_new.csv']));
    
    % Convert into array and header for ease.
    data_all_in_header = data_tbl.Properties.VariableNames;
    data_all_in = table2array(data_tbl);
    
    % Get indices of subjects whose white matter values are all NaN.
    idx_notnan = ~all(isnan(data_all_in(:, 8:end)), 2);
    
    % Get easy index for age group.
    group = data_tbl.group_age(idx_notnan);
    
    % Remove subjects whose white matter values are all NaN from the table and from the array.
    data_tbl = data_tbl(idx_notnan, 8:end);
    data_all_in = data_all_in(idx_notnan, 8:end);
    data_all_in_header = data_tbl.Properties.VariableNames(:);
    
    % Display.
    disp([wm_measure_here{w}]);
    if strcmp(wm_measure_here{w}, 'fa')
        titlestring = 'Fractional Anisotropy (FA)';
    elseif strcmp(wm_measure_here{w}, 'od')
        titlestring = 'Orientation Dispersion (OD)';
    elseif strcmp(wm_measure_here{w}, 'icvf')
        titlestring = 'Neurite Density (ND)';
    elseif strcmp(wm_measure_here{w}, 'isovf')
        titlestring = 'Isotropic Volume Fraction (ISOVF)';
    end
    
    
    % Get the organizational indices.
    organizeT = [1 4 5 8 9 10 11 14 15 18 19 20];
    
    %% INDIVIDUAL TRACTS
    
    % Bootstrapping for each tract individually.
    for t = 1:length(organizeT)
        
        % Get data to plot.
        children = data_all_in(group ~= 3, organizeT(t));
        adults = data_all_in(group == 3, organizeT(t));
        
        % Bootstrap to get error bars because unequal sample sizes.
        for r = 1:10000
            
            this_c = randsample(children, size(children, 1), true);
            this_a = randsample(adults, size(adults, 1), true);
            
            this_mean_c = nanmean(this_c);
            this_mean_a = nanmean(this_a);
            
            diff_real(r) = this_mean_a - this_mean_c;
            
            clear this_c this_a this_mean_c this_mean_a
            
        end
        
        % creat the null distribution
        null_dis = cat(1, children, adults);
        
        for r = 1:10000
            
            this_c = randsample(null_dis, size(children, 1), true);
            this_a = randsample(null_dis, size(adults, 1), true);
            
            this_mean_c = nanmean(this_c);
            this_mean_a = nanmean(this_a);
            
            diff_null(r) = this_mean_a - this_mean_c;
            
            clear this_c this_a this_mean_c this_mean_a
            
        end
        
        ci = prctile(diff_real, [100*alphastat/2, 100*(1-alphastat/2)]);

        % Things to keep for plotting.
        tractname{t} = data_all_in_header{organizeT(t)};
        diff_real_ind(t) = nanmean(diff_real);
        ci_real_ind(t, :) = ci;
        
        tci = prctile(diff_null, 100*(1-alphastat));
        tci_real_ind(t) = tci;
        
        clear ci
        
    end
    
    
    %% VERTICAL GROUP
    
    % Get the organizational indices.
    dorsalT = [4 5 8 9 10 14 15 18 19 20];
    
    % Get data to plot.
    children = data_all_in(group ~= 3, dorsalT);
    adults = data_all_in(group == 3, dorsalT);
    
    % Bootstrap to get error bars because unequal sample sizes.
    for r = 1:10000
        
        for tract = 1:10
            
            this_c(:, tract) = randsample(children(:, tract), size(children, 1), true);
            this_a(:, tract) = randsample(adults(:, tract), size(adults, 1), true);
            
        end
        
        this_mean_c = nanmean(this_c(:));
        this_mean_a = nanmean(this_a(:));
        
        diff_real(r) = this_mean_a - this_mean_c;
        
        clear this_c this_a this_mean_c this_mean_a
        
    end
    
    % creat the null distribution
    null_dis = cat(1, children, adults);
    
    for r = 1:10000
        
        for tract = 1:10
            
            this_c(:, tract) = randsample(null_dis(:, tract), size(children, 1), true);
            this_a(:, tract) = randsample(null_dis(:, tract), size(adults, 1), true);
            
        end
        
        this_mean_c = nanmean(this_c(:));
        this_mean_a = nanmean(this_a(:));
        
        diff_null(r) = this_mean_a - this_mean_c;
        
        clear this_c this_a this_mean_c this_mean_a
        
    end
    
    ci = prctile(diff_real, [100*alphastat/2, 100*(1-alphastat/2)]);

    % Things to keep for plotting.
   d_diff = diff_real;
   diff_real_dorsal= nanmean(diff_real);
    ci_real_dorsal = ci;
    
    tci = prctile(diff_null, 100*(1-alphastat));
    tci_real_dorsal = tci;
    
end

% Plot
figure;
hold on;

% dorsal
b = barh([0 0 0 0 0 0 0 0 0 0 0 0 diff_real_dorsal]);
b.BarWidth = .8;
b.FaceColor = verticalcolor; b.FaceAlpha = .8;
b.EdgeColor = verticalcolor; b.EdgeAlpha = .8;

b = barh([0 0 0 0 0 0 0 0 0 0 0 diff_real_ind(1) 0]);
b.BarWidth = .3;
b.FaceColor = mdlfangcolor; b.FaceAlpha = alphablend;
b.EdgeColor = mdlfangcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 0 0 0 0 0 diff_real_ind(2) 0 0]);
b.BarWidth = .3;
b.FaceColor = mdlfsplcolor; b.FaceAlpha = alphablend;
b.EdgeColor = mdlfsplcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 0 0 0 0 diff_real_ind(3) 0 0 0]);
b.BarWidth = .3;
b.FaceColor = tpccolor; b.FaceAlpha = alphablend;
b.EdgeColor = tpccolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 0 0 0 diff_real_ind(4) 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = vofcolor; b.FaceAlpha = alphablend;
b.EdgeColor = vofcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 0 0 diff_real_ind(5) 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = pArccolor; b.FaceAlpha = alphablend;
b.EdgeColor = pArccolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 0 diff_real_ind(7) 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = mdlfangcolor; b.FaceAlpha = alphablend;
b.EdgeColor = mdlfangcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 diff_real_ind(8) 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = mdlfsplcolor; b.FaceAlpha = alphablend;
b.EdgeColor = mdlfsplcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 diff_real_ind(9) 0 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = tpccolor; b.FaceAlpha = alphablend;
b.EdgeColor = tpccolor; b.EdgeAlpha = alphablend;

b = barh([0 0 0 diff_real_ind(10) 0 0 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = vofcolor; b.FaceAlpha = alphablend;
b.EdgeColor = vofcolor; b.EdgeAlpha = alphablend;

b = barh([0 0 diff_real_ind(11) 0 0 0 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = pArccolor; b.FaceAlpha = alphablend;
b.EdgeColor = pArccolor; b.EdgeAlpha = alphablend;

b = barh([0 diff_real_ind(12) 0 0 0 0 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = pArccolor; b.FaceAlpha = 0.1;
b.EdgeColor = pArccolor; b.EdgeAlpha = 0.1;

b = barh([diff_real_ind(6) 0 0 0 0 0 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = pArccolor; b.FaceAlpha = 0.1;
b.EdgeColor = pArccolor; b.EdgeAlpha = 0.1;

% Confidence intervals
plot([ci_real_dorsal(1) ci_real_dorsal(2)], [13 13], 'Color', verticalcolor)
plot([ci_real_ind(1, 1) ci_real_ind(1, 2)], [12 12], 'Color', mdlfangcolor)
plot([ci_real_ind(2, 1) ci_real_ind(2, 2)], [11 11], 'Color', mdlfsplcolor)
plot([ci_real_ind(3, 1) ci_real_ind(3, 2)], [10 10], 'Color', tpccolor)
plot([ci_real_ind(4, 1) ci_real_ind(4, 2)], [9 9], 'Color', vofcolor)
plot([ci_real_ind(5, 1) ci_real_ind(5, 2)], [8 8], 'Color', pArccolor)
plot([ci_real_ind(7, 1) ci_real_ind(7, 2)], [7 7], 'Color', mdlfangcolor)
plot([ci_real_ind(8, 1) ci_real_ind(8, 2)], [6 6], 'Color', mdlfsplcolor)
plot([ci_real_ind(9, 1) ci_real_ind(9, 2)], [5 5], 'Color', tpccolor)
plot([ci_real_ind(10, 1) ci_real_ind(10, 2)], [4 4], 'Color', vofcolor)
plot([ci_real_ind(11, 1) ci_real_ind(11, 2)], [3 3], 'Color', pArccolor)
plot([ci_real_ind(12, 1) ci_real_ind(12, 2)], [2 2], 'Color', aslantcolor)
plot([ci_real_ind(6, 1) ci_real_ind(6, 2)], [1 1], 'Color', aslantcolor)

% % Plot means of null distributions
% cicolor = [224 128 128]/255;
% plot([tci_real_dorsal tci_real_dorsal], [9.5 10.5], 'Color', cicolor)
% plot([tci_real_ind(8) tci_real_ind(8)], [8.75 9.25], 'Color', cicolor)
% plot([tci_real_ind(7) tci_real_ind(7)], [7.75 8.25], 'Color', cicolor)
% plot([tci_real_ind(6) tci_real_ind(6)], [6.75 7.25], 'Color', cicolor)
% plot([tci_real_ind(5) tci_real_ind(5)], [5.75 6.25], 'Color', cicolor)
% plot([tci_real_ventral tci_real_ventral], [4.5 5.5], 'Color', cicolor)
% plot([tci_real_ind(4) tci_real_ind(4)], [3.75 4.25], 'Color', cicolor)
% plot([tci_real_ind(3) tci_real_ind(3)], [2.75 3.25], 'Color', cicolor)
% plot([tci_real_ind(2) tci_real_ind(2)], [1.75 2.25], 'Color', cicolor)
% plot([tci_real_ind(1) tci_real_ind(1)], [0.75 1.25], 'Color', cicolor)

% Cover the yaxis - aesthetic.
plot([0 0], [0 11], 'k')

% Separate posterior from anterior.
plot([0 1], [2.5 2.5], 'k')

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = ytickvalues;
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = [tractname([7 1 12 11 10 9 8 6 5 4 3 2]), 'Vertical (avg)'];
yax.FontName = fontname;
yax.FontSize = fontsizey;
yax.FontSmoothing = fontsmoothing;
yax.Color = fontcolor;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [lo(w) hi(w)];
xax.TickValues = [lo(w) mid(w) hi(w)];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.TickLabels = {num2str(lo(w), '%1.2f'), '', num2str(hi(w), '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsizex;
xax.FontAngle = fontangle;
xax.FontSmoothing = fontsmoothing;

a = gca;
a.TitleFontWeight = 'normal';
box off
xlabel({['Difference ' titlestring ':']; 'Adults - Children'}, 'FontName', fontname, 'FontSize', fontsizex, 'FontAngle', fontangle, 'Color', fontcolor, 'FontSmoothing', fontsmoothing);
a.XLabel.FontSize = fontsizex;
pbaspect([1 1 1]);

% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots', ['plot_gpmeans_' wm_measure_here{w} '_sosdeniosed_' hemisphere '_diff_fig1_vertical2']), '-dpng', '-r600')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_gpmeans_' wm_measure_here{w} '_sosdenoised_' hemisphere '_diff_fig1_vertical2']), '-depsc', '-r600')
    
end

hold off;


