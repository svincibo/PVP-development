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

dorsalcolor = [236 176 32]/255; %burnt yellow
slf12color = [204 148 29]/255;
slf3color = [255 229 173]/255;

ventralcolor= [14 114 184]/255; % blue
ifofcolor = [13 73 143]/255;
ilfcolor = [94 188 255]/255;

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
ytickvalues = 1:10;
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
    organizeT = [12 13 2 3 17 16 7 6];
    
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
    
    
    %% DORSAL GROUP
    
    % Get the organizational indices.
    dorsalT = [17 16 7 6];
    
    % Get data to plot.
    children = data_all_in(group ~= 3, dorsalT);
    adults = data_all_in(group == 3, dorsalT);
    
    % Bootstrap to get error bars because unequal sample sizes.
    for r = 1:10000
        
        for tract = 1:4
            
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
        
        for tract = 1:4
            
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
    
    %% VENTRAL GROUP
    
    % Get the organizational indices.
    ventralT = [12 13 2 3];
    
    % Get data to plot.
    children = data_all_in(group ~= 3, ventralT);
    adults = data_all_in(group == 3, ventralT);
    
    % Bootstrap to get error bars because unequal sample sizes.
    for r = 1:10000
        
        for tract = 1:4
            
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
        
        for tract = 1:4
            
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
    v_diff = diff_real;
    diff_real_ventral= nanmean(diff_real);
    ci_real_ventral = ci;
    
    tci = prctile(diff_null, 100*(1-alphastat));
    tci_real_ventral = tci;
    
end

% Plot
figure;
hold on;

% dorsal
b = barh([0 0 0 0 0 0 0 0 0 diff_real_dorsal]);
b.BarWidth = .8;
b.FaceColor = dorsalcolor; b.FaceAlpha = .8;
b.EdgeColor = dorsalcolor; b.EdgeAlpha = .8;

b = barh([0 0 0 0 0 0 diff_real_ind(6) 0 diff_real_ind(8) 0]);
b.BarWidth = .3;
b.FaceColor = slf12color; b.FaceAlpha = alphablend;
b.EdgeColor = slf12color; b.EdgeAlpha = alphablend;

b = barh([0 0 0 0 0 diff_real_ind(5) 0 diff_real_ind(7) 0 0]);
b.BarWidth = .3;
b.FaceColor = slf3color; b.FaceAlpha = alphablend;
b.EdgeColor = slf3color; b.EdgeAlpha = alphablend;

% ventral
b = barh([0 0 0 0 diff_real_ventral 0 0 0 0 0]);
b.BarWidth = .8;
b.FaceColor = ventralcolor; b.FaceAlpha = .8;
b.EdgeColor = ventralcolor; b.EdgeAlpha = .8;

b = barh([0 diff_real_ind(2) 0 diff_real_ind(4) 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = ilfcolor; b.FaceAlpha = alphablend;
b.EdgeColor = ilfcolor; b.EdgeAlpha = alphablend;

b = barh([diff_real_ind(1) 0 diff_real_ind(3) 0 0 0 0 0 0 0]);
b.BarWidth = .3;
b.FaceColor = ifofcolor; b.FaceAlpha = alphablend;
b.EdgeColor = ifofcolor; b.EdgeAlpha = alphablend;

% Confidence intervals
plot([ci_real_dorsal(1) ci_real_dorsal(2)], [10 10], 'Color', dorsalcolor)
plot([ci_real_ind(8, 1) ci_real_ind(8, 2)], [9 9], 'Color', slf12color)
plot([ci_real_ind(7, 1) ci_real_ind(7, 2)], [8 8], 'Color', slf3color)
plot([ci_real_ind(6, 1) ci_real_ind(6, 2)], [7 7], 'Color', slf12color)
plot([ci_real_ind(5, 1) ci_real_ind(5, 2)], [6 6], 'Color', slf3color)
plot([ci_real_ventral(1) ci_real_ventral(2)], [5 5], 'Color', ventralcolor)
plot([ci_real_ind(4, 1) ci_real_ind(4, 2)], [4 4], 'Color', ilfcolor)
plot([ci_real_ind(3, 1) ci_real_ind(3, 2)], [3 3], 'Color', ifofcolor)
plot([ci_real_ind(2, 1) ci_real_ind(2, 2)], [2 2], 'Color', ilfcolor)
plot([ci_real_ind(1, 1) ci_real_ind(1, 2)], [1 1], 'Color', ifofcolor)

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

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = ytickvalues;
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = [tractname(1:4), 'Ventral (avg)', tractname(5:8), 'Dorsal (avg)'];
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
    
    print(fullfile(rootDir, 'plots', ['plot_gpmeans_' wm_measure_here{w} '_sosdenoised_' hemisphere '_diff_fig1']), '-dpng', '-r600')
    print(fullfile(rootDir, 'plots', 'eps', ['plot_gpmeans_' wm_measure_here{w} '_sosdenoised_' hemisphere '_diff_fig1']), '-depsc', '-r600')
    
end

hold off;

%% DORSAL - VENTRAL DIFFERENCE TEST

% test for difference betwen s_dorsal and s_ventral (bootstrapped samples of adult-child differences in dorsal and ventral fa).

organizeT = [12 13 2 3 17 16 7 6];

% Get data to plot.
children = data_all_in(group ~= 3, organizeT);
adults = data_all_in(group == 3, organizeT);

% Get real distribution, for previous bootstraps.
mu_diff_real = nanmean(d_diff - v_diff);
sigma_diff_real = nanstd(d_diff - v_diff);

% Get null distribution.
null_dis = cat(1, children, adults);

for r = 1:10000
    
    % Dorsal
    for tract = 1:4
        
        this_cd(:, tract) = randsample(null_dis(:, tract), size(children, 1), true);
        this_ad(:, tract) = randsample(null_dis(:, tract), size(adults, 1), true);
        
    end
    
    % Ventral
    for tract = 1:4
        
        this_cv(:, tract) = randsample(null_dis(:, tract), size(children, 1), true);
        this_av(:, tract) = randsample(null_dis(:, tract), size(adults, 1), true);
        
    end
    
    this_mean_d = nanmean(this_ad(:)) - nanmean(this_cd(:));
    this_mean_v = nanmean(this_av(:)) - nanmean(this_cv(:));
    
    diff_null(r) = this_mean_d - this_mean_v;
    
    clear this_cd this_ad this_cv this_av this_mean_d this_mean_v
    
end

mu_diff_null = nanmean(diff_null);
sigma_diff_null = nanstd(diff_null);

z = (mu_diff_real - mu_diff_null)./sigma_diff_null;
p = 1-normcdf(abs(z), 0, 1);
