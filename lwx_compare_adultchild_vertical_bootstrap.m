clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

alphastat = 0.01;
multcompcorrection = 'no';

hemisphere = 'both';% left, right, both

%% READ IN DATA AND ORGANIZE.

% Read in 'final' data. This data should have all subjects removed, cleaned, etc.
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));

% SELECT the measurements of the tracts that I care about.
if strcmp(hemisphere, 'left')
    
    % Get indices of columns that correspond to tracts of interest.
    toi = find(strcmp(d.Properties.VariableNames, 'leftpArc') + strcmp(d.Properties.VariableNames, 'leftTPC') +  ...
        + strcmp(d.Properties.VariableNames, 'leftMDLFspl') + strcmp(d.Properties.VariableNames, 'leftMDLFang'));
    
elseif strcmp(hemisphere, 'right')
    
    % Get indices of columns that correspond to tracts of interest.
    toi = find(strcmp(d.Properties.VariableNames, 'rightpArc') + strcmp(d.Properties.VariableNames, 'rightTPC') +  ...
        + strcmp(d.Properties.VariableNames, 'rightMDLFspl') + strcmp(d.Properties.VariableNames, 'rightMDLFang'));
    
elseif strcmp(hemisphere, 'both')
    
    % Get indices of columns that correspond to tracts of interest.
    toi = find(strcmp(d.Properties.VariableNames, 'leftpArc') + strcmp(d.Properties.VariableNames, 'rightpArc')  ...
        + strcmp(d.Properties.VariableNames, 'leftTPC') + strcmp(d.Properties.VariableNames, 'rightTPC') ...
        + strcmp(d.Properties.VariableNames, 'leftMDLFspl') + strcmp(d.Properties.VariableNames, 'rightMDLFspl') ...
        + strcmp(d.Properties.VariableNames, 'leftMDLFang') + strcmp(d.Properties.VariableNames, 'rightMDLFang'));
    
end

%% INDIVIDUAL TRACTS

% Bootstrapping for each tract individually.
for t = 1:length(toi)
    
    % Get data to plot.
    children = table2array(d(d.group_age3 ~= 3, toi(t)));
    adults = table2array(d(d.group_age3 == 3, toi(t)));
    
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
    tractname{t} = d.Properties.VariableNames{toi(t)};
    diff_real_ind(t) = nanmean(diff_real);
    ci_real_ind(t, :) = ci;
    
    tci = prctile(diff_null, 100*(1-alphastat));
    tci_real_ind(t) = tci;
    disp([num2str(t) ' ' tractname{t}])
    
    clear ci
    
end


%% VERTICAL GROUP

% Get data to plot.
children = table2array(d(d.group_age3 ~= 3, toi));
adults = table2array(d(d.group_age3 == 3, toi));

% Bootstrap to get error bars because unequal sample sizes.
for r = 1:10000
    
    for tract = 1:length(toi)
        
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
    
    for tract = 1:length(toi)
        
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
diff_real_vertical= nanmean(diff_real);
ci_real_vertical = ci;

tci = prctile(diff_null, 100*(1-alphastat));
tci_real_vertical = tci;

%% PLOT

dorsalcolor = [0 128 128]/255; %teal
leftcolor = [105 105 105]/255; %dark gray
rightcolor = [169 169 169]/255; %light gray

lo = -.01;
mid = 0.05;
hi = 0.10;

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
save_figures = 'yes';
alphablend = .5;

figure;
hold on;

ytickvalues = 1:9;
ylim_lo = min(ytickvalues)-.5; ylim_hi = max(ytickvalues)+.5;

% dorsal
b = barh([0 0 0 0 0 0 0 0 diff_real_vertical]);
b.BarWidth = .8;
b.FaceColor = dorsalcolor; b.FaceAlpha = .8;
b.EdgeColor = dorsalcolor; b.EdgeAlpha = .8;

b = barh([ 0 0 0 0 diff_real_ind(2) diff_real_ind(2) diff_real_ind(3) diff_real_ind(4) 0]); %left
b.BarWidth = .3;
b.FaceColor = leftcolor; b.FaceAlpha = alphablend;
b.EdgeColor = leftcolor; b.EdgeAlpha = alphablend;

b = barh([diff_real_ind(5) diff_real_ind(6) diff_real_ind(7) diff_real_ind(8) 0 0 0 0 0]); %right
b.BarWidth = .3;
b.FaceColor = rightcolor; b.FaceAlpha = alphablend;
b.EdgeColor = rightcolor; b.EdgeAlpha = alphablend;

% Confidence intervals
plot([ci_real_vertical(1) ci_real_vertical(2)], [9 9], 'Color', dorsalcolor)
plot([ci_real_ind(4, 1) ci_real_ind(4, 2)], [8 8], 'Color', leftcolor)
plot([ci_real_ind(3, 1) ci_real_ind(3, 2)], [7 7], 'Color', leftcolor)
plot([ci_real_ind(2, 1) ci_real_ind(2, 2)], [6 6], 'Color', leftcolor)
plot([ci_real_ind(1, 1) ci_real_ind(1, 2)], [5 5], 'Color', leftcolor)
plot([ci_real_ind(8, 1) ci_real_ind(8, 2)], [4 4], 'Color', rightcolor)
plot([ci_real_ind(7, 1) ci_real_ind(7, 2)], [3 3], 'Color', rightcolor)
plot([ci_real_ind(6, 1) ci_real_ind(6, 2)], [2 2], 'Color', rightcolor)
plot([ci_real_ind(5, 1) ci_real_ind(5, 2)], [1 1], 'Color', rightcolor)

ticklabels = [tractname([5, 6, 7, 8, 1, 2, 3, 4]), 'Vertical (avg)'];

% Cover the yaxis - aesthetic.
plot([0 0], [0 max(ytickvalues)+1], 'k')

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = ytickvalues;
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = ticklabels;
yax.FontName = fontname;
yax.FontSize = fontsizey;
yax.FontSmoothing = fontsmoothing;
yax.Color = fontcolor;

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [lo hi];
xax.TickValues = [lo+.01 mid hi];
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.TickLabels = {num2str(0, '%1.2f'), '', num2str(hi, '%1.2f')};
xax.FontName = fontname;
xax.FontSize = fontsizex;
xax.FontAngle = fontangle;
xax.FontSmoothing = fontsmoothing;

a = gca;
titlestring = 'Fractional Anisotropy (FA)';
a.TitleFontWeight = 'normal';
box off
xlabel({['Difference ' titlestring ':']; 'Adults - Children'}, 'FontName', fontname, 'FontSize', fontsizex, 'FontAngle', fontangle, 'Color', fontcolor, 'FontSmoothing', fontsmoothing);
a.XLabel.FontSize = fontsizex;
pbaspect([1 1 1]);

% Write.
print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_adultchild_vertical_' hemisphere '_singleshell']), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_adultchild_vertical_' hemisphere '_singleshell']), '-depsc', '-r600')

hold off;

