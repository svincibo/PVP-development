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
    toi = find(strcmp(d.Properties.VariableNames, 'leftSLF1And2') + strcmp(d.Properties.VariableNames, 'leftSLF3')  ...
        + strcmp(d.Properties.VariableNames, 'leftILF') + strcmp(d.Properties.VariableNames, 'leftIFOF'));
    
    % Get indices of columns that correspond to dorsal tracts.
    dorsal = find(strcmp(d.Properties.VariableNames, 'leftSLF1And2') + strcmp(d.Properties.VariableNames, 'leftSLF3'));
    
    % Get indices of columns that correspond to ventral tracts.
    ventral = find(strcmp(d.Properties.VariableNames, 'leftILF') + strcmp(d.Properties.VariableNames, 'leftIFOF'));
    
elseif strcmp(hemisphere, 'right')
    
    % Get indices of columns that correspond to tracts of interest.
    toi = find(strcmp(d.Properties.VariableNames, 'rightSLF1And2') + strcmp(d.Properties.VariableNames, 'rightSLF3')  ...
        + strcmp(d.Properties.VariableNames, 'rightILF') + strcmp(d.Properties.VariableNames, 'rightIFOF'));
    
    % Get indices of columns that correspond to dorsal tracts.
    dorsal = find(strcmp(d.Properties.VariableNames, 'rightSLF1And2') + strcmp(d.Properties.VariableNames, 'rightSLF3') );
    
    % Get indices of columns that correspond to ventral tracts.
    ventral = find(strcmp(d.Properties.VariableNames, 'rightILF') + strcmp(d.Properties.VariableNames, 'rightIFOF'));
    
elseif strcmp(hemisphere, 'both')
    
    % Get indices of columns that correspond to tracts of interest.
    toi = find(strcmp(d.Properties.VariableNames, 'leftSLF1And2') + strcmp(d.Properties.VariableNames, 'rightSLF1And2') ...
        + strcmp(d.Properties.VariableNames, 'leftSLF3') + strcmp(d.Properties.VariableNames, 'rightSLF3') ...
        + strcmp(d.Properties.VariableNames, 'leftILF') + strcmp(d.Properties.VariableNames, 'rightILF') ...
        + strcmp(d.Properties.VariableNames, 'leftIFOF') + strcmp(d.Properties.VariableNames, 'rightIFOF'));
    
    % Get indices of columns that correspond to dorsal tracts.
    dorsal = find(strcmp(d.Properties.VariableNames, 'leftSLF1And2') + strcmp(d.Properties.VariableNames, 'rightSLF1And2') ...
        + strcmp(d.Properties.VariableNames, 'leftSLF3') + strcmp(d.Properties.VariableNames, 'rightSLF3'));
    
    % Get indices of columns that correspond to ventral tracts.
    ventral = find(strcmp(d.Properties.VariableNames, 'leftILF') + strcmp(d.Properties.VariableNames, 'rightILF') ...
        + strcmp(d.Properties.VariableNames, 'leftIFOF') + strcmp(d.Properties.VariableNames, 'rightIFOF'));
    
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


%% DORSAL GROUP

% Get data to plot.
children = table2array(d(d.group_age3 ~= 3, dorsal));
adults = table2array(d(d.group_age3 == 3, dorsal));

% Bootstrap to get error bars because unequal sample sizes.
for r = 1:10000
    
    for tract = 1:length(dorsal)
        
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
    
    for tract = 1:length(dorsal)
        
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

% Get data to plot.
children = table2array(d(d.group_age3 ~= 3, ventral));
adults = table2array(d(d.group_age3 == 3, ventral));

% Bootstrap to get error bars because unequal sample sizes.
for r = 1:10000
    
    for tract = 1:length(ventral)
        
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
    
    for tract = 1:length(ventral)
        
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

%% PLOT

dorsalcolor = [236 176 32]/255; %burnt yellow
slf12color = [204 148 29]/255;
slf3color = [255 229 173]/255;

ventralcolor= [14 114 184]/255; % blue
ifofcolor = [13 73 143]/255;
ilfcolor = [94 188 255]/255;

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

if strcmp(hemisphere, 'both')
    
    ytickvalues = 1:10;
    ylim_lo = min(ytickvalues)-.5; ylim_hi = max(ytickvalues)+.5;
    
    % dorsal
    b = barh([0 0 0 0 0 0 0 0 0 diff_real_dorsal]);
    b.BarWidth = .8;
    b.FaceColor = dorsalcolor; b.FaceAlpha = .8;
    b.EdgeColor = dorsalcolor; b.EdgeAlpha = .8;
    
    b = barh([0 0 0 0 0 0 diff_real_ind(7) 0 diff_real_ind(3) 0]); %SLF1And2
    b.BarWidth = .3;
    b.FaceColor = slf12color; b.FaceAlpha = alphablend;
    b.EdgeColor = slf12color; b.EdgeAlpha = alphablend;
    
    b = barh([0 0 0 0 0 diff_real_ind(8) 0 diff_real_ind(4) 0 0]); %SLF3
    b.BarWidth = .3;
    b.FaceColor = slf3color; b.FaceAlpha = alphablend;
    b.EdgeColor = slf3color; b.EdgeAlpha = alphablend;
    
    % ventral
    b = barh([0 0 0 0 diff_real_ventral 0 0 0 0 0]);
    b.BarWidth = .8;
    b.FaceColor = ventralcolor; b.FaceAlpha = .8;
    b.EdgeColor = ventralcolor; b.EdgeAlpha = .8;
    
    b = barh([0 diff_real_ind(6) 0 diff_real_ind(2) 0 0 0 0 0 0]); %ILF
    b.BarWidth = .3;
    b.FaceColor = ilfcolor; b.FaceAlpha = alphablend;
    b.EdgeColor = ilfcolor; b.EdgeAlpha = alphablend;
    
    b = barh([diff_real_ind(5) 0 diff_real_ind(1) 0 0 0 0 0 0 0]); %IFOF
    b.BarWidth = .3;
    b.FaceColor = ifofcolor; b.FaceAlpha = alphablend;
    b.EdgeColor = ifofcolor; b.EdgeAlpha = alphablend;
    
    % Confidence intervals
    plot([ci_real_dorsal(1) ci_real_dorsal(2)], [10 10], 'Color', dorsalcolor)
    plot([ci_real_ind(3, 1) ci_real_ind(3, 2)], [9 9], 'Color', slf12color)
    plot([ci_real_ind(4, 1) ci_real_ind(4, 2)], [8 8], 'Color', slf3color)
    plot([ci_real_ind(7, 1) ci_real_ind(7, 2)], [7 7], 'Color', slf12color)
    plot([ci_real_ind(8, 1) ci_real_ind(8, 2)], [6 6], 'Color', slf3color)
    plot([ci_real_ventral(1) ci_real_ventral(2)], [5 5], 'Color', ventralcolor)
    plot([ci_real_ind(2, 1) ci_real_ind(2, 2)], [4 4], 'Color', ilfcolor)
    plot([ci_real_ind(1, 1) ci_real_ind(1, 2)], [3 3], 'Color', ifofcolor)
    plot([ci_real_ind(6, 1) ci_real_ind(6, 2)], [2 2], 'Color', ilfcolor)
    plot([ci_real_ind(5, 1) ci_real_ind(5, 2)], [1 1], 'Color', ifofcolor)

    ticklabels = [tractname([5, 6, 1, 2]), 'Ventral (avg)', tractname([8, 7, 4, 3]), 'Dorsal (avg)'];
    
elseif strcmp(hemisphere, 'left') || strcmp(hemisphere, 'right')
    
    ytickvalues = 1:6;
    ylim_lo = min(ytickvalues)-.5; ylim_hi = max(ytickvalues)+.5;
    
    % dorsal
    b = barh([0 0 0 0 0 diff_real_dorsal]);
    b.BarWidth = .8;
    b.FaceColor = dorsalcolor; b.FaceAlpha = .8;
    b.EdgeColor = dorsalcolor; b.EdgeAlpha = .8;
    
    b = barh([0 0 0 diff_real_ind(4) diff_real_ind(3) 0]); %SLF1And2, SLF3
    b.BarWidth = .3;
    b.FaceColor = slf12color; b.FaceAlpha = alphablend;
    b.EdgeColor = slf12color; b.EdgeAlpha = alphablend;
    
    % ventral
    b = barh([0 0 diff_real_ventral 0 0 0]);
    b.BarWidth = .8;
    b.FaceColor = ventralcolor; b.FaceAlpha = .8;
    b.EdgeColor = ventralcolor; b.EdgeAlpha = .8;
    
    b = barh([diff_real_ind(1) diff_real_ind(2) 0 0 0]); %ILF, IFOF
    b.BarWidth = .3;
    b.FaceColor = ilfcolor; b.FaceAlpha = alphablend;
    b.EdgeColor = ilfcolor; b.EdgeAlpha = alphablend;
    
    % Confidence intervals
    plot([ci_real_dorsal(1) ci_real_dorsal(2)], [6 6], 'Color', dorsalcolor)
    plot([ci_real_ind(3, 1) ci_real_ind(3, 2)], [5 5], 'Color', slf12color)
    plot([ci_real_ind(4, 1) ci_real_ind(4, 2)], [4 4], 'Color', slf12color)
    plot([ci_real_ventral(1) ci_real_ventral(2)], [3 3], 'Color', ventralcolor)
    plot([ci_real_ind(2, 1) ci_real_ind(2, 2)], [2 2], 'Color', ilfcolor)
    plot([ci_real_ind(1, 1) ci_real_ind(1, 2)], [1 1], 'Color', ilfcolor)
    
    ticklabels = [tractname([1, 2]), 'Ventral (avg)', tractname([4, 3]), 'Dorsal (avg)'];
    
end

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
print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_adultchild_' hemisphere '_singleshell']), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_adultchild' hemisphere '_singleshell']), '-depsc', '-r600')

hold off;

%% DORSAL - VENTRAL DIFFERENCE TEST

% test for difference betwen s_dorsal and s_ventral (bootstrapped samples of adult-child differences in dorsal and ventral fa).

% Get data to plot.
children = table2array(d(d.group_age3 ~= 3, toi));
adults = table2array(d(d.group_age3 == 3, toi));

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