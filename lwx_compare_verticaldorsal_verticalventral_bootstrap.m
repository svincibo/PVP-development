clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

multcompcorrection = 'no';

hemisphere = 'both2';% left, right, both, both2

%% READ IN DATA AND ORGANIZE.

% Read in 'final' data. This data should have all subjects removed, cleaned, etc.
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));

% Include only children.
d = d(d.group_age3 ~=3, :);

% SELECT the measurements of the tracts that I care about.
if strcmp(hemisphere, 'left')
    
    tpc = d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC_z')));
    pArc = d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc_z')));
    mdlfspl = d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl_z')));
    mdlfang = d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang_z')));
    vof = d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF_z')));
    aslant = d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant_z')));
    slf12 = d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2_z')));
    slf3 = d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3_z')));
    ilf = d(:, find(strcmp(d.Properties.VariableNames, 'leftILF_z')));
    ifof = d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF_z')));
    
elseif strcmp(hemisphere, 'right')
    
    tpc = d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC_z')));
    pArc = d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc_z')));
    mdlfspl = d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl_z')));
    mdlfang = d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang_z')));
    vof = d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF_z')));
    aslant = d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant_z')));
    slf12 = d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2_z')));
    slf3 = d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3_z')));
    ilf = d(:, find(strcmp(d.Properties.VariableNames, 'rightILF_z')));
    ifof = d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF_z')));
    
elseif strcmp(hemisphere, 'both')
    
    tpc = d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC_z') | strcmp(d.Properties.VariableNames, 'rightTPC_z')));
    pArc = d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc_z') | strcmp(d.Properties.VariableNames, 'rightpArc_z')));
    mdlfspl = d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl_z') | strcmp(d.Properties.VariableNames, 'rightMDLFspl_z')));
    mdlfang = d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang_z') | strcmp(d.Properties.VariableNames, 'rightMDLFang_z')));
    vof = d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF_z') | strcmp(d.Properties.VariableNames, 'rightVOF_z')));
    aslant = d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant_z') | strcmp(d.Properties.VariableNames, 'rightAslant_z')));
    slf12 = d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2_z') | strcmp(d.Properties.VariableNames, 'rightSLF1And2_z')));
    slf3 = d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3_z') | strcmp(d.Properties.VariableNames, 'rightSLF3_z')));
    ilf = d(:, find(strcmp(d.Properties.VariableNames, 'leftILF_z') | strcmp(d.Properties.VariableNames, 'rightILF_z')));
    ifof = d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF_z') | strcmp(d.Properties.VariableNames, 'rightIFOF_z')));
    
elseif strcmp(hemisphere, 'both2')
    
    tpc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC_z')))));
    pArc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc_z')))));
    mdlfspl = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl_z')))));
    mdlfang = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang_z')))));
    vof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF_z')))));
    aslant = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant_z')))));
    slf12 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2_z')))));
    slf3 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3_z')))));
    ilf = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightILF_z')))));
    ifof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF_z')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF_z')))));    
    
end

% Compute correlations among tracts.
if strcmp(hemisphere, 'both2')
    [rho, pmat] = corr(cat(2, aslant, slf12, slf3, mdlfang, mdlfspl, tpc, pArc, ilf, ifof, vof), 'rows', 'pairwise');

else
    [rho, pmat] = corr(table2array(cat(2, aslant, slf12, slf3, mdlfang, mdlfspl, tpc, pArc, ilf, ifof, vof)), 'rows', 'pairwise');
end

%% Bootstrap testing for difference between dorsal-vertical and dorsal-ventral correlations.

% z_aslant, z_slf12, z_slf3, z_mdlfang, z_mdlfspl, z_tpc, z_pArc, z_ilf, z_ifof, z_vof
G = [0 1 1 3 3 3 3 2 2 0];

% Get real distribution of vv-vd difference: Ventral-vertical correlation greater than dorsal-vertical correlation.
vv = rho(find(G==2), find(G==3));
vd = rho(find(G==1), find(G==3));
diff = vv - vd;
mu_diff=nanmean(diff(:));
sigma_diff=nanstd(diff(:));

% Get null distribution of vv-vd difference, where the correlation category (i.e., vertical-ventral, vertical-dorsal) labelling is broken.
null_dis = [vv(:); vd(:)];
for i = 1:1000
    
    % Randomly select a permutation with replacement.
    this_vv = randsample(null_dis, size(vv(:), 1), true);
    
    % Randomly select a permutation with replacement.
    this_vd = randsample(null_dis, size(vd(:), 1), true);
    
    % Get the correlation at that location.
    diff_null(i) = mean(this_vv) - mean(this_vd);
    
    clear this_vv this_vd
    
end
mu_diff_null=nanmean(diff_null);
sigma_diff_null=nanstd(diff_null);

% Test for significance.
z = (mu_diff - mu_diff_null)./sigma_diff_null;
p = 1-normcdf(abs(z), 0, 1);

%% PLOT CORRELATION PLOT
f1 = figure(1);
f1.InvertHardcopy = 'on';
hold on;

if strcmp(multcompcorrection, 'yes')
    
    mask = pmat < (.05/44);
    imagesc(rho.*mask, 'AlphaData', mask ~= 0 & eye(size(mask)) ~= 1);
    
else
    
    mask = ones(size(pmat));
    imagesc(rho.*mask, 'AlphaData', mask ~= 0);
    
end

% load(fullfile(rootDir, 'supportFiles', 'redblue_colormap.mat'))
% colormap(redblue./255)
colormap(parula)
cb = colorbar; caxis([-1 1]); cb.TickLength = 0.0000000001;
set(get(cb,'ylabel'),'string','correlation');
set(gca,'color', .75*[1 1 1]);
a = gca;
a.XLim = [0 size(rho, 2)]+.5;
a.YLim = [0 size(rho, 1)]+.5;
a.XTick = 1:size(rho, 2);
a.YTick = 1:size(rho, 1);
a.TickLength = [0 0];
a.YTickLabel = {'aslant', 'slf12', 'slf3', 'mdlfang', 'mdlfspl', 'tpc', 'pArc', 'ilf', 'ifof', 'vof'};
a.XTickLabel = {'aslant', 'slf12', 'slf3', 'mdlfang', 'mdlfspl', 'tpc', 'pArc', 'ilf', 'ifof', 'vof'};
a.XTickLabelRotation = 45;
title('Fractional Anisotropy');
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_corr_singleshell_wm_' hemisphere '_' multcompcorrection]), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_facorr_singleshell_wm_' hemisphere '_' multcompcorrection]), '-depsc', '-r600')

hold off;

%% PLOT DISTRIBUTIONS

figure(2)
hold on;
linewidth = 0.5;
linestyle = '-';
fontname = 'Arial';
fontsizex = 16; fontsizey = 10;
fontangle = 'italic';
fontcolor = [0 0 0];
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
dorsalcolor = [0.9290 0.6940 0.1250]; %burnt yellow
ventralcolor= [0 0.4470 0.7410]; % blue
ymax = .25;

% --- Correlations of dorsal with vertical
dots1 = rho(find(G==1), find(G==3));
h1 = histfit(dots1(:), numel(dots1));
N1=sum(h1(1).YData);
h1(1).YData=h1(1).YData/N1;
h1(2).YData=h1(2).YData/N1;
h1(1).FaceAlpha = 0; h1(1).FaceColor = 'w';
h1(1).EdgeAlpha = 0; h1(1).EdgeColor = 'w';
% h1(1).EdgeAlpha = 0.2; h1(1).EdgeColor = dorsalcolor;
h1(2).Color = dorsalcolor;
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel', yt/numel(dots1(:)))
% ylim([0 1])

% --- Correlations of ventral with vertical
dots2 = rho(find(G==2), find(G==3));
h2 = histfit(dots2(:), numel(dots2));
N2=sum(h2(1).YData);
h2(1).YData=h2(1).YData/N2;
h2(2).YData=h2(2).YData/N2;
h2(1).FaceAlpha = 0; h2(1).FaceColor = 'w';
h2(1).EdgeAlpha = 0; h2(1).EdgeColor = 'w';
% h2(1).EdgeAlpha = 0.2; h2(1).EdgeColor = ventralcolor;
h2(2).Color = ventralcolor;

% DORSAL DIST Compute the mean and standard deviation of the actual data.
mu=mean(dots1(:));
sigma=std(dots1(:));
% Put up lines and patch to indicate the mean, and mean +/- one standard deviation.
x = [mu-sigma mu-sigma mu+sigma mu+sigma];
y = [0 max(ylim) max(ylim) 0];
patch(x, y, dorsalcolor, 'FaceAlpha', .2, 'FaceColor', dorsalcolor, 'EdgeColor', dorsalcolor, 'EdgeAlpha', .2)
line([mu, mu], [0 ymax], 'Color', dorsalcolor, 'LineWidth', 1);

% VENTRAL DIST Compute the mean and standard deviation of the actual data.
mu=mean(dots2(:));
sigma=std(dots2(:));
% Put up lines and patch to indicate the mean, and mean +/- one standard deviation.
x = [mu-sigma mu-sigma mu+sigma mu+sigma];
y = [0 max(ylim) max(ylim) 0];
patch(x, y, ventralcolor, 'FaceAlpha', .2, 'FaceColor', ventralcolor, 'EdgeColor', ventralcolor, 'EdgeAlpha', .2)
line([mu, mu], [0 ymax], 'Color', ventralcolor, 'LineWidth', 1);

ylim([0 ymax]);

% Add text with z and p-value.
title(['zscore = ' num2str(z) ', p = ', num2str(p)]);

% % xaxis
xax = get(gca, 'xaxis');
xax.Limits = [-0.2 1];
xax.TickValues = [-0.2 0 0.5 1];
xax.TickDirection = 'out';
xax.TickLabels = {'-0.2', '0', '0.5', '1'};
xax.FontAngle = fontangle;
xax.Label.String = 'Correlation with Vertical White Matter';
yax.Label.FontAngle = fontangle;
xax.Label.FontName = fontname;
xax.Label.FontSize = fontsizex;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 ymax];
yax.TickValues = [0 ymax/2 ymax];
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.TickLabels = {'0', '0.1', '0.2'};
yax.Label.String = {'Relative Frequency'};
yax.Label.FontName = fontname;
yax.Label.FontSize = fontsizey;

box off;

pbaspect([6 2 6])

% Write.
print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_verticaldorsal_verticalventral_hist_' hemisphere '_btsrp']), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_verticaldorsal_verticalventral_hist_' hemisphere '_btsrp']), '-depsc')

hold off;
