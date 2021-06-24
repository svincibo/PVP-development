clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

multcompcorrection = 'no';

hemisphere = 'both2';% left, right, both, both2
group = 'children'; % adults, children

dorsalcolor = [236 176 32]/255; %burnt yellow
ventralcolor= [14 114 184]/255; % blue

%% READ IN DATA AND ORGANIZE.

% Read in 'final' data. This data should have all subjects removed, cleaned, etc.
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));

if strcmp(group, 'children')
    
    d = d(d.group_age3 ~= 3, :);
    
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
    
elseif strcmp(group, 'adults')
    
    % Include only adults; adult data needs to be z-scored.
    d = d(d.group_age3 == 3, :);
        
    % SELECT the measurements of the tracts that I care about.
    if strcmp(hemisphere, 'left')
        
        tpc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC')))); tpc = (tpc-nanmean(tpc))/nanstd(tpc);
        pArc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc')))); pArc = (pArc-nanmean(pArc))/nanstd(pArc);
        mdlfspl = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl')))); mdlfspl = (mdlfspl-nanmean(mdlfspl))/nanstd(mdlfspl);
        mdlfang = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang')))); mdlfang = (mdlfang-nanmean(mdlfang))/nanstd(mdlfang);
        vof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF')))); vof = (vof-nanmean(vof))/nanstd(vof);
        aslant = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant')))); aslant = (aslant-nanmean(aslant))/nanstd(aslant);
        slf12 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2')))); slf12 = (slf12-nanmean(slf12))/nanstd(slf12);
        slf3 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3')))); slf3 = (slf3-nanmean(slf3))/nanstd(slf3);
        ilf = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF')))); ilf = (ilf-nanmean(ilf))/nanstd(ilf);
        ifof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF')))); ifof = (ifof-nanmean(ifof))/nanstd(ifof);
        
    elseif strcmp(hemisphere, 'right')
        
        tpc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC')))); tpc = (tpc-nanmean(tpc))/nanstd(tpc);
        pArc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc')))); pArc = (pArc-nanmean(pArc))/nanstd(pArc);
        mdlfspl = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl')))); mdlfspl = (mdlfspl-nanmean(mdlfspl))/nanstd(mdlfspl);
        mdlfang = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang')))); mdlfang = (mdlfang-nanmean(mdlfang))/nanstd(mdlfang);
        vof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF')))); vof = (vof-nanmean(vof))/nanstd(vof);
        aslant = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant')))); aslant = (aslant-nanmean(aslant))/nanstd(aslant);
        slf12 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2')))); slf12 = (slf12-nanmean(slf12))/nanstd(slf12);
        slf3 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3')))); slf3 = (slf3-nanmean(slf3))/nanstd(slf3);
        ilf = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightILF')))); ilf = (ilf-nanmean(ilf))/nanstd(ilf);
        ifof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF')))); ifof = (ifof-nanmean(ifof))/nanstd(ifof);
        
    elseif strcmp(hemisphere, 'both')
        
        tpc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC') | strcmp(d.Properties.VariableNames, 'rightTPC')))); tpc = (tpc-nanmean(tpc))/nanstd(tpc);
        pArc = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc') | strcmp(d.Properties.VariableNames, 'rightpArc')))); pArc = (pArc-nanmean(pArc))/nanstd(pArc);
        mdlfspl = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl') | strcmp(d.Properties.VariableNames, 'rightMDLFspl')))); mdlfspl = (mdlfspl-nanmean(mdlfspl))/nanstd(mdlfspl);
        mdlfang = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang') | strcmp(d.Properties.VariableNames, 'rightMDLFang')))); mdlfang = (mdlfang-nanmean(mdlfang))/nanstd(mdlfang);
        vof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF') | strcmp(d.Properties.VariableNames, 'rightVOF')))); vof = (vof-nanmean(vof))/nanstd(vof);
        aslant = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant') | strcmp(d.Properties.VariableNames, 'rightAslant')))); aslant = (aslant-nanmean(aslant))/nanstd(aslant);
        slf12 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2') | strcmp(d.Properties.VariableNames, 'rightSLF1And2')))); slf12 = (slf12-nanmean(slf12))/nanstd(slf12);
        slf3 = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3') | strcmp(d.Properties.VariableNames, 'rightSLF3')))); slf3 = (slf3-nanmean(slf3))/nanstd(slf3);
        ilf = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF') | strcmp(d.Properties.VariableNames, 'rightILF')))); ilf = (ilf-nanmean(ilf))/nanstd(ilf);
        ifof = table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF') | strcmp(d.Properties.VariableNames, 'rightIFOF')))); ifof = (ifof-nanmean(ifof))/nanstd(ifof);
        
    elseif strcmp(hemisphere, 'both2')
        
        tpc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftTPC')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightTPC'))))); tpc = (tpc-nanmean(tpc))/nanstd(tpc);
        pArc = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftpArc')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightpArc'))))); pArc = (pArc-nanmean(pArc))/nanstd(pArc);
        mdlfspl = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFspl')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFspl'))))); mdlfspl = (mdlfspl-nanmean(mdlfspl))/nanstd(mdlfspl);
        mdlfang = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftMDLFang')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightMDLFang'))))); mdlfang = (mdlfang-nanmean(mdlfang))/nanstd(mdlfang);
        vof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftVOF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightVOF'))))); vof = (vof-nanmean(vof))/nanstd(vof);
        aslant = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftAslant')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightAslant'))))); aslant = (aslant-nanmean(aslant))/nanstd(aslant);
        slf12 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF1And2')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF1And2'))))); slf12 = (slf12-nanmean(slf12))/nanstd(slf12);
        slf3 = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftSLF3')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightSLF3'))))); slf3 = (slf3-nanmean(slf3))/nanstd(slf3);
        ilf = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftILF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightILF'))))); ilf = (ilf-nanmean(ilf))/nanstd(ilf);
        ifof = cat(1, table2array(d(:, find(strcmp(d.Properties.VariableNames, 'leftIFOF')))), table2array(d(:, find(strcmp(d.Properties.VariableNames, 'rightIFOF'))))); ifof = (ifof-nanmean(ifof))/nanstd(ifof);
        
    end
    
    % Compute correlations among tracts.
    [rho, pmat] = corr(cat(2, aslant, slf12, slf3, mdlfang, mdlfspl, tpc, pArc, ilf, ifof, vof), 'rows', 'pairwise');
    
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

disp(['z = ' num2str(z) ', p = ' num2str(p)]);

disp(['mean vv = ' num2str(mean(vv, 'all')) ', meanvd = ' num2str(mean(vd, 'all'))]);

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

load(fullfile(rootDir, 'supportFiles', 'redblue_colormap.mat'))
colormap(flipud(redblue(1:128, :)./255));
cb = colorbar; caxis([0 1]); cb.TickLength = 0.0000000001;
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

print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_corr_singleshell_wm_' hemisphere '_' multcompcorrection '_' group]), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_facorr_singleshell_wm_' hemisphere '_' multcompcorrection '_' group]), '-depsc', '-r600')

hold off;

%% PLOT DISTRIBUTIONS

figure(2)
hold on;
linewidth = 0.5;
linestyle = '-';
fontname = 'Arial';
fontsizex = 16; fontsizey = 16; fontsizez = 16;
fontangle = 'italic';
fontcolor = [0 0 0];
fontsmoothing = 'off';
yticklength = 0.05; xticklength = 0.05; zticklength = 0.05;
% dorsalcolor = [224 83 114]/255; %salmon
% ventralcolor= [189 80 199]/255; % purple
% verticalcolor = [161 95 84]/255; % brown
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

disp(['stdev vd = ' num2str(sigma)])

% VENTRAL DIST Compute the mean and standard deviation of the actual data.
mu=mean(dots2(:));
sigma=std(dots2(:));
% Put up lines and patch to indicate the mean, and mean +/- one standard deviation.
x = [mu-sigma mu-sigma mu+sigma mu+sigma];
y = [0 max(ylim) max(ylim) 0];
patch(x, y, ventralcolor, 'FaceAlpha', .2, 'FaceColor', ventralcolor, 'EdgeColor', ventralcolor, 'EdgeAlpha', .2)
line([mu, mu], [0 ymax], 'Color', ventralcolor, 'LineWidth', 1);

disp(['stdev vv = ' num2str(sigma)])

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
xax.Label.FontAngle = fontangle;
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
print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_verticaldorsal_verticalventral_hist_' hemisphere '_btsrp_' group]), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_verticaldorsal_verticalventral_hist_' hemisphere '_btsrp_' group]), '-depsc')

hold off;

% %% Multidimensional Scaling
% 
% % Calculate distance matrix from correlation matrix.
% D = pdist(rho(2:end-1, 2:end-1), 'Euclidean');
% 
% % Select only the VV and VD dissimilarity scores and labels.
% tractnames = {'slf12', 'slf3', 'mdlfang', 'mdlfspl', 'tpc', 'pArc', 'ilf', 'ifof'};
% 
% % Calculate eigenvalues.
% [Y,eigvals] = cmdscale(D);
% 
% % Calculate normalized eigenvalues.
% neigvals = eigvals/max(abs(eigvals));
% 
% % Display eigenvalues and normalized eigenvalues.
% disp('eigenvals       eigvals/max(abs(eigvals))')
% disp([eigvals neigvals]);
% 
% % Calculate error.
% maxerr = max(abs(D - pdist(Y(:,1))))/max(D);
% disp(['Max relative error for 1D: ' num2str(maxerr) '.']);
% 
% maxerr = max(abs(D - pdist(Y(:,1:2))))/max(D);
% disp(['Max relative error for 2D: ' num2str(maxerr) '.']);
% 
% maxerr = max(abs(D - pdist(Y(:,1:3))))/max(D) ;
% disp(['Max relative error for 3D: ' num2str(maxerr) '.']);
% 
% maxerr = max(abs(D - pdist(Y)))/max(D) ;
% disp(['Max relative error for full reconstruction: ' num2str(maxerr) '.']);
% 
% linestyle = 'none';
% marker = 'o';
% markersize = 16;
% lolimit = -0.50; hilimit = 0.77;
% 
% figure(3)
% hold on;
% %ventral
% v_idx = find(strcmp(tractnames, 'ilf') | strcmp(tractnames, 'ifof'));
% %vertical
% vp_idx = find(strcmp(tractnames, 'pArc') | strcmp(tractnames, 'tpc') | strcmp(tractnames, 'mdlfspl') | strcmp(tractnames, 'mdlfang'));
% %dorsal 
% d_idx = find(strcmp(tractnames, 'slf12') | strcmp(tractnames, 'slf3'));
% 
% p = plot3(Y(v_idx, 1), Y(v_idx, 2), Y(v_idx, 3), Y(vp_idx, 1), Y(vp_idx, 2), Y(vp_idx, 3), Y(d_idx, 1), Y(d_idx, 2), Y(d_idx, 3));
% p(1).LineStyle = linestyle;
% p(1).Marker = marker;
% p(1).MarkerSize = markersize;
% p(1).MarkerFaceColor = ventralcolor;
% p(1).MarkerEdgeColor = ventralcolor;
% 
% p(2).LineStyle = linestyle;
% p(2).Marker = marker;
% p(2).MarkerSize = markersize;
% p(2).MarkerFaceColor = verticalcolor;
% p(2).MarkerEdgeColor = verticalcolor;
% 
% p(3).LineStyle = linestyle;
% p(3).Marker = marker;
% p(3).MarkerSize = markersize;
% p(3).MarkerFaceColor = dorsalcolor;
% p(3).MarkerEdgeColor = dorsalcolor;
% 
% hold on;
% % tractnames = {'SLF12', 'SLF3', 'MDLFang', 'MDLFspl', 'TPC', 'pArc', 'ILF', 'IFOF'};
% % t = text(Y(:,1)-.05,Y(:,2),Y(:,3)-.08,tractnames, 'Color', 'k', 'FontSize', 16);
% % t(1).Color = dorsalcolor; t(2).Color = dorsalcolor;
% % t(3).Color = verticalcolor; t(4).Color = verticalcolor; t(5).Color = verticalcolor; t(6).Color = verticalcolor;
% % t(7).Color = ventralcolor; t(8).Color = ventralcolor;
% 
% axis equal;
% box off;
% view(-20, -5);
% 
% xlimhi = hilimit; xlimlo = lolimit;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% xax.TickValues = [xlimlo ((xlimlo+xlimhi)/2)+(xlimlo/2) (xlimlo+xlimhi)/2 ((xlimlo+xlimhi)/2)+(xlimhi/2) xlimhi];
% xax.TickDirection = 'out';
% xax.TickLabels = {num2str(xlimlo, '%2.2f'), '', num2str((xlimlo+xlimhi)/2, '%2.2f'), '', num2str(xlimhi, '%2.2f')};
% xax.FontAngle = fontangle;
% xax.FontSize = fontsizex;
% labstr = {['First Dimension,\ e = ' num2str(neigvals(1), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% xax.Label.String = labstr;
% yax.Label.FontAngle = fontangle;
% xax.Label.FontName = fontname;
% xax.Label.FontSize = fontsizex;
% 
% % yaxis
% ylimhi = hilimit; ylimlo = lolimit;
% yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
% yax.TickValues = [ylimlo ((ylimlo+ylimhi)/2)+(ylimlo/2) (ylimlo+ylimhi)/2 ((ylimlo+ylimhi)/2)+(ylimhi/2) ylimhi];
% yax.TickDirection = 'out';
% yax.TickLabels = {num2str(ylimlo, '%2.2f'), '', num2str((ylimlo+ylimhi)/2, '%2.2f'), '', num2str(ylimhi, '%2.2f')};
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.FontSize = fontsizey;
% labstr = {['Second Dimension,\ e = ' num2str(neigvals(2), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% yax.Label.String = labstr;
% yax.Label.FontName = fontname;
% yax.Label.FontSize = fontsizey;
% 
% % zaxis
% zlimhi = hilimit; zlimlo = lolimit;
% zax = get(gca,'zaxis');
% zax.Limits = [zlimlo zlimhi];
% zax.TickValues = [zlimlo ((zlimlo+zlimhi)/2)+(zlimlo/2) (zlimlo+zlimhi)/2 ((zlimlo+zlimhi)/2)+(zlimhi/2) zlimhi];
% zax.TickDirection = 'out';
% zax.TickLabels = {num2str(zlimlo, '%2.2f'), '', num2str((zlimlo+zlimhi)/2, '%2.2f'), '', num2str(zlimhi, '%2.2f')};
% zax.TickDirection = 'out';
% zax.TickLength = [zticklength zticklength];
% zax.FontSize = fontsizez;
% labstr = {['Third Dimension,\ e = ' num2str(neigvals(3), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% zax.Label.String = labstr;
% zax.Label.FontName = fontname;
% zax.Label.FontSize = fontsizez;
% 
% % Write.
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_mds_' hemisphere '_children_3D']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_mds_' hemisphere '_children_3D']), '-depsc')
% 
% figure(4)
% hold on;
% limits = 0.50;
% 
% p = plot3(Y(v_idx, 1), Y(v_idx, 2), Y(v_idx, 3), Y(vp_idx, 1), Y(vp_idx, 2), Y(vp_idx, 3), Y(d_idx, 1), Y(d_idx, 2), Y(d_idx, 3));
% p(1).LineStyle = linestyle;
% p(1).Marker = marker;
% p(1).MarkerSize = markersize;
% p(1).MarkerFaceColor = ventralcolor;
% p(1).MarkerEdgeColor = ventralcolor;
% 
% p(2).LineStyle = linestyle;
% p(2).Marker = marker;
% p(2).MarkerSize = markersize;
% p(2).MarkerFaceColor = verticalcolor;
% p(2).MarkerEdgeColor = verticalcolor;
% 
% p(3).LineStyle = linestyle;
% p(3).Marker = marker;
% p(3).MarkerSize = markersize;
% p(3).MarkerFaceColor = dorsalcolor;
% p(3).MarkerEdgeColor = dorsalcolor;
% 
% hold on;
% % tractnames = {'SLF12', 'SLF3', 'MDLFang', 'MDLFspl', 'TPC', 'pArc', 'ILF', 'IFOF'};
% % t = text(Y(:,1)-.05,Y(:,2),Y(:,3)-.08,tractnames, 'Color', 'k', 'FontSize', 16);
% % t(1).Color = dorsalcolor; t(2).Color = dorsalcolor;
% % t(3).Color = verticalcolor; t(4).Color = verticalcolor; t(5).Color = verticalcolor; t(6).Color = verticalcolor;
% % t(7).Color = ventralcolor; t(8).Color = ventralcolor;
% 
% view(0, 90); % to display only first and second dimension
% xlimhi = 1.1; xlimlo = -limits;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% xax.TickValues = [xlimlo ((xlimlo+xlimhi)/2)+(xlimlo/2) (xlimlo+xlimhi)/2 ((xlimlo+xlimhi)/2)+(xlimhi/2) xlimhi];
% xax.TickDirection = 'out';
% xax.TickLabels = {num2str(xlimlo, '%2.2f'), '', num2str((xlimlo+xlimhi)/2, '%2.2f'), '', num2str(xlimhi, '%2.2f')};
% xax.FontAngle = fontangle;
% xax.FontSize = fontsizex;
% labstr = {['First Dimension,\ e = ' num2str(neigvals(1), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% xax.Label.String = labstr;
% yax.Label.FontAngle = fontangle;
% xax.Label.FontName = fontname;
% xax.Label.FontSize = fontsizex;
% 
% % yaxis
% ylimhi = 0.7; ylimlo = -limits;
% yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
% yax.TickValues = [ylimlo ((ylimlo+ylimhi)/2)+(ylimlo/2) (ylimlo+ylimhi)/2 ((ylimlo+ylimhi)/2)+(ylimhi/2) ylimhi];
% yax.TickDirection = 'out';
% yax.TickLabels = {num2str(ylimlo, '%2.2f'), '', num2str((ylimlo+ylimhi)/2, '%2.2f'), '', num2str(ylimhi, '%2.2f')};
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.FontSize = fontsizey;
% labstr = {['Second Dimension,\ e = ' num2str(neigvals(2), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% yax.Label.String = labstr;
% yax.Label.FontName = fontname;
% yax.Label.FontSize = fontsizey;
% 
% % zaxis
% zlimhi = limits; zlimlo = -limits;
% zax = get(gca,'zaxis');
% zax.Limits = [zlimlo zlimhi];
% zax.TickValues = [zlimlo ((zlimlo+zlimhi)/2)+(zlimlo/2) (zlimlo+zlimhi)/2 ((zlimlo+zlimhi)/2)+(zlimhi/2) zlimhi];
% zax.TickDirection = 'out';
% zax.TickLabels = {num2str(zlimlo, '%2.2f'), '', num2str((zlimlo+zlimhi)/2, '%2.2f'), '', num2str(zlimhi, '%2.2f')};
% zax.TickDirection = 'out';
% zax.TickLength = [zticklength zticklength];
% zax.FontSize = fontsizez;
% labstr = {['Third Dimension,\ e = ' num2str(neigvals(3), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% zax.Label.String = labstr;
% zax.Label.FontName = fontname;
% zax.Label.FontSize = fontsizez;
% % Write.
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_mds_' hemisphere '_children_2D']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_mds_' hemisphere '_children_2D']), '-depsc')
% 
% hold off;
% 
% figure(5)
% bar(eigvals/max(abs(eigvals)), 'FaceColor', 'k', 'EdgeColor', 'k', 'FaceAlpha', .5, 'EdgeAlpha', .5)
% xlimhi = 8.5; xlimlo = 0.5;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% xax.TickDirection = 'out';
% xax.FontAngle = fontangle;
% xax.FontSize = fontsizex;
% labstr = {'Eigenvalue Number'};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% xax.Label.String = labstr;
% yax.Label.FontAngle = fontangle;
% xax.Label.FontName = fontname;
% xax.Label.FontSize = fontsizex;
% 
% % yaxis
% ylimhi = 1; ylimlo = 0;
% yax = get(gca,'yaxis');
% yax.Limits = [ylimlo ylimhi];
% yax.TickValues = [ylimlo (ylimlo+ylimhi)/2 ylimhi];
% yax.TickDirection = 'out';
% yax.TickLabels = {num2str(ylimlo, '%2.2f'), num2str((ylimlo+ylimhi)/2, '%2.2f'), num2str(ylimhi, '%2.2f')};
% yax.TickDirection = 'out';
% yax.TickLength = [yticklength yticklength];
% yax.FontSize = fontsizey;
% labstr = {'Normalized Eigenvalue'};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% yax.Label.String = labstr;
% yax.Label.FontName = fontname;
% yax.Label.FontSize = fontsizey;
% 
% box off;
% 
% % Write.
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_mds_' hemisphere '_children_eig']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_mds_' hemisphere '_children_eig']), '-depsc')



