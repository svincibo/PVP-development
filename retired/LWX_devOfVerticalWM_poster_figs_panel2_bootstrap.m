clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

multcompcorrection = 'no';

hemisphere = 'left';% left, right, both
if strcmp(hemisphere, 'right')
    opt = 10;
elseif strcmp(hemisphere, 'left')
    opt = 0;
elseif strcmp(hemisphere, 'both')
    opt = 20;
% else
%     opt = 0;
end

excludegroup = 3;
if excludegroup == 3
    include = 'childrenOnly';
elseif excludegroup == [1 2]
    include = 'adultsonly';
else
    include = 'all';
end

wm_measure = {'fa'}; %, 'md', 'ad', 'rd', 'od', 'icvf', 'isovf'};
beh_measure = {'age', 'lit', 'vm', 'fm'};
load(fullfile(rootDir, 'supportFiles', 'redblue_colormap.mat'))

%% Tractography

% Read in 'final' data. This data should have all subjects removed, cleaned, etc.
d = readtable(fullfile(rootDir, 'supportFiles', 'LWX_data_forMatlab_fa_singleshell.csv'));

% Include only children.
d = d(d.group_age3 ~=3, :);

% Get index matrices for hypothesis-driven grouping of WM tracts.
if strcmp(hemisphere, 'left')
    
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices of vertical tracts.
        TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftTPC'); 
        pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftpArc');
        MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl');
        MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFang');
        VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftVOF');
        
        aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftAslant');
        
        slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2');
        slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF3');
        
        ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftILF');
        ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftIFOF');
        
        % Indices of vertical tracts, z.
        z_TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftTPC_z'); 
        z_pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftpArc_z');
        z_MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl_z');
        z_MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFang_z');
        z_VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftVOF_z');
        
        z_aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftAslant_z');
        
        z_slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2_z');
        z_slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF3_z');
        
        z_ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftILF_z');
        z_ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftIFOF_z');
        
    end
    
elseif strcmp(hemisphere, 'right')
    
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices of vertical tracts.
        TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightTPC');
        pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightpArc');
        MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl');
        MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightMDLFang');
        VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightVOF');
        
        aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightAslant');
        
        slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2');
        slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightSLF3');
        
        ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightILF');
        ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightIFOF');
        
        % Indices of vertical tracts, z.
        z_TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightTPC_z'); 
        z_pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightpArc_z');
        z_MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl_z');
        z_MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightMDLFang_z');
        z_VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightVOF_z');
        
        z_aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightAslant_z');
        
        z_slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2_z');
        z_slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightSLF3_z');
        
        z_ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightILF_z');
        z_ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'rightIFOF_z');
        
    end
    
elseif strcmp(hemisphere, 'both')
    
    for k = 1:length(d.Properties.VariableNames)
        
        % Indices of vertical tracts.
        TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftTPC') || strcmp(d.Properties.VariableNames{k}, 'rightTPC');
        pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftpArc') || strcmp(d.Properties.VariableNames{k}, 'rightpArc');
        MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl');
        MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFang') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFang');
        VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftVOF') || strcmp(d.Properties.VariableNames{k}, 'rightVOF');
        
        aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftAslant') || strcmp(d.Properties.VariableNames{k}, 'rightAslant');
        
        slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2') || strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2');
        slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF3') || strcmp(d.Properties.VariableNames{k}, 'rightSLF3');
        
        ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftILF') || strcmp(d.Properties.VariableNames{k}, 'rightILF');
        ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftIFOF') || strcmp(d.Properties.VariableNames{k}, 'rightIFOF');
        
                % Indices of vertical tracts, z.
        z_TPC_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftTPC_z') || strcmp(d.Properties.VariableNames{k}, 'rightTPC_z');
        z_pArc_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftpArc_z') || strcmp(d.Properties.VariableNames{k}, 'rightpArc_z');
        z_MDLFspl_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFspl_z') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFspl_z');
        z_MDLFang_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftMDLFang_z') || strcmp(d.Properties.VariableNames{k}, 'rightMDLFang_z');
        z_VOF_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftVOF_z') || strcmp(d.Properties.VariableNames{k}, 'rightVOF_z');
        
        z_aslant_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftAslant_z') || strcmp(d.Properties.VariableNames{k}, 'rightAslant_z');
        
        z_slf12_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF1And2_z') || strcmp(d.Properties.VariableNames{k}, 'rightSLF1And2_z');
        z_slf3_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftSLF3_z') || strcmp(d.Properties.VariableNames{k}, 'rightSLF3_z');
        
        z_ilf_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftILF_z') || strcmp(d.Properties.VariableNames{k}, 'rightILF_z');
        z_ifof_idx(k) = strcmp(d.Properties.VariableNames{k}, 'leftIFOF_z') || strcmp(d.Properties.VariableNames{k}, 'rightIFOF_z');
        
    end
    
end

% SELECT the measurements of the tracts that I care (h and v) and
% the subjects that I care about (non-adults).
% Categorize into h or v. Convert all zeros to NaN.
tpc = d(:, TPC_idx); z_tpc = d(:, z_TPC_idx); 
pArc = d(:, pArc_idx); z_pArc =  d(:, z_pArc_idx); 
mdlfspl = d(:, MDLFspl_idx); z_mdlfspl = d(:, z_MDLFspl_idx); 
mdlfang = d(:, MDLFang_idx); z_mdlfang = d(:, z_MDLFang_idx); 
vof = d(:, VOF_idx); z_vof = d(:, z_VOF_idx);
aslant = d(:, aslant_idx); z_aslant = d(:, z_aslant_idx); 
slf12 = d(:, slf12_idx); z_slf12 = d(:, z_slf12_idx); 
slf3 = d(:, slf3_idx); z_slf3 = d(:, z_slf3_idx);
ilf = d(:, ilf_idx); z_ilf = d(:, z_ilf_idx);
ifof = d(:, ifof_idx); z_ifof = d(:, z_ifof_idx);

% Subselect only data for children and for the covariates of interest.
beh = array2table(cat(2, d.c_lit, d.c_vm, d.c_fm));

% Get measure-specific z-scores.
z_beh = array2table(cat(2, d.c_lit_z, d.c_vm_z, d.c_fm_z));

%% FIGURE 1

f1 = figure(1);
[rho, p] = partialcorr(table2array(cat(2, z_slf12, z_slf3, z_aslant, z_mdlfspl, z_mdlfang, z_tpc, z_pArc, z_vof, z_ilf, z_ifof, z_beh)), ...
    cat(2, d.cov_age, d.cov_sex), 'rows', 'pairwise');
 
p_wm = p(1:10, 1:10);
rho_wm = rho(1:10, 1:10);

if strcmp(multcompcorrection, 'yes')
    mask = p_wm < (.05/44);
else
    mask = ones(size(p_wm));
end
if strcmp(multcompcorrection, 'yes')
imagesc(rho_wm.*mask, 'AlphaData', mask ~= 0 & eye(size(mask)) ~= 1);
else
  imagesc(rho_wm.*mask, 'AlphaData', mask ~= 0);  
end
f1.InvertHardcopy = 'on';
hold on;
colormap(redblue./255)
cb = colorbar; caxis([-1 1]); cb.TickLength = 0.0000000001;
set(get(cb,'ylabel'),'string','correlation');
set(gca,'color', .75*[1 1 1]);
a = gca;
a.XTick = 1:size(rho_wm, 2);
a.YTick = 1:size(rho_wm, 1);
a.TickLength = [0 0];
a.YTickLabel = {'slf12', 'slf3', 'aslant', 'mdlfspl', 'mdlfang', 'tpc', 'pArc', 'vof', 'ilf', 'ifof'};
a.XTickLabel = {'slf12', 'slf3', 'aslant', 'mdlfspl', 'mdlfang', 'tpc', 'pArc', 'vof', 'ilf', 'ifof'};
a.XTickLabelRotation = 45;
title('Fractional Anisotropy');
pbaspect([1 1 1])

% figure(w+4)
% [rho, p] = partialcorr(cat(2, z_ilf, z_ifof, z_vof, z_tpc, z_pArc, z_mdlfspl, z_mdlfang, z_aslant, z_slf12, z_slf3, z_beh, age), sex);
% if strcmp(multcompcorrection, 'yes')
%     mask = p < (.05/45);
% else
%     mask = ones(size(p));
% end
% imagesc(rho.*mask);
% hold on;
% colorbar; caxis([-1 1]);
% a = gca;
% a.XTick = 1:15;
% a.YTick = 1:15;
% a.YTickLabel = {'ilf', 'ifof', 'vof', 'tpc', 'pArc', 'mdlfspl', 'mdlfang', 'aslant', 'slf12', 'slf3', 'literacy', 'visualmotor', 'finemotor', 'age', 'sex'};
% a.XTickLabel = {'ilf', 'ifof', 'vof', 'tpc', 'pArc', 'mdlfspl', 'mdlfang', 'aslant', 'slf12', 'slf3',  'literacy', 'visualmotor', 'finemotor', 'age', 'sex'};
% a.XTickLabelRotation = 45;
% title([wm_measure{w}]);

print(fullfile(rootDir, 'plots-singleshell', ['plot_mdlEstimates_fa_pcorr_singleshell_wm_' hemisphere '_' multcompcorrection]), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_mdlEstimates_fa_pcorr_singleshell_wm_' hemisphere '_' multcompcorrection]), '-depsc', '-r600')

hold off;

%% FIGURE 2

f2=figure(2); 
 
p_wmbeh = p(1:13, 11:13);
rho_wmbeh = rho(1:13, 11:13);

if strcmp(multcompcorrection, 'yes')
    mask = p_wmbeh < (.05/32);
else
    mask = ones(size(p_wmbeh));
end
if strcmp(multcompcorrection, 'yes')  
    imagesc(rho_wmbeh.*mask, 'AlphaData', mask ~= 0 & eye(size(mask)) ~= 1 & rho_wmbeh~=1);  
else
  imagesc(rho_wmbeh.*mask, 'AlphaData', mask ~= 0);  
end
f2.InvertHardcopy = 'on';
hold on;
colormap(redblue./255)
cb = colorbar; caxis([-1 1]); cb.TickLength = 0.0000000001;
set(get(cb,'ylabel'),'string','correlation');
set(gca,'color', .75*[1 1 1]);
a = gca;
a.XTick = 1:size(rho_wmbeh, 2);
a.YTick = 1:size(rho_wmbeh, 1);
a.TickLength = [0 0];
a.YTickLabel = {'slf12', 'slf3', 'aslant', 'mdlfspl', 'mdlfang', 'tpc', 'pArc', 'vof', 'ilf', 'ifof', 'literacy', 'visual-motor', 'fine-motor'};
a.XTickLabel = {'literacy', 'visual-motor', 'fine-motor'};
a.XTickLabelRotation = 45;
title('Fractional Anisotropy');
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', ['plot_mdlEstimates_fa_pcorr_singleshell_beh_' hemisphere '_' multcompcorrection]), '-dpng', '-r600')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_mdlEstimates_fa_pcorr_singleshell_beh_' hemisphere '_' multcompcorrection]), '-depsc', '-r600')

hold off;

% %% FIGURE 3
% %close all
% figure(3)
markercolor = [0 0.50196 0.50196]; % teal
capsize = 0;
marker = 'o';
linewidth = 0.5;
linestyle = '-';
markersize = 100;
fontname = 'Arial';
fontsizex = 16; fontsizey = 10;
fontangle = 'italic';
fontcolor = [0 0 0];
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
save_figures = 'yes';
alpha = .8;
% 
% % Plot correlation between TPC and literacy.
% x = z_beh(:, 3); behname = 'finemotor';
% y = z_tpc; tractname = 'tpc';
% z = cat(2, age, sex);
% s = scatter(y, x);
% s.SizeData = markersize;
% s.MarkerFaceColor = markercolor;
% s.MarkerEdgeColor = markercolor;
% hold on;
% [r_out, p_out, ~, ~] = plotpartialcorr2(y, x, z, markercolor);
% legend({['r = ' num2str(r_out, '%.3f')], ['p = ' num2str(p_out, '%.3f')]}, 'Location', 'northwest');
% legend('boxoff');
% 
% % yaxis
% xax = get(gca, 'yaxis');
% xax.Limits = [-2.05 2.05];
% xax.TickValues = [-2 0 2];
% xax.TickDirection = 'out';
% xax.TickLength = [yticklength yticklength];
% xax.TickLabels = {'-2', '', '2'};
% xax.Label.String = [behname ' (zscored)'];
% % xax.Label.FontAngle = fontangle;
% xax.Label.FontName = fontname;
% xax.Label.FontSize = fontsizey;
% 
% % xaxis
% yax = get(gca,'xaxis');
% yax.Limits = [-2.05 2.05];
% yax.TickValues = [-2.0 0 2.0];
% yax.TickDirection = 'out';
% yax.TickLength = [xticklength xticklength];
% yax.TickLabels = {'-2.0', '', '2.0'};
% yax.Label.String = {['FA in ' tractname ' (zscored)']};
% yax.FontAngle = fontangle;
% yax.Label.FontName = fontname;
% yax.Label.FontSize = fontsizex;
% 
% % general
% a = gca;
% %     a.TitleFontWeight = 'normal';
% box off
% 
% pbaspect([1 1 1])
% 
% %     pos=get(gca,'Position');
% %     pos1=pos-[0 .02 0 0];
% %     set(gca,'Position', pos1);
% 
% % Write.
% if strcmp(save_figures, 'yes')
%     
%     print(fullfile(rootDir, ['plots-' shell 'shell'], ['plot_anova_' wm_measure{w} '_' shell 'shell_' tractname '_' behname '_' hemisphere]), '-dpng')
%     print(fullfile(rootDir, ['plots-' shell 'shell'], 'eps', ['plot_anova_' wm_measure{w} '_ ' shell 'shell_' tractname '_' behname '_' hemisphere]), '-depsc')
%     
% end
% 
% hold off;
% 
%% FIGURE 4

% Bootstrap testing for difference between dorsal-vertical and dorsal-ventral correlations.

%z_slf12, z_slf3, z_aslant, z_mdlfspl, z_mdlfang, z_tpc, z_pArc, z_vof, z_ilf, z_ifof
% G = [1 1 0 3 3 3 3 3 2 2];
G = [1 1 0 3 3 3 3 0 2 2];

% Get real distribution: Ventral-vertical correlation greater than dorsal-vertical correlation.
vv = rho(find(G==2), find(G==3));
vd = rho(find(G==1), find(G==3));
dots_diff = vv - vd;
mu_diff=nanmean(dots_diff(:));
sigma_diff=nanstd(dots_diff(:));

% Get null distribution where tract to vertical, dorsal, ventral mapping is broken.
null_dis = [vv(:); vd(:)]; 
for iRand = 1:1000
    
    % Randomly select a permutation with replacement.
    this_vv = randsample(null_dis, size(vv(:), 1), true);
    
    % Randomly select a permutation with replacement.
    this_vd = randsample(null_dis, size(vd(:), 1), true);
    
    % Get the correlation at that location.
    diff_null(iRand) = mean(this_vv) - mean(this_vd);
    
    clear this_vv this_vd
    
end
mu_diff_null=nanmean(diff_null);
sigma_diff_null=nanstd(diff_null);

% Test for significance.
z = (mu_diff - mu_diff_null)./sigma_diff_null;
p = 1-normcdf(abs(z), 0, 1);

figure(4)
dorsalcolor = [0.9290 0.6940 0.1250]; %burnt yellow
ventralcolor= [0 0.4470 0.7410]; % blue
markersize = 50;
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
hold on;

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

% DORSAL SCATTER Compute the mean and standard deviation of the actual data.
mu=mean(dots1(:));
sigma=std(dots1(:));
% Put up lines and patch to indicate the mean, and mean +/- one standard deviation.
x = [mu-sigma mu-sigma mu+sigma mu+sigma];
y = [0 max(ylim) max(ylim) 0];
patch(x, y, dorsalcolor, 'FaceAlpha', .2, 'FaceColor', dorsalcolor, 'EdgeColor', dorsalcolor, 'EdgeAlpha', .2)
line([mu, mu], [0 ymax], 'Color', dorsalcolor, 'LineWidth', 1); 

% % Plot frequencies as scatter instead of bars.
% s1 = scatter(h1(1).XData, h1(1).YData);
% s1.SizeData = markersize;
% s1.MarkerFaceColor = dorsalcolor;
% s1.MarkerEdgeColor = dorsalcolor;

% VENTRAL SCATTER Compute the mean and standard deviation of the actual data.
mu=mean(dots2(:));
sigma=std(dots2(:));
% Put up lines and patch to indicate the mean, and mean +/- one standard deviation.
x = [mu-sigma mu-sigma mu+sigma mu+sigma];
y = [0 max(ylim) max(ylim) 0];
patch(x, y, ventralcolor, 'FaceAlpha', .2, 'FaceColor', ventralcolor, 'EdgeColor', ventralcolor, 'EdgeAlpha', .2)
line([mu, mu], [0 ymax], 'Color', ventralcolor, 'LineWidth', 1); 

ylim([0 ymax]);

% % Plot frequencies as scatter instead of bars.
% s2 = scatter(h2(1).XData, h2(1).YData);
% s2.SizeData = markersize;
% s2.MarkerFaceColor = ventralcolor;
% s2.MarkerEdgeColor = ventralcolor;

% Add legend.
% legend({['r = ' num2str(r_out, '%.3f')], ['p = ' num2str(p_out, '%.3f')]}, 'Location', 'northwest');
% legend('boxoff');

% Add text with z and p-value.
title(['zscore = ' num2str(z) ', p = ', num2str(p)]);
 
% % xaxis
xax = get(gca, 'xaxis');
xax.Limits = [-0.2 1];
xax.TickValues = [-0.2 0 0.5 1];
xax.TickDirection = 'out';
% xax.TickLength = [xticklength xticklength];
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
% 
% % general
% a = gca;
% %     a.TitleFontWeight = 'normal';
% box off
% 
pbaspect([6 2 6])
% 
% Write.
if strcmp(save_figures, 'yes')
    
    print(fullfile(rootDir, 'plots-singleshell', ['plot_anova_fa_singleshell_ventraldorsal_hist_' hemisphere '_btsrp']), '-dpng')
    print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_anova_af_singleshell_ventraldorsal_hist_' hemisphere '_btsrp']), '-depsc')
    
end
hold off;
