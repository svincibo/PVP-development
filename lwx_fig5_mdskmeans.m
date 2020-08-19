clear all; close all; clc
format shortG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

group = 'adults'; % adults, children

hemisphere = 'both2'; % left, right, both, both2

n_clust = 3; % for kmeans, set n_clust = 3 for the separation of pathways hypothesis, set n_clust = 2 for the VTP closer to VH than DH hypothesis
n_rep = 10000; % for kmeans, number of clustering steps taken; chooses best based on total sum of distances between points and cluster centroids

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set working directories.
rootDir = '/Volumes/240/lwx/';

dorsalcolor = [224 83 114]/255; %salmon
ventralcolor= [189 80 199]/255; % purple
verticalcolor = [161 95 84]/255; % brown
darkgray = [150 150 150]/255;
gray = [175 175 175]/255;
lightgray = [200 200 200]/255;

linestyle = 'none';
marker = 'o'; markersize = 16;
fontname = 'Arial'; fontsizex = 16; fontsizey = 16; fontsizez = 16; fontangle = 'italic'; fontcolor = [0 0 0]; fontsmoothing = 'off';
yticklength = 0.05; xticklength = 0.05; zticklength = 0.05;

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

%% Multidimensional Scaling

% Calculate distance matrix from correlation matrix.
D = pdist(rho(2:end-1, 2:end-1), 'Euclidean');

% Select only the VV and VD dissimilarity scores and labels.
tractnames = {'slf12', 'slf3', 'mdlfang', 'mdlfspl', 'tpc', 'pArc', 'ilf', 'ifof'};

% Calculate eigenvalues.
[Y,eigvals] = cmdscale(D);

% Calculate normalized eigenvalues.
neigvals = eigvals/max(abs(eigvals));

% Display eigenvalues and normalized eigenvalues.
disp('eigenvals       eigvals/max(abs(eigvals))')
disp([eigvals neigvals]);

% Calculate error.
maxerr = max(abs(D - pdist(Y(:,1))))/max(D);
disp(['Max relative error for 1D: ' num2str(maxerr) '.']);

maxerr = max(abs(D - pdist(Y(:,1:2))))/max(D);
disp(['Max relative error for 2D: ' num2str(maxerr) '.']);

maxerr = max(abs(D - pdist(Y(:,1:3))))/max(D) ;
disp(['Max relative error for 3D: ' num2str(maxerr) '.']);

maxerr = max(abs(D - pdist(Y)))/max(D) ;
disp(['Max relative error for full reconstruction: ' num2str(maxerr) '.']);

linestyle = 'none';
marker = 'o';
markersize = 16;
lolimit = -0.50; hilimit = 0.77;

%ventral
v_idx = find(strcmp(tractnames, 'ilf') | strcmp(tractnames, 'ifof'));
%vertical
vp_idx = find(strcmp(tractnames, 'pArc') | strcmp(tractnames, 'tpc') | strcmp(tractnames, 'mdlfspl') | strcmp(tractnames, 'mdlfang'));
%dorsal
d_idx = find(strcmp(tractnames, 'slf12') | strcmp(tractnames, 'slf3'));

% figure(1)
% hold on;
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
% tractnames = {'SLF12', 'SLF3', 'MDLFang', 'MDLFspl', 'TPC', 'pArc', 'ILF', 'IFOF'};
% t = text(Y(:,1)-.05,Y(:,2),Y(:,3)-.08,tractnames, 'Color', 'k', 'FontSize', 16);
% t(1).Color = dorsalcolor; t(2).Color = dorsalcolor;
% t(3).Color = verticalcolor; t(4).Color = verticalcolor; t(5).Color = verticalcolor; t(6).Color = verticalcolor;
% t(7).Color = ventralcolor; t(8).Color = ventralcolor;
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
% labstr = {['First Dimension,\ e = ' num2str(eigvals(1), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% xax.Label.String = labstr;
% xax.Label.FontAngle = fontangle;
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
% labstr = {['Second Dimension,\ e = ' num2str(eigvals(2), '%2.3f')]};
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
% labstr = {['Third Dimension,\ e = ' num2str(eigvals(3), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% zax.Label.String = labstr;
% zax.Label.FontName = fontname;
% zax.Label.FontSize = fontsizez;
% 
% % Write.
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_' hemisphere '_' group '_mds_3D']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_' hemisphere '_' group '_mds_3D']), '-depsc')
% 
% figure(2)
% hold on;
% limits = 0.50;
% 
% plot(Y(v_idx, 1), Y(v_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
%     'MarkerFaceColor', ventralcolor, 'MarkerEdgeColor', ventralcolor);
% plot(Y(vp_idx, 1), Y(vp_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
%     'MarkerFaceColor', verticalcolor, 'MarkerEdgeColor', verticalcolor);
% plot(Y(d_idx, 1), Y(d_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
%     'MarkerFaceColor', dorsalcolor, 'MarkerEdgeColor', dorsalcolor);
% 
% tractnames = {'SLF12', 'SLF3', 'MDLFang', 'MDLFspl', 'TPC', 'pArc', 'ILF', 'IFOF'};
% t = text(Y(:,1)-0.05, Y(:,2)+0.08, tractnames, 'Color', 'k', 'FontSize', 16);
% t(1).Color = dorsalcolor; t(2).Color = dorsalcolor;
% t(3).Color = verticalcolor; t(4).Color = verticalcolor; t(5).Color = verticalcolor; t(6).Color = verticalcolor;
% t(7).Color = ventralcolor; t(8).Color = ventralcolor;
% 
% xlimhi = 1.1; xlimlo = -limits;
% xax = get(gca, 'xaxis');
% xax.Limits = [xlimlo xlimhi];
% xax.TickValues = [xlimlo ((xlimlo+xlimhi)/2)+(xlimlo/2) (xlimlo+xlimhi)/2 ((xlimlo+xlimhi)/2)+(xlimhi/2) xlimhi];
% xax.TickDirection = 'out';
% xax.TickLabels = {num2str(xlimlo, '%2.2f'), '', num2str((xlimlo+xlimhi)/2, '%2.2f'), '', num2str(xlimhi, '%2.2f')};
% xax.FontAngle = fontangle;
% xax.FontSize = fontsizex;
% labstr = {['First Dimension,\ e = ' num2str(eigvals(1), '%2.3f')]};
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
% labstr = {['Second Dimension,\ e = ' num2str(eigvals(2), '%2.3f')]};
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
% labstr = {['Third Dimension,\ e = ' num2str(eigvals(3), '%2.3f')]};
% labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
% zax.Label.String = labstr;
% zax.Label.FontName = fontname;
% zax.Label.FontSize = fontsizez;
% % Write.
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_' hemisphere '_' group '_mds_2D']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_' hemisphere '_' group '_mds_2D']), '-depsc')
% 
% hold off;
% 
% figure(4)
% bar(eigvals, 'FaceColor', 'k', 'EdgeColor', 'k', 'FaceAlpha', .5, 'EdgeAlpha', .5)
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
% ylimhi = 1.4; ylimlo = 0;
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
% print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_' hemisphere '_' group '_mds_eig']), '-dpng')
% print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_' hemisphere '_' group '_mds_eig']), '-depsc')
% 
%% kmeans clustering

% The first two eigenvalues explain the majority of the variance in the data, so we will apply kmeans clustering to the two vectors predicted by
% the two primary eigenvalues (n_dim = 2) to determine if the tracts in each pathway cluster together.

n_dim = 2; % can only handle n_dim = 2 at this point (and likely that is all that is needed for this project)

% for plotting
xlimhi = 0.9; xlimlo = -0.55;
ylimhi = 0.7; ylimlo = -0.55;

% Perform kmeans clustering on the predicted data.
[clust_idx, C, sumd, D] = kmeans(Y(:, 1:n_dim), n_clust, 'Replicates', n_rep, 'Options', statset('Display', 'final')); % kmeans uses squared Euclidean distances by default

% Make a mesh; compute the distance from each centroid, C, to points on the grid (procedure found here: https://www.mathworks.com/help/stats/kmeans.html).
x = xlimlo:0.01:xlimhi; 
y = ylimlo:0.01:ylimhi;
[xG,yG] = meshgrid(x, y);
XGrid = [xG(:), yG(:)]; % Defines a fine grid on the plot
idx2Region = kmeans(XGrid, n_clust, 'MaxIter', 1, 'Start', C);

% Calculate silhouette values.
s = silhouette(Y(:, 1:n_dim), clust_idx, 'sqEuclidean');

% Display silhouette values and means per cluster.
disp(tractnames); disp(s');
for n = 1:n_clust
    disp(['Mean silhouette value for cluster ' num2str(n) ' (1:darkgray, 2:gray, 3:lightgray), including ' num2str(length(find(clust_idx == n))) ' points: ' num2str(mean(s(find(clust_idx == n)))) '.']);
end

figure(5);
hold on;
markersize = 11; 
if n_clust == 3
    
    % Add grid.
    g = gscatter(XGrid(:,1), XGrid(:,2), idx2Region, [darkgray; gray; lightgray], '.');
    g(1).MarkerSize = markersize; g(2).MarkerSize = markersize; g(3).MarkerSize = markersize;
    
    % Add centroids.
    scatter(C(1, 1), C(1, 2), 'Marker', marker, 'MarkerFaceColor', darkgray, 'MarkerEdgeColor', darkgray, 'SizeData', 100)
    scatter(C(2, 1), C(2, 2), 'Marker', marker, 'MarkerFaceColor', gray, 'MarkerEdgeColor', gray, 'SizeData', 100)
    scatter(C(3, 1), C(3, 2), 'Marker', marker, 'MarkerFaceColor', lightgray, 'MarkerEdgeColor', lightgray, 'SizeData', 100)
    
elseif n_clust == 2
    
    % Add grid.
    g = gscatter(XGrid(:,1), XGrid(:,2), idx2Region, [gray; lightgray], '.');
    g(1).MarkerSize = markersize; g(2).MarkerSize = markersize;
    
    % Add centroids.
    scatter(C(1, 1), C(1, 2), 'Marker', marker, 'MarkerFaceColor', gray, 'MarkerEdgeColor', gray, 'SizeData', 100)
    scatter(C(2, 1), C(2, 2), 'Marker', marker, 'MarkerFaceColor', lightgray, 'MarkerEdgeColor', lightgray, 'SizeData', 100)
    
end

% Add points.
markersize = 16;
plot(Y(v_idx, 1), Y(v_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
    'MarkerFaceColor', ventralcolor, 'MarkerEdgeColor', ventralcolor);
plot(Y(vp_idx, 1), Y(vp_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
    'MarkerFaceColor', verticalcolor, 'MarkerEdgeColor', verticalcolor);
plot(Y(d_idx, 1), Y(d_idx, 2), 'LineStyle', linestyle, 'Marker', marker, 'MarkerSize', markersize, ...
    'MarkerFaceColor', dorsalcolor, 'MarkerEdgeColor', dorsalcolor);

% tractnames = {'SLF12', 'SLF3', 'MDLFang', 'MDLFspl', 'TPC', 'pArc', 'ILF', 'IFOF'};
% t = text(Y(:,1)-0.05,Y(:,2)+0.08,tractnames, 'Color', 'k', 'FontSize', 16);
% t(1).Color = dorsalcolor; t(2).Color = dorsalcolor;
% t(3).Color = verticalcolor; t(4).Color = verticalcolor; t(5).Color = verticalcolor; t(6).Color = verticalcolor;
% t(7).Color = ventralcolor; t(8).Color = ventralcolor;

xax = get(gca, 'xaxis');
xax.Limits = [xlimlo xlimhi];
xax.TickValues = [xlimlo ((xlimlo+xlimhi)/2)+(xlimlo/2) (xlimlo+xlimhi)/2 ((xlimlo+xlimhi)/2)+(xlimhi/2) xlimhi];
xax.TickDirection = 'out';
xax.TickLabels = {num2str(xlimlo, '%2.2f'), '', num2str((xlimlo+xlimhi)/2, '%2.2f'), '', num2str(xlimhi, '%2.2f')};
xax.FontAngle = fontangle;
xax.FontSize = fontsizex;
labstr = {['First Dimension,\ e = ' num2str(eigvals(1), '%2.3f')]};
labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
xax.Label.String = labstr;
xax.Label.FontAngle = fontangle;
xax.Label.FontName = fontname;
xax.Label.FontSize = fontsizex;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylimlo ylimhi];
yax.TickValues = [ylimlo ((ylimlo+ylimhi)/2)+(ylimlo/2) (ylimlo+ylimhi)/2 ((ylimlo+ylimhi)/2)+(ylimhi/2) ylimhi];
yax.TickDirection = 'out';
yax.TickLabels = {num2str(ylimlo, '%2.2f'), '', num2str((ylimlo+ylimhi)/2, '%2.2f'), '', num2str(ylimhi, '%2.2f')};
yax.TickDirection = 'out';
yax.TickLength = [yticklength yticklength];
yax.FontSize = fontsizey;
labstr = {['Second Dimension,\ e = ' num2str(eigvals(2), '%2.3f')]};
labstr = cellfun(@(x) strrep(x, '\', '\newline'), labstr, 'UniformOutput', false);
yax.Label.String = labstr;
yax.Label.FontName = fontname;
yax.Label.FontSize = fontsizey;

legend off
box off

% Write.
print(fullfile(rootDir, 'plots-singleshell', ['plot_fa_singleshell_' hemisphere '_' group '_mds_kmeans_k' num2str(n_clust)]), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_fa_singleshell_' hemisphere '_' group '_mds_kmeans_k' num2str(n_clust)]), '-depsc')

hold off;
























