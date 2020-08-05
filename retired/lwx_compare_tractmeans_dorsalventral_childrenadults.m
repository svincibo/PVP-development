clear all; close all; clc
format shortG

% yc_color = [0.6350 0.0780 0.1840]; %red
% oc_color = [0 0.4470 0.7410]; %blue
% a_color = [0.41176 0.41176 0.41176]; %gray

alphastat = 0.01;

wmoi = {'fa'};
shell = 'singleshell'; % singleshell, multishell, multishell2sub
hemisphere = 'both';% left, right, both

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
if strcmp(shell, 'singleshell')
    
    blprojectid = 'proj-5e849e65952fef3dcd7a1700'; % singleshell
    
elseif strcmp(shell, 'multishell')
    
    blprojectif = 'proj-5a74d1e26ed91402ce400cca'; % multishell, denoised using model trained on 1 subs
    
elseif strcmp(shell, 'multishell2sub')
    
    blprojectif = 'proj-5ef88795c67a0de16e2a017f'; % multishell, denoised using model trained on 2 subs
    
end

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

for w = 1:length(wmoi)
    
    % Read in wm measurement data.
    load(fullfile(rootDir, ['supportFiles/LWX_data_fa_' shell '.mat']))
    wm = data_tbl;
    
    % Should outliers be removed? If so, which subIDs?
    remove_outliers = 'yes';
    if strcmp(remove_outliers, 'yes')
        
        % Identify outliers to be removed - conservative removal.
        %     outlier = [108 126 214 318];%
        % 108, snr is below 2 SD of group mean
        % 126, dwi image has major distortions, visual inspection
        % 214, major motion artifacts, visual inspection
        % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
        
        % Identify outliers to be removed - liberal removal.
        outlier = [108 116 119 125 126 206 214 303 317 318];
        % 116, FD > 2
        % 119, FD > 2
        % 125, FD > 2
        % 206, FD > 2
        % 303, SNR < 15
        % 317, FD > 2
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(wm.subID, outlier);
        
        % Remove outliers.
        wm = wm(~idx_outlier, :);
        
    end
    
    % Get indices of subjects whose white matter values are all NaN.
    idx_notnan = ~all(isnan(table2array(wm(:, 11:end))), 2);
    
    % Remove subjects whose white matter values are all NaN.
    wm = wm(idx_notnan, :);
    
    % Get table header for ease.
    wm_header = wm.Properties.VariableNames;
    
    if strcmp(hemisphere, 'both')
        
        % Get indices of columns that correspond to tracts of interest.
        toi = find(strcmp(wm_header, 'leftSLF1And2') + strcmp(wm_header, 'rightSLF1And2') ...
            + strcmp(wm_header, 'leftSLF3') + strcmp(wm_header, 'rightSLF3') ...
            + strcmp(wm_header, 'leftILF') + strcmp(wm_header, 'rightILF') ...
            + strcmp(wm_header, 'leftIFOF') + strcmp(wm_header, 'rightIFOF'));
        
        % Get indices of columns that correspond to dorsal tracts.
        dorsal = find(strcmp(wm_header, 'leftSLF1And2') + strcmp(wm_header, 'rightSLF1And2') ...
            + strcmp(wm_header, 'leftSLF3') + strcmp(wm_header, 'rightSLF3'));
        
        % Get indices of columns that correspond to ventral tracts.
        ventral = find(strcmp(wm_header, 'leftILF') + strcmp(wm_header, 'rightILF') ...
            + strcmp(wm_header, 'leftIFOF') + strcmp(wm_header, 'rightIFOF'));
        
    elseif strcmp(hemisphere, 'left')
        
        % Get indices of columns that correspond to tracts of interest.
        toi = find(strcmp(wm_header, 'leftSLF1And2') + strcmp(wm_header, 'leftSLF3')  ...
            + strcmp(wm_header, 'leftILF') + strcmp(wm_header, 'leftIFOF'));
        
        % Get indices of columns that correspond to dorsal tracts.
        dorsal = find(strcmp(wm_header, 'leftSLF1And2') + strcmp(wm_header, 'leftSLF3'));
        
        % Get indices of columns that correspond to ventral tracts.
        ventral = find(strcmp(wm_header, 'leftILF') + strcmp(wm_header, 'leftIFOF'));
        
    elseif strcmp(hemisphere, 'right')
        
        % Get indices of columns that correspond to tracts of interest.
        toi = find(strcmp(wm_header, 'rightSLF1And2') + strcmp(wm_header, 'rightSLF3')  ...
            + strcmp(wm_header, 'rightILF') + strcmp(wm_header, 'rightIFOF'));
        
        % Get indices of columns that correspond to dorsal tracts.
        dorsal = find(strcmp(wm_header, 'rightSLF1And2') + strcmp(wm_header, 'rightSLF3') );
        
        % Get indices of columns that correspond to ventral tracts.
        ventral = find(strcmp(wm_header, 'rightILF') + strcmp(wm_header, 'rightIFOF'));
        
    end
    
    
    %% INDIVIDUAL TRACTS
    
    % Bootstrapping for each tract individually.
    for t = 1:length(toi)
        
        % Get data to plot.
        children = table2array(wm(wm.gp_age ~= 3, toi(t)));
        adults = table2array(wm(wm.gp_age == 3, toi(t)));
        
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
        tractname{t} = wm_header{toi(t)};
        diff_real_ind(t) = nanmean(diff_real);
        ci_real_ind(t, :) = ci;
        
        tci = prctile(diff_null, 100*(1-alphastat));
        tci_real_ind(t) = tci;
        
        clear ci
        
    end
    
    
    %% DORSAL GROUP
    
    % Get data to plot.
    children = table2array(wm(wm.gp_age ~= 3, dorsal));
    adults = table2array(wm(wm.gp_age == 3, dorsal));
    
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
    children = table2array(wm(wm.gp_age ~= 3, ventral));
    adults = table2array(wm(wm.gp_age == 3, ventral));
    
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
    
    %     Plot
    
    dorsalcolor = [236 176 32]/255; %burnt yellow
    slf12color = [204 148 29]/255;
    slf3color = [255 229 173]/255;
    
    ventralcolor= [14 114 184]/255; % blue
    ifofcolor = [13 73 143]/255;
    ilfcolor = [94 188 255]/255;
    
    %'fa', 'ad', 'md', 'rd'
    lo = [-.01 0 0 .65];
    mid = [0.00 .03 .06 .75];
    hi = [0.10 .06 .12 .85];
    
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
    titlestring = 'Fractional Anisotropy (FA)';
    a.TitleFontWeight = 'normal';
    box off
    xlabel({['Difference ' titlestring ':']; 'Adults - Children'}, 'FontName', fontname, 'FontSize', fontsizex, 'FontAngle', fontangle, 'Color', fontcolor, 'FontSmoothing', fontsmoothing);
    a.XLabel.FontSize = fontsizex;
    pbaspect([1 1 1]);
    
    % Write.
    if strcmp(save_figures, 'yes')
        
        print(fullfile(rootDir, ['plots-' shell], ['plot_gpmeans_' hemisphere '_' shell '_' wmoi{w} '_diff_fig1']), '-dpng', '-r600')
        print(fullfile(rootDir, ['plots-' shell], 'eps', ['plot_gpmeans' hemisphere '_' shell '_' wmoi{w} '_diff_fig1']), '-depsc', '-r600')
        
    end
    
    hold off;
    
    %% DORSAL - VENTRAL DIFFERENCE TEST
    
    % test for difference betwen s_dorsal and s_ventral (bootstrapped samples of adult-child differences in dorsal and ventral fa).
    
    % Get data to plot.
    children = table2array(wm(wm.gp_age ~= 3, toi));
    adults = table2array(wm(wm.gp_age == 3, toi));
    
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
    
end

