clear all; close all; clc
format shortG

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

wmoi = {'fa'};

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in behavioral data.
beh_in = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Read in wm measurement data.
load([rootDir 'supportFiles/LWX_data_fa_singleshell.mat'])
wm_in = data_tbl;

% Read in snr data.
snr_in = readtable([rootDir 'supportFiles/lwx_data_snr_singleshell.csv.csv'], 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - conservative removal.
    %         outlier = [108 126 214 318];
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
    
end

% Find subjects in wm_in that also appear in snr_in.
idx = ismember(wm_in.subID, snr_in.subID);

% Remove subjects from the wm_in data who do not appear in the snr_in data.
wm_in = wm_in(idx, :);

% Find subjects in snr_in that also appear in wm_in.
idx = ismember(snr_in.subID, wm_in.subID);

% Remove subjects from the snr_in data who do not appear in the wm_in data.
snr_in = snr_in(idx, :);

% Concatenate all variables of interest for ease.
d = cat(2, wm_in, table(snr_in.m));

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(d.subID, outlier);
    
    % Remove outliers.
    d = d(~idx_outlier, :);
    
else
    
    d = d;
    
end

% Get goup membership for ease.
group = d.gp_age;

% Reassign snr_in to snr for ease.
snr = d.Var1;

% Get table header for ease.
d_header = d.Properties.VariableNames;

% Plot on correlation plot for each tract of interest.
count = 0;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
xticklength = 0;
alphablend = .8;
for t = 1:size(d, 2)
    
    % Plot only columns in d that are tracts of interest measurements.
    if strcmp(d_header{t}, 'leftSLF1And2') || strcmp(d_header{t}, 'rightSLF1And2') ...
            || strcmp(d_header{t}, 'leftIFOF') || strcmp(d_header{t}, 'rightIFOF') ...
            || strcmp(d_header{t}, 'leftILF') || strcmp(d_header{t}, 'rightILF') ...
            || strcmp(d_header{t}, 'leftSLF3') || strcmp(d_header{t}, 'rightSLF3') ...
            || strcmp(d_header{t}, 'leftAslant') || strcmp(d_header{t}, 'rightAslant') ...
            || strcmp(d_header{t}, 'leftTPC') || strcmp(d_header{t}, 'rightTPC') ...
            || strcmp(d_header{t}, 'leftpArc') || strcmp(d_header{t}, 'rightpArc') ...
            || strcmp(d_header{t}, 'leftMDLFspl') || strcmp(d_header{t}, 'rightMDLFspl') ...
            || strcmp(d_header{t}, 'leftVOF') || strcmp(d_header{t}, 'rightVOF') ...
            || strcmp(d_header{t}, 'leftMDLFang') || strcmp(d_header{t}, 'rightMDLFang')
        
        % Update counter.
        count = count + 1;
        
        figure(count)
        hold on
        
        % Take only non-NaN.
        wm_tp = table2array(d(~isnan(table2array(d(:, t))), t));
        snr_tp = snr(~isnan(table2array(d(:, t))));
        group_tp = group(~isnan(table2array(d(:, t))));
       
        % Plot.
        gscatter(snr_tp, wm_tp, group_tp, [yc_color; oc_color; a_color], '.', 30)
        hold on;
        c1 = polyfit(snr_tp, wm_tp, 1);
        x1 = linspace(0, max(snr_tp)+5);
        f1 = polyval(c1, x1);
        plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', 'k')
        % hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2);
        % x2 = (1:size(f1, 2))'.*.30;
        % hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], a_color);
        % set(hp3, 'facecolor', a_color, 'edgecolor', 'none', 'facealpha', .2);
        
        tbltest = array2table(cat(2, snr_tp, wm_tp), 'VariableNames', {'snr', wmoi{1}});
        spec = 'snr ~ fa';
        out = fitlm(tbltest, spec)
        
        %  xaxis
        xax = get(gca, 'xaxis');
        xax.Limits = [10 max(snr_tp)+5];
        xax.TickValues = 10:5:max(snr_tp)+5;
        xax.TickDirection = 'out';
        xax.TickLength = [xticklength xticklength];
        xax.FontName = fontname;
        xax.FontSize = fontsize;
        xax.FontAngle = fontangle;
        
        % yaxis
        yax = get(gca,'yaxis');
%         yax.Limits = [0 max(wm_tp)+0.5];
        yax.Limits = [0.4 0.6];

        yax.TickValues = 1:1:ceil(max(wm_tp));
        yax.TickDirection = 'out';
        yax.TickLabels = 1:1:ceil(max(wm_tp));
        yax.FontName = fontname;
        yax.FontSize = 8;
        
        % general
        a = gca;
        %     a.TitleFontWeight = 'normal';
        box off
        
        legend({'Younger Children', 'Older Children', 'Adults', 'fit', 'sd'}, 'Location', 'northwest');
        legend box off
        
        a.XLabel.String = 'Signal-to-Noise Ratio (SNR)';
        a.XLabel.FontSize = fontsize;
        a.YLabel.String = 'Fractional Anisotropy (FA)';
        a.YLabel.FontSize = fontsize;
        pbaspect([1 1 1])
        
        title(d_header{t})
        
        print(fullfile(rootDir, 'plots', ['plot_snr_vs_' wmoi{1} '_' d_header{t}]), '-dpng')
        
        hold off
        
        clear tp snr_tp group_tp
        
    end % if strcmp
    
end % end for t




