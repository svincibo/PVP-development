% This script reads in SNR values and plots them according to session (pre-
% vs post-training) and group (expert=3, beginner=2, control=1).

clear all; close all; clc
format shortG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get subIDs present in beh.
beh_subjects = beh.SubjectID;

% Get groups from beh.
beh_group = beh.group_age;

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
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    % 116, FD > 2
    % 119, FD > 2
    % 125, FD > 2
    % 206, FD > 2 
    % 303, SNR < 15 
    % 317, FD > 2
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh_subjects  == str2num(grp_contents(s).name(5:7)))))
        
        % Only collect values for subjects that have diffusion data (noted by whether or not an snr file was created).
        t = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/product.json'));
        
        if numel(t) ~= 0
            %             exist(fullfile(grp_contents(s).folder, grp_contents(s).name, 'dt-raw.tag-snr-cc.id-5e8d9d8afa1b637b4f2f9350/product.json'))
            % Display current sub ID.
            disp(grp_contents(s).name)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
            % Remove the '.' and '..' files.
            sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get training group.
            group(sub_count) = beh_group(find((beh_subjects == str2num(grp_contents(s).name(5:7)))));
            
            clear data_snr_temp get_temp t
            
        end % if exist
        
    end % end if isempty
    
end % end t

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    b0 = b0(~idx_outlier);
    m = m(~idx_outlier);
    group = group(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

% Write out table for anova.
t_out = array2table(cat(2, subID', group', m', b0'), 'VariableNames', {'subID', 'group', 'm', 'b0'});
writetable(t_out, fullfile(rootDir, 'supportFiles', 'lwx_data_snr_singleshell.csv'));

% Group differences test
disp('Are there b0 SNR differences among groups?')
[p, tableout, stats] = anova1(b0, group, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between younger and older children?')
[h, p, ci stats] = ttest2(b0(group == 1), b0(group == 2));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between older children and adults?')
[h, p, ci stats] = ttest2(b0(group == 2), b0(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between younger children and adults?')
[h, p, ci stats] = ttest2(b0(group == 1), b0(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there b0 SNR difference between children and adults?')
[h, p, ci stats] = ttest2(b0(group ~= 3), b0(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

disp('Are there weighted SNR differences among groups?')
[p, tableout, stats] = anova1(m, group, 'off');
disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between younger and older children?')
[h, p, ci stats] = ttest2(m(group == 1), m(group == 2));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between older children and adults?')
[h, p, ci stats] = ttest2(m(group == 2), m(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between younger children and adults?')
[h, p, ci stats] = ttest2(m(group == 1), m(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Posthoc two-samples t-tests
disp('Is there weighted SNR difference between children and adults?')
[h, p, ci stats] = ttest2(m(group ~= 3), m(group == 3));
disp(['t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p) '.'])

% Group differences plot: b0
snr = b0;
figure(1)
hold on;
capsize = 0;
marker = 'o';
linewidth = 1.5;
linestyle = 'none';
markersize = 10;
fontname = 'Arial';
fontsize = 16;
fontangle = 'italic';
yticklength = 0;
xticklength = 0.05;
xtickvalues = [1 2 3];
alphablend = .8;
ylim_lo = 0;
ylim_hi = 50;

a_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
c_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

% Children
b1 = bar(1, nanmean(snr(group ~= 3)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group ~= 3)) - nanstd(snr(group ~= 3)) nanmean(snr(group ~= 3)) + nanstd(snr(group ~= 3))], 'Color', c_color)
% Adults
b2 = bar(2, nanmean(snr(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', a_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Children', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, b0 volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', ['plot_barplot_snr_b0_bygroup_outliers' remove_outliers]), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_barplot_snr_b0_bygroup_outliers' remove_outliers]), '-depsc')

hold off;
clear snr

% Group differences plot: weighted
snr = m;
figure(2)
hold on;

% Adults
b1 = bar(1, nanmean(snr(group == 1)), 'FaceColor', c_color, 'EdgeColor', c_color, 'FaceAlpha', alphablend);
plot([1 1], [nanmean(snr(group ~= 3)) - nanstd(snr(group ~= 3)) nanmean(snr(group ~= 3)) + nanstd(snr(group ~= 3))], 'Color', c_color)
% Adults
b2 = bar(2, nanmean(snr(group == 3)), 'FaceColor', a_color, 'EdgeColor', a_color, 'FaceAlpha', alphablend);
plot([2 2], [nanmean(snr(group == 3)) - nanstd(snr(group == 3)) nanmean(snr(group == 3)) + nanstd(snr(group == 3))], 'Color', a_color)

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0.5 2.5];
xax.TickValues = [1 2];
xax.TickDirection = 'out';
xax.TickLength = [yticklength yticklength];
xlabels = {'Children', 'Adults'};
xlabels = cellfun(@(x) strrep(x, ' ', '\newline'), xlabels, 'UniformOutput', false);
xax.TickLabels = xlabels;
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [ylim_lo ylim_hi];
yax.TickValues = [ylim_lo (ylim_lo+ylim_hi)/2 ylim_hi];
yax.TickDirection = 'out';
yax.TickLength = [xticklength xticklength];
yax.TickLabels = {num2str(ylim_lo, '%1.0f'), '', num2str(ylim_hi, '%1.0f')};
yax.FontName = fontname;
yax.FontSize = fontsize;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

a.YLabel.String = 'SNR, weighted volumes';

a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', ['plot_barplot_snr_weighted_bygroup_outliers' remove_outliers]), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', ['plot_barplot_snr_weighted_bygroup_outliers' remove_outliers]), '-depsc')

hold off;
