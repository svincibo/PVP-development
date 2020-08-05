% This script reads in streamline count values (i.e., nfibers) and
% checks that the number of streamlines in each tract is correlated across
% subjects between sessions and checks that there are no significant
% differences in streamline count within tracts between groups at either session.
% A box plot is provided to help identify unusually low streamline counts within
% a particular tract and to compare sessions.

clear all; close all; clc
format shortG

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

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
    
else
    
    outlier = [];
    
end

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

% Read in behavioral data.
beh_data_in_tbl = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% Load in each tract's tractography measures for this subject.
sub_count = 0;
for t = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))))
        
        % Get contents of the directory where the tract measures for this subject are stored.
        sub_contents_tractstats = dir(fullfile(grp_contents(t).folder, grp_contents(t).name, '/dt-neuro-tractmeasures*/*.csv'));
        
        % Remove the '.' and '..' files.
        sub_contents_tractstats = sub_contents_tractstats(arrayfun(@(x) x.name(1), sub_contents_tractstats) ~= '.');
        
        if ~isempty(sub_contents_tractstats)
            
            % Display current sub ID.
            disp(grp_contents(t).name)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Read in data for this subject and this tract.
            data_tbl_in = readtable(fullfile(sub_contents_tractstats.folder, sub_contents_tractstats.name));
            
            % Convert into header for ease.
            data_all_in_header = data_tbl_in.TractName;
            
            % Get index matrices for hypothesis-driven grouping of WM tracts.
            for k = 1:length(data_all_in_header)
                
                % Indices of horizontal tracts.
                toi_idx(k) = strcmp(data_all_in_header{k}, 'leftSLF1And2') || strcmp(data_all_in_header{k}, 'rightSLF1And2') ...
                    || strcmp(data_all_in_header{k}, 'leftSLF3') || strcmp(data_all_in_header{k}, 'rightSLF3') ...
                    || strcmp(data_all_in_header{k}, 'leftILF') || strcmp(data_all_in_header{k}, 'rightILF') ...
                    || strcmp(data_all_in_header{k}, 'leftAslant') || strcmp(data_all_in_header{k}, 'rightAslant') ...
                    || strcmp(data_all_in_header{k}, 'leftTPC') || strcmp(data_all_in_header{k}, 'rightTPC') ...
                    || strcmp(data_all_in_header{k}, 'leftpArc') || strcmp(data_all_in_header{k}, 'rightpArc') ...
                    || strcmp(data_all_in_header{k}, 'leftMDLFspl') || strcmp(data_all_in_header{k}, 'rightMDLFspl') ...
                    || strcmp(data_all_in_header{k}, 'leftVOF') || strcmp(data_all_in_header{k}, 'rightVOF') ...
                    || strcmp(data_all_in_header{k}, 'leftMDLFang') || strcmp(data_all_in_header{k}, 'rightMDLFang') ...
                    || strcmp(data_all_in_header{k}, 'leftIFOF') || strcmp(data_all_in_header{k}, 'rightIFOF');
                
            end % end k
            
            % Get positions of tracts of interest in data_tbl.
            toi = find(toi_idx == 1);
            
            for t2 = 1:length(toi)
                
                % Read in mean stat.
                m(sub_count, t2) = data_tbl_in.StreamlineCount(toi(t2));
                
                % Grab tract name for grouping variable.
                tract{sub_count, t2} = char(data_tbl_in.TractName(toi(t2)));
                
                if m(sub_count, t2) <= 100
                    
%                     % Change entry to NaN so that it will not be included.
%                     m(t, t2) = NaN;
                      
                    % Alert user if tract has less than 300 streamlines.
                    disp(['Check data. Number of streamlines is ' num2str(m(sub_count, t2)) ' for ' tract{sub_count, t2} ' in subject ' grp_contents(t).name(5:7) '.'])
                    
                end % end if
                
            end % end t
            
            % Grab subID.
            sub(sub_count) = str2num(grp_contents(t).name(5:7));
            
            % Get age group.
            group(sub_count) = beh_data_in_tbl.group_age(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))));
            
            % Get age in months.
            age(sub_count) = beh_data_in_tbl.group_age(find((beh_data_in_tbl.SubjectID == str2num(grp_contents(t).name(5:7)))));
            
            clear data_temp toi toi_idx
            
        end % end if empty
        
    end % end if empty
    
end % end t

% Find empty cells.
idx = find(cellfun(@isempty,tract));

% Enter 'empty' in empty cells.
tract(idx) = {'empty'};

% Enter NaN for m in empty cells.
m(idx) = NaN;

% Get a list of unique tract names.
list_tract = tract(1, :);
list_tract = list_tract(~strcmp(list_tract, 'empty'));

% % Determine which subIDs appear in both WM and BEH.
sub_tract_beh = intersect(sub, beh_data_in_tbl.SubjectID);
%
% % Get indices of subjects who appear in both WM and BEH.
sub_idx_wm = ismember(sub, sub_tract_beh);
% sub_idx_beh = ismember(beh_data_in_tbl.No, sub_tract_beh);

% Select only subjects who appear in both WM and BEH.
% Concatenate into one data array and one header array.
% Remove redundant subID columns.
data_out = cat(2, sub',  group', age', m(sub_idx_wm, :));
data_out_header = [{'subID', 'group', 'cov_age'}, list_tract];

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(data_out(:, find(strcmp(data_out_header, {'subID'}))), outlier);
    
    % Remove outliers.
    data_out = data_out(~idx_outlier, :);
    
end

data_tbl = array2table(data_out, 'VariableNames', data_out_header);

% Save all variables.
streamlinecounts = data_tbl;
% % save([rootDir 'supportFiles/lwx_data_streamlinecount.mat'], 'streamlinecounts')

% Write out table.
writetable(data_tbl, fullfile(rootDir, 'supportFiles', 'lwx_data_streamlinecount.csv'));

% Group differences test
d = m(sub_idx_wm, :);
for tn = 1:size(d, 2)
    
%     disp(list_tract{tn});
    
    disp(['Are there streamline count differences between age groups for the ' list_tract{tn} '?'])
    tracttotest = d(:, tn);
    [p, tableout, stats] = anova1(tracttotest, group, 'off');
    disp(['F(' num2str(tableout{2, 3}) ', ' num2str(tableout{3, 3}) ') = ' num2str(tableout{2, 5}) ', p = ' num2str(tableout{2, 6}) '.'])
    
end % end tn

% Look at boxplot of streamline count.

% Convert into shorter format.
figure(1)
boxplot(d(group == 1, :), 'colors', yc_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
hold on
boxplot(d(group == 2, :), 'colors', oc_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
boxplot(d(group == 3, :), 'colors', a_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
xtickangle(90)
ylabel('Track Name')
xlabel('Streamline Count')
set(gca, 'XScale', 'log')
xlim([1 100000])
legend({'Younger Children', 'Older Children', 'Adults'}, 'Location', 'northwest')
legend box off
box off
print(fullfile(rootDir, 'plots-singleshell', 'plot_singleshell_boxplot_streamlinecount_allgroups'), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', 'plot_singleshell_boxplot_streamlinecount_allgroups'), '-depsc')

figure(2)
boxplot(d(group == 1, :), 'colors', yc_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
xtickangle(90)
ylabel('Track Name')
xlabel('Streamline Count')
set(gca, 'XScale', 'log')
xlim([1 100000])
title('Younger Children')
box off
print(fullfile(rootDir, 'plots-singleshell', 'plot_singleshell_boxplot_streamlinecount_yc'), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', 'plot_singleshell_boxplot_streamlinecount_yc'), '-depsc')

figure(3)
boxplot(d(group == 2, :), 'colors', oc_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
xtickangle(90)
title('Older Children')
ylabel('Track Name')
xlabel('Streamline Count')
set(gca, 'XScale', 'log')
xlim([1 100000])
box off
print(fullfile(rootDir, 'plots-singleshell', 'plot_singleshell_boxplot_streamlinecount_oc'), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', 'plot_singleshell_boxplot_streamlinecount_oc'), '-depsc')

figure(4)
boxplot(d(group == 3, :), 'colors', a_color, 'Labels', list_tract, 'PlotStyle', 'compact', 'Widths', 1, 'Orientation', 'horizontal')
xtickangle(90)
title('Adults')
ylabel('Track Name')
xlabel('Streamline Count')
set(gca, 'XScale', 'log')
xlim([1 100000])

box off
print(fullfile(rootDir, 'plots-singleshell', 'plot_singleshell_boxplot_streamlinecount_a'), '-dpng')
print(fullfile(rootDir, 'plots-singleshell', 'eps', 'plot_singleshell_boxplot_streamlinecount_a'), '-depsc')
