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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behavioral/demographic data.
    if ~isempty(find((beh_subjects == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
        % Remove the '.' and '..' files.
        sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
        
        if ~isempty(sub_contents_snr)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            m(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get standard deviation of SNR in X, Y, and Z directions.
            sd(sub_count) = std([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get training group.
            group(sub_count) = beh_group(find((beh_subjects == str2num(grp_contents(s).name(5:7)))));
            
            clear data_snr_temp get_temp
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

% Fix group so that it is just children and adult.
group(find(group == 2)) = 1; % change older children to "children = 1""
group(find(group == 3)) = 2; % then change adults to "adults = 2"

% Concatenate and sort according to group.
toplot = cat(2, subID', group', m', sd', b0');
toplot = sortrows(toplot, [2 3 1], 'ascend');

% cheap trick
clear subID group m sd b0
subID = toplot(:, 1)';
group = toplot(:, 2)';
m = toplot(:, 3)';
sd = toplot(:, 4)';
b0 = toplot(:, 5)';

% SNR plot
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
xticklength = 0;
alphablend = .8;

c_color = [0 0 0]; %[75 75 75]/255; % gray [.146 0 0]; % light black
a_color = [204 0 204]/255; %pink [178 34 34]/255; % firebrick red [0 .73 .73]; % turquoise

gscatter(m, 1:length(m), group, [a_color; c_color], '.', 20)
hold on;
gscatter(b0, 1:length(b0), group, [a_color; c_color], 'x', 8)
plot([15 15], [0 length(subID)+0.5], ':k')

for p = 1:length(m)
    
    if group(p) == 1
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', a_color)
        
    elseif group(p) == 2
        
        plot([m(p) - abs(sd(p)) m(p) + abs(sd(p))], [p p], 'Color', c_color)
              
    end
    
end

% xaxis
xax = get(gca, 'xaxis');
xax.Limits = [min([b0 m + sd])-5 max([b0 m + sd])+5];
xax.TickValues = floor(min([b0 m + sd])):10:ceil(max([b0 m + sd]));
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 length(subID)+0.5];
yax.TickValues = 1:length(subID);
yax.TickDirection = 'out';
yax.TickLabels = subID;
yax.FontName = fontname;
yax.FontSize = 8;

% general
a = gca;
%     a.TitleFontWeight = 'normal';
box off

legend({'Children', 'Adults'}, 'Location', 'southeast');
legend box off

a.XLabel.String = 'Signal-to-Noise Ratio (SNR)';
a.XLabel.FontSize = fontsize;
a.XLabel.FontAngle = fontangle;

a.YLabel.String = 'Participant ID';
a.YLabel.FontSize = fontsize;
a.YLabel.FontAngle = fontangle;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots-singleshell', 'plot_snr'), '-dpng')

hold off;