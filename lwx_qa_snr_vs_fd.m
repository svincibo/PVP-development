clear all; close all; clc
format shortG

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

blprojectid = 'proj-5e61139282b37f2cfe8fdb28';

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed - conservative removal.
%         outlier = [108 126 212 214 318];
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 212, physical anomaly that precludes tracking of vertical tracks, visual inspection
    % 214, major motion artifacts, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Identify outliers to be removed - liberal removal.
    outlier = [108 116 119 125 126 206 212 214 303 317 318];
    % 116, FD > 2
    % 119, FD > 2
    % 125, FD > 2
    % 206, FD > 2 
    % 303, SNR < 15 
    % 317, FD > 2
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

% Get contents of the directory where the tract measures for this subject are stored.
grp_contents = dir(fullfile(rootDir, blprojectid));

% Remove the '.' and '..' files.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');

% Keep only names that are subject folders.
grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');

% For each subject.
sub_count = 0;
for s = 1:size(grp_contents, 1)
    
    % Only collect values for subjects that have both MRI and behaviora/demographic data.
    if ~isempty(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))))
        
        % Display current sub ID.
        disp(grp_contents(s).name)
        
        % Get contents of the directory where the SNR values for this subject are stored.
        sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
        % Remove the '.' and '..' files.
        sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
        
        % Only collect values for subjects that have diffusion data (noted by whether or not an snr file was created).
        t = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/product.json'));
        
        if ~isempty(sub_contents_motion) && ~isempty(t)
            
            % Update subject counter for when not all subjects are used/needed.
            sub_count = sub_count + 1;
            
            %% MOTION
            
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_motion = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-neuro-dtiinit*/*ecXform.mat'));
            % Remove the '.' and '..' files.
            sub_contents_motion = sub_contents_motion(arrayfun(@(x) x.name(1), sub_contents_motion) ~= '.');
            
            % Get motion parameters for this subject.
            load(fullfile(sub_contents_motion.folder,sub_contents_motion.name));
            
            % Select only the translation/rotation parameters.
            mot = vertcat(xform(:).ecParams);
            mot = mot(:, 1:6); % xyz translation xyz rotation
            
            % Motion parameters represent the translation/rotation that must occur
            % to return the image at timepoint tn to the place that it was at timepoint
            % t0. Thus, to calculate instantaneouse parameters, we need a moving
            % difference.
            for m = 1:size(mot, 2)
                
                % Get moving difference for each parameter. Append row of zeros for t1, per convention (Power et al., 2014).
                % This step creates the fd timecourses for each motion parameter for each run in the motion struct and represents instantaneous movement.
                movingdifference(:, m) = [0 ; diff(mot(:, m), 1, 1)]';
                
            end
            
            % Get an overall fd for all 6 parameters for each run.
            % This step creates the fd summary statistic for all 6 parameters for each timepoint in a run for each run (e.g., scalar FD timecourse).
            motion(sub_count, :) = sum(abs(movingdifference), 2)';
            
            %% SNR
            % Get contents of the directory where the SNR values for this subject are stored.
            sub_contents_snr = dir(fullfile(grp_contents(s).folder, grp_contents(s).name, '/dt-raw.tag-snr*/*product.json'));
            % Remove the '.' and '..' files.
            sub_contents_snr = sub_contents_snr(arrayfun(@(x) x.name(1), sub_contents_snr) ~= '.');
            
            % Get SNR for this subject.
            data_snr_temp = jsondecode(fileread([sub_contents_snr.folder filesep sub_contents_snr.name]));
            
            % Get SNR in b0 images.
            b0(sub_count) = str2num(data_snr_temp.SNRInB0_X_Y_Z{1});
            
            % Get mean SNR in X, Y, and Z directions.
            snr(sub_count) = mean([str2num(data_snr_temp.SNRInB0_X_Y_Z{2}), str2num(data_snr_temp.SNRInB0_X_Y_Z{3}), str2num(data_snr_temp.SNRInB0_X_Y_Z{4})]);
            
            %% DEMOGRAPHICS
            
            % Get subID.
            subID(sub_count) = str2num(grp_contents(s).name(5:7));
            
            % Get ylabel.
            lab{sub_count} = grp_contents(s).name;
            
            % Get session.
            age(sub_count) = beh.Age_months(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
            % Get training group.
            group(sub_count) = beh.group_age(find((beh.SubjectID == str2num(grp_contents(s).name(5:7)))));
            
        end % if ~isempty
        
    end % end if exist
    
end % end s

fd = mean(motion, 2);
sd = std(motion, 0, 2);
snr = snr';
group = group';

% Remove outliers.
if strcmp(remove_outliers, 'yes') && exist('outlier')
    
    % Get index for outliers to be removed.
    idx_outlier = ismember(subID, outlier);
    
    % Remove outliers.
    subID = subID(~idx_outlier);
    fd = fd(~idx_outlier);
    snr = snr(~idx_outlier);
    group = group(~idx_outlier);
    
    % Set figure note.
    ttlstr = 'SNR outlier removed.';
    
else
    
    % Set figure note.
    ttlstr = 'SNR outlier retained.';
    
end

% FD plot
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

gscatter(snr, fd, group, [yc_color; oc_color; a_color], '.', 30)
hold on;
c1 = polyfit(snr, fd, 1);
x1 = linspace(0,length(snr));
f1 = polyval(c1,x1);
plot(x1, f1, 'LineWidth', 3, 'LineStyle', '-', 'Color', 'k')
% hi = f1 + std(f1, 0, 2); lo = f1 - std(f1, 0, 2); 
% x2 = (1:size(f1, 2))'.*.30;
% hp3 = patch([x2; x2(end:-1:1); x2(1)], [lo'; hi(end:-1:1)'; lo(1)], a_color);
% set(hp3, 'facecolor', a_color, 'edgecolor', 'none', 'facealpha', .2);

tbltest = array2table(cat(2, snr, fd), 'VariableNames', {'snr', 'fd'});
spec = 'snr ~fd';
out = fitlm(tbltest, spec)

%  xaxis
xax = get(gca, 'xaxis');
xax.Limits = [0 max(snr)+5];
xax.TickValues = 0:5:max(snr)+5;
xax.TickDirection = 'out';
xax.TickLength = [xticklength xticklength];
xax.FontName = fontname;
xax.FontSize = fontsize;
xax.FontAngle = fontangle;

% yaxis
yax = get(gca,'yaxis');
yax.Limits = [0 max(fd)+0.5];
yax.TickValues = 1:1:ceil(max(fd));
yax.TickDirection = 'out';
yax.TickLabels = 1:1:ceil(max(fd));
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
a.YLabel.String = 'Framewise Displacement (FD)';
a.YLabel.FontSize = fontsize;
pbaspect([1 1 1])

print(fullfile(rootDir, 'plots', 'plot_snr_vs_fd'), '-dpng')

hold off;