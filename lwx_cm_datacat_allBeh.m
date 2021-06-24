% This script reads in FA, OD, and ICVF measures (from Brad Caron's
% TractProfiles App) for each of the tracts generated (from Dan Bullock's
% White Matter Segmentation App). It also reads in tract statistics (e.g.,
% number of streamlines for each tract (from Dan Bullock's Check Tract
% Quality App). It also reads in behavioral data collected as part of the
% LWX study.

% Read in tractprofiles and plot mean and standard deviations tract.

clear all; close all; clc
format shortG

measure = {'fa'}; % 'volume', 'fa', 'gmd', 'md'
roi = 'lobe'; %'lobe';

fontname = 'Arial';
fontsizex = 16; fontsizey = 12;
fontangle = 'italic';
fontsmoothing = 'off';
yticklength = 0;
xticklength = 0.05;
linewidth = .15;
alpha = .5;
save_figures = 'yes';

yc_color = [0.6350 0.0780 0.1840]; %red
oc_color = [0 0.4470 0.7410]; %blue
a_color = [0.41176 0.41176 0.41176]; %gray

% Set working directories.
rootDir = '/Volumes/240/lwx/';

% Get bl project foldername.
blprojectid = 'proj-5e849e65952fef3dcd7a1700';

% Read in lobe codes.
lobe = readtable([rootDir 'supportFiles/LUT_glasser.csv'], 'TreatAsEmpty', {'.', 'na'});

% Read in parcellation codes.
parcel = readtable([rootDir 'supportFiles/key_glasser.csv'], 'TreatAsEmpty', {'.', 'na'});

% Read in behavioral data.
beh = readtable([rootDir 'supportFiles/LWX_all_groupings.csv'], 'TreatAsEmpty', {'.', 'na'});

% Should outliers be removed? If so, which subIDs?
remove_outliers = 'yes';
if strcmp(remove_outliers, 'yes')
    
    % Identify outliers to be removed.
    %     outlier = [108 126 318];%
    % 108, snr is below 2 SD of group mean
    % 126, dwi image has major distortions, visual inspection
    % 318, snr is below 2 SD of group mean and dwi image has major distortions, visual inspection
    
    % Full set of outliers, including images removed based on visual inspection alone, being conservative to keep only the best data.
    outlier = [108 116 125 126 203 206 212 214 315 316 318];
    
end

% Should we include children only or all subjects?
include = 'all'; % options: childrenonly, all

% Parse table into array and header to make things easier.
beh_data_in = table2array(beh);
beh_data_in_header = beh.Properties.VariableNames;

% Identify outliers to be removed.
% 128 because WM measure z-scores are consistenly above z = +/-4.5 for all tracts
% 315 because strange WM in occipital lobe leading to no right or left VOF
% 318 because SNR is extremely low (snr = 3.925 with a z=-3.14)relative to others (range in snr = [7.4, 17.5]).
% outlier = [128 315 318];

for w = 1:length(measure)
    
    wm_measure = measure{w};
    
    %% TRACTOGRAPHY.
    
    % Get contents of the directory where the tract measures for this subject are stored.
    grp_contents = dir(fullfile(rootDir, blprojectid));
    
    % Remove the '.' and '..' files.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) ~= '.');
    
    % Keep only names that are subject folders.
    grp_contents = grp_contents(arrayfun(@(x) x.name(1), grp_contents) == 's');
    
    scount = 0;
    % Load in measures for this subject.
    for i = 1:size(grp_contents, 1)
        
        if strcmp(measure, 'fa') || strcmp(measure, 'md')
            
            % Get contents of the directory where the tract measures for this subject are stored.
            if strcmp(roi, 'lobe')
                
                sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-cortex_mapping_stats.id*', 'parc_MEAN.csv'));
                
            elseif strcmp(roi, 'tractendpoints')
                
                sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-cortex_mapping_stats.tag-tract_endpoints*', 'tracts_MEAN.csv'));
            
            end
        
        elseif strcmp(measure, 'volume')
            
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-SupraTentorial*', 'rois.csv'));
            
        elseif strcmp(measure, 'gmd')
            
            % Get contents of the directory where the tract measures for this subject are stored.
            sub_contents_measure = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-gray_matter_density*', 'parc_MEAN.csv'));
            
        end
        
        sub_contents_quality = dir(fullfile(grp_contents(i).folder, grp_contents(i).name,  '*parc-stats.tag-cortex_mapping_stats.tag-tract_endpoints*', 'tracts_COUNT_NONZERO.csv'));

        % Remove the '.' and '..' files.
        sub_contents_measure = sub_contents_measure(arrayfun(@(x) x.name(1), sub_contents_measure) ~= '.');
        sub_contents_quality = sub_contents_quality(arrayfun(@(x) x.name(1), sub_contents_quality) ~= '.');

        if ~isempty(sub_contents_measure)
            
            scount = scount + 1;
            
            % Grab subID.
            sub(scount) = str2num(grp_contents(i).name(end-2:end));
            
            % Display current sub ID.
            disp(grp_contents(scount).name)
            
            % Read in microstructural data for this subject.
            data_temp = readtable(fullfile(sub_contents_measure.folder, sub_contents_measure.name));
            
            % Read in quality data for this subject.
            quality_temp = readtable(fullfile(sub_contents_quality.folder, sub_contents_quality.name));
            
            % Change all measurements with less than 100 non-zero vertices to NaN.
            lo_quality_idx = find(table2array(quality_temp(:, 4)) < 500);
            data_temp{lo_quality_idx, 4:8} = NaN;
            
            % Change all measurements with snr <10 to NaN.
            lo_snr_idx = find(data_temp.snr<5);
            data_temp{lo_snr_idx, 4:8} = NaN;
            
            % Deal with idiosyncracies of the volume csv.
            if strcmp(measure, 'volume')
                
                % Remove the last 14 because those are subcortical.
                data_temp = data_temp(1:end-14, :);
                
                % Rename ROI_names.
                data_temp.ROI_name = str2double(erase(data_temp.ROI_name, 'ROI_'));
                for r = 1:size(data_temp, 1)
                    
                    name_temp(r) = parcel.ROI_name(find(parcel.ROI_number == data_temp.ROI_name(r)));
                    
                end
                name_temp = erase(name_temp, {'lh.', 'rh.', '.label'});
                data_temp.ROI_name = name_temp';
                
            end
            
            % If data_out exists, append; if not, create.
            if i == 1
                
                % Create data_out array.
                data_out = data_temp;
                
            else
                
                % Concatenate this array with the previous subject's array.
                data_out = cat(1, data_out, data_temp);
                
            end
            
            clear data_temp
            
        end
        
        clear sub_contents_tractprofiles
        
    end % end i

    if strcmp(roi, 'tractendpoints')
        
        % Select only ROIs from tracts that I care about.
        for t = 1:length(data_out.structureID)
            
            idx(t) = contains(data_out.structureID{t}, 'leftSLF1And2') || contains(data_out.structureID{t}, 'rightSLF1And2') ...
                || contains(data_out.structureID{t}, 'leftSLF3') || contains(data_out.structureID{t}, 'rightSLF3') ...
                || contains(data_out.structureID{t}, 'leftAslant') || contains(data_out.structureID{t}, 'rightAslant') ...
                || contains(data_out.structureID{t}, 'leftILF') || contains(data_out.structureID{t}, 'rightILF') ...
                || contains(data_out.structureID{t}, 'leftIFOF') || contains(data_out.structureID{t}, 'rightIFOF') ...
                || contains(data_out.structureID{t}, 'leftTPC') || contains(data_out.structureID{t}, 'rightTPC') ...
                || contains(data_out.structureID{t}, 'leftpArc') || contains(data_out.structureID{t}, 'rightpArc') ...
                || contains(data_out.structureID{t}, 'leftMDLFspl') || contains(data_out.structureID{t}, 'rightMDLFspl') ...
                || contains(data_out.structureID{t}, 'leftVOF') || contains(data_out.structureID{t}, 'rightVOF') ...
                || contains(data_out.structureID{t}, 'leftMDLFang') || contains(data_out.structureID{t}, 'rightMDLFang');
            
        end
        data_out = data_out(idx, :);
    end
    
      % Get a list of unique cortex names.
    if strcmp(measure, 'fa') || strcmp(measure, 'md') || strcmp(measure, 'gmd')
        list_roi = unique(data_out.structureID);
    elseif strcmp(measure, 'volume')
        list_roi = unique(data_out.ROI_name);
    end

    for y = 1:size(data_out, 1)
        
        name_temp = erase(data_out.structureID{y}, {'lh.', 'rh.', '_gaussian_1mm', '_FiberEndpoint'});
        data_out.structureID{y} = name_temp;
        
    end
    
    % Get measurements for each ROI (reorganizing so that each column is an roi).
    clear list_roi
    list_roi = unique(data_out.structureID);
    list_sub = unique(data_out.subjectID);
    for r = 1:size(list_roi, 1)
        
        for s = 1:length(list_sub)
                   
            % Find index of this tract for this subject.
            idx_rs = find(strcmp(data_out.structureID, list_roi{r}) & data_out.subjectID == list_sub(s));
               
            % Account for missing rois.
            if idx_rs
                
                wm_measures(s, r) = mean(data_out.fa(idx_rs));
                
            else
                
                wm_measures(s, r) = NaN;
                
            end
            
            
        end
        
    end
    
    % Convert all zeros to NaN;
    wm_measures(wm_measures==0) = NaN;

    % Append subID.
    wm_measures = cat(2, transpose(sub), wm_measures);
    wm_header = [{'subID'}, transpose(list_roi)];
    
    %     Create grouping and behavioral vectors.
    beh_out = cat(2, beh.SubjectID, beh.group_age, beh.Age_months, beh.group_lit, beh.c_lit, ...
        beh.group_vm, beh.c_vm, beh.group_fm, beh.c_fm, beh.Sex, beh.rBeeryVMI, beh.rBeeryVP, beh.rBeeryMC, beh.rgPegs_dom_forward, beh.rgPegs_nondom_forward, ...
        beh.rWJIV_LetterWordIdentification, beh.rWJIV_Spelling, beh.rWJIV_WordAttack, beh.rWJIV_SpellingOfSounds);
    
    % Determine which subIDs appear in both WM and BEH.
    sub_wm_beh = intersect(wm_measures(:, find(strcmp(wm_header, 'subID'))), beh.SubjectID);
    
    % Get indices of subjects who appear in both WM and BEH.
    sub_idx_wm = ismember(wm_measures(:, find(strcmp(wm_header, 'subID'))), sub_wm_beh);
    sub_idx_beh = ismember(beh.SubjectID, sub_wm_beh);
    
    % Select only subjects who appear in both WM and BEH.
    % Concatenate into one data array and one header array.
    % Remove redundant subID columns.
    data_all = cat(2, beh_out(sub_idx_beh, :), wm_measures(sub_idx_wm, find(strcmp(wm_header, 'subID'))+1:end));
    data_all_header = [{'subID',  'gp_age', 'cov_age', 'gp_lit', 'c_lit', 'gp_vm', 'c_vm', ...
        'gp_fm', 'c_fm', 'cov_sex', 'BeeryVMI', 'BeeryVP', 'BeeryMC', 'gPegs_dom_forward', 'gPegs_nondom_forward', ...
        'WJIV_LetterWordIdentification', 'WJIV_Spelling', 'WJIV_WordAttack', 'WJIV_SpellingOfSounds'}, wm_header{find(strcmp(wm_header, 'subID'))+1:end}];
    
    % Remove outliers.
    if strcmp(remove_outliers, 'yes') && exist('outlier')
        
        % Get index for outliers to be removed.
        idx_outlier = ismember(data_all(:, find(strcmp(data_all_header, {'subID'}))), outlier);
        
        % Remove outliers.
        data_all = data_all(~idx_outlier, :);
        
    end

    % Create the output table.
    data_tbl = array2table(cat(2, data_all), 'VariableNames', data_all_header);
    
    % Save all variables.
    save(fullfile(rootDir, 'supportFiles', ['LWX_CM_data_' wm_measure '_singleshell.mat']), 'data_tbl')
    
    % Save as a CSV file.
    writetable(data_tbl, fullfile(rootDir, 'supportFiles', ['LWX_CM_data_' wm_measure '_singleshell.csv']))
    
    % Reset for next loop.
    %     clearvars -except w rootDir beh_data_in_tbl beh_data_in_header beh blprojectid remove_outliers w_measures outlier
    
end
