%% function extract_ED_and_SD_from_ROI( which_ROI , subjnum, re_run  )
%
% This function extracts participants' data from some brain mask, then
% calculates euclidean distance (ED) and signed difference (SD).
%
% INPUTS:
%
% which_ROI [string] - which ROI to extract from? These are the ROIs in our
%                      Fig 3.
%                      Options: rV1, rLOC, rVO (visual ROIs 1-3)
%                               rBA41, rBA42, rBA22 (auditory ROIs 1-3)
%                               rBA3, rBA1, rBA2 (somatosensory ROIs 1-3)
%  
% subjnum [numeric] - this is the raw subject number in the master_dataset
%                     column rawSubjID      
%
% ----------------------------------------------------------
% Trial-by-trial predictions of subjective time from human brain activity,
% Sherman et al. 2020, biorxiv
% m.sherman@sussex.ac.uk
% ----------------------------------------------------------

function extract_ED_and_SD_from_mask( which_ROI , subjnum  , rerun, mode)

if nargin == 2; rerun = 0; end

%% change your paths!

path_ROI_masks = 'D:/ERICK/MD/ROIs for analysis'; % where are your masks kept?
path_fMRI      = 'C:/Users/erick/Documents/mestrado/subjects'; % this is where the 04_CISC... folders are kept with the preprocessed slice data
path_data = 'D:/ERICK/MD';
path_savedat   = fullfile(path_data, 'extracted_voxels', int2str(subjnum));
if ~isfolder(path_savedat)
    mkdir(path_savedat); % this is the folder where you'll save the data
end
path_savedat

%% parameters
TR = 0.8;

%% load the behavioural data
behaviour = readtable([path_data '/master_dataset.csv']);
behaviour = behaviour(behaviour.rawSubjID == subjnum,:);

%% which brain area are we extracting data from?
mask_file = [path_ROI_masks '/' which_ROI '.nii'];
functional_atlas = [path_ROI_masks '/rGlasser.nii'];

%% get the co-ordinates of that ROI
mask                          = spm_read_vols(spm_vol(mask_file),1);
fMask = spm_read_vols(spm_vol(functional_atlas),1);

%% get participant's folder
subj_str    = sprintf('%02d',subjnum);
% subj_str
subj_folder = dir([path_fMRI '/' subj_str '_CISC*']);
% subj_folder

if ~isempty(subj_folder) 
   disp('A pasta existe');
else
    disp('A pasta NÃO existe')
   return
end

cisc_file   = subj_folder.name;
subj_folder = [path_fMRI '/' cisc_file '/'];

sfile = [path_savedat '/' which_ROI '_corr.mat'];
disp(sfile)

if rerun || ~exist(sfile)
    
%% get participant's EPIs
EPI_folders = {};
run_label = {};
for iblock = 1:4
    epi_folder = dir([subj_folder 'fMRI_block_' num2str(iblock) '_0*']);
%     data_folder = [subj_folder data_folder.name '/'];
    epi_folder
%     epi_folder = dir([data_folder 'fMRI_block_' num2str(iblock) '_preprocessed*']);
%     epi_folder
    if ~isempty(epi_folder)
        EPI_folders = [EPI_folders; [subj_folder epi_folder.name] ];
        run_label = [run_label; iblock];
    end
end
EPI_folders
n_runs      = numel(EPI_folders);

%% for each folder, get the EPIs
EPI_files = cell(n_runs,1);
for irun = 1:n_runs
    epi_files        = dir([ EPI_folders{irun}, '/sw_uf*.nii']);
%     epi_files
    EPI_files{irun}  = arrayfun(@(x) [EPI_folders{irun} '/' x.name],epi_files,'UniformOutput',false);
end
% EPI_files

%% convert onsets to TRs
onsets_TR = cell(n_runs,1);
offsets_TR = cell(n_runs,1);
voxsum = cell(n_runs,1);
ed = cell(n_runs,1);
fc = cell(n_runs,1);

for irun = 1:n_runs

    clc;
%     EPI_folders{irun}
    sprintf('run %d of %d...',[irun n_runs])
    label = run_label{irun};
    which_ROI
    
    %behavioural_data = ['C:\Users\maxine\Dropbox\Projects\Sackler postdoc\fMRI & Time (Warrick)\data\' cisc_file '\eye_and_behaviour_' num2str(irun) '.mat'];
    %load(behavioural_data);
    
    % keep the behavioural data
    subjdata{irun} = behaviour(behaviour.run==label,:);
    
    % extract the onsets + offsets
    onsets  = subjdata{irun}.movieStartSecs;
    offsets = subjdata{irun}.movieStartSecs + subjdata{irun}.veridicalDuration;
    
    % convert to TR
    for itrial = 1:numel(offsets)
        onsets_TR{itrial} = floor(onsets(itrial)/TR);
        offsets_TR{itrial} = floor(offsets(itrial)/TR);
    end
    
%     subjdata{irun}
    subjdata{irun} = addvars(subjdata{irun}, onsets_TR, 'NewVariableNames', 'onsets_TR');
    subjdata{irun} = addvars(subjdata{irun}, offsets_TR, 'NewVariableNames', 'offsets_TR');
    subjdata{irun}
    t_corr = table2array(readtable([path_savedat '/time_corr.csv']));
%     t_corr
    t_corr = t_corr(irun, :);
%     t_corr(itrial)
    
    if strcmp(which_ROI, 'rGlasser')
         % loop trials and extract data
        for itrial = 1:numel(offsets)
            
            t_on  = subjdata{1, irun}.onsets_TR{itrial} - t_corr(itrial);
            t_off = subjdata{1, irun}.offsets_TR{itrial} - t_corr(itrial);
            
            % take those EPIs
    %         EPI_files{irun} 
            tEPI = EPI_files{irun}(t_on:t_off);
            
            avg_data = nan( 360 , numel(tEPI) );
%             idx = load('D:\ERICK\MD\ROIs for analysis\rGlasser_visual_index.txt');

            for node = 1:360
%                 node = idx(iter);
                % initialise data for the trial
                trial_data = nan( numel(tEPI) , numel(find(fMask==node)) );
            
                % for each EPI, extract the voxels
                for iepi = 1:numel(tEPI)
                    
                    % get the voxels of the current EPI
        %             iepi
                    current_epi = spm_read_vols(spm_vol(tEPI{iepi}));
                    
                    % get the signal from all voxels in the mask
                    trial_data(iepi,:) = current_epi( find(fMask == node) );
                end
    
                avg_data(node, :) = mean(trial_data, 2);
    
            end

%         %% compute euclidian distance
%         ed{itrial} =  sum(abs( diff(trial_data)),2);
%         
%         %% compute signed difference 
%         voxsum{itrial}  =  sum( diff(trial_data) ,2);

%         %% compute FC
%         fc_comp = triu(corr(avg_data'), 1);
%         fc_up = reshape(fc_comp, [], 1);
%         mascara_zeros = fc_up ~= 0;
%         fc{itrial} = fc_up(mascara_zeros);

            %% save time-series
            vd = table2array(subjdata{irun}(itrial, 'humanReport'));
    %         disp(vd);
            path_ts   = [path_savedat '/' 'time_series' '/' 'smoothed'];
            if ~isfolder(path_ts)
                mkdir(path_ts); % this is the folder where you'll save the data
            end
            tsfile = [path_ts '/' 'rGlasserTS_run' num2str(label) '_trial' num2str(itrial) '_dur=' num2str(vd) '.mat'];
            save(tsfile, 'avg_data');
        end
    else
         % loop trials and extract data
        for itrial = 1:numel(offsets)
            
            fc{itrial} = 'N/A';
            t_on  = subjdata{1, irun}.onsets_TR{itrial} - t_corr(itrial);
            t_off = subjdata{1, irun}.offsets_TR{itrial} - t_corr(itrial);
            
            % take those EPIs
    %         EPI_files{irun} 
            tEPI = EPI_files{irun}(t_on:t_off);
            
            % initialise data for the trial
            trial_data = nan( numel(tEPI) , numel(find(mask>0)) );
            
            % for each EPI, extract the voxels
            for iepi = 1:numel(tEPI)
                
                % get the voxels of the current EPI
    %             iepi
                current_epi = spm_read_vols(spm_vol(tEPI{iepi}));
                
                % get the signal from all voxels in the mask
                trial_data(iepi,:) = current_epi( find(mask > 0) );
            end
            
            %% compute euclidian distance
            ed{itrial} =  sum(abs( diff(trial_data)),2);
            
            %% compute signed difference 
            voxsum{itrial}  =  sum( diff(trial_data) ,2);
        end
    end
    
    subjdata{irun} = addvars(subjdata{irun}, voxsum, 'NewVariableNames', 'voxdiff');
    subjdata{irun} = addvars(subjdata{irun}, ed, 'NewVariableNames', 'ED');
    subjdata{irun} = addvars(subjdata{irun}, fc, 'NewVariableNames', 'FC');
    
end
save(sfile,'subjdata');
disp(['subj ' num2str(subjnum) ': saved in ' sfile])
end

end