function extract_ED_and_SD_from_mask( which_ROI , subjnum  , rerun, mode)

if nargin == 2; rerun = 0; end

%% change your paths!

path_ROI_masks = 'D:/ERICK/MD/ROIs for analysis'; % where are your masks kept?
path_fMRI      = 'C:/Users/erick/Documents/mestrado/subjects'; % this is where the 04_CISC... folders are kept with the preprocessed slice data
path_data = 'D:/ERICK/MD';
base_path = 'C:/Users/erick/Documents/mestrado';
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
subj_folder = dir([path_fMRI '/' subj_str '_CISC*']);

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

% if rerun || ~exist(sfile)
    
%% get participant's EPIs
EPI_folders = {};
run_label = {};
for iblock = 1:4
    epi_folder = dir([subj_folder 'fMRI_block_' num2str(iblock) '_0*']);
%     data_folder = [subj_folder data_folder.name '/'];
    epi_folder
%     epi_folder = dir([data_folder 'fMRI_block_' num2str(iblock) '_preprocessed*']);
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
    EPI_files{irun}  = arrayfun(@(x) [EPI_folders{irun} '/' x.name],epi_files,'UniformOutput',false);
end

%% convert onsets to TRs
onsets_TR = cell(n_runs,1);
offsets_TR = cell(n_runs,1);
voxsum = cell(n_runs,1);
ed = cell(n_runs,1);
fc = cell(n_runs,1);

for irun = 1:n_runs

    clc;
    sprintf('run %d of %d...',[irun n_runs])
    label = run_label{irun};
    which_ROI
    
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
    
    subjdata{irun} = addvars(subjdata{irun}, onsets_TR, 'NewVariableNames', 'onsets_TR');
    subjdata{irun} = addvars(subjdata{irun}, offsets_TR, 'NewVariableNames', 'offsets_TR');
    subjdata{irun}
    t_corr = table2array(readtable([path_savedat '/time_corr.csv']));
    t_corr = t_corr(irun, :);
    
    if strcmp(which_ROI, 'rGlasser')
         
        % take those EPIs
        tEPI = EPI_files{irun};
        path_tEPI   = fullfile(base_path, 'extracted_voxels', int2str(subjnum), 'whole_TS', 'smoothed');
        if ~isfolder(path_tEPI)
            mkdir(path_tEPI)
        end
        
        avg_data = nan( 360 , numel(tEPI) );
%             idx = load('D:\ERICK\MD\ROIs for analysis\rGlasser_visual_index.txt');

        for node = 1:360
%                 node = idx(iter);
            % initialise data for the trial
            full_data = nan( numel(tEPI) , numel(find(fMask==node)));

            % for each EPI, extract the voxels
            for iepi = 1:numel(tEPI)
                
                % get the voxels of the current EPI
                current_epi = spm_read_vols(spm_vol(tEPI{iepi}));
                
                % get the signal from all voxels in the mask
                full_data(iepi,:) = current_epi( find(fMask == node) );
            end

            avg_data(node, :) = mean(full_data, 2);

        end

        tsfile = [path_tEPI '/' 'rGlasserTS_run' num2str(label) '.mat'];
        save(tsfile, 'avg_data');

    else
        path_ts   = [path_savedat '/' 'time_series' '/' 'smoothed'];
        if ~isfolder(path_ts)
            mkdir(path_ts); % this is the folder where you'll save the data
        end

        for itrial = 1:numel(offsets)
            
            t_on  = subjdata{1, irun}.onsets_TR{itrial} - t_corr(itrial);
            t_off = subjdata{1, irun}.offsets_TR{itrial} - t_corr(itrial);
       
            tEPI = EPI_files{irun}(t_on:t_off);
            
            avg_data = nan( 1 , numel(tEPI) );

            trial_data = nan( numel(tEPI) , numel(find(fMask>0)) );
    
            % for each EPI, extract the voxels
            for iepi = 1:numel(tEPI)
                
                % get the voxels of the current EPI
                current_epi = spm_read_vols(spm_vol(tEPI{iepi}));
                
                % get the signal from all voxels in the mask
                trial_data(iepi,:) = current_epi( find(fMask > 0) );
            end
    
            avg_data(1, :) = mean(trial_data, 2);

            %% save time-series
            vd = table2array(subjdata{irun}(itrial, 'humanReport'));
    
            tsfile = [path_ts '/' which_ROI '_TS_run' num2str(label) '_trial' num2str(itrial) '_dur=' num2str(vd) '.mat'];
            save(tsfile, 'avg_data');
        end
        end
    end
% end
end
