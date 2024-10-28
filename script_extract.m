%preprocess

% clc
% subj = 45;
% preprocess_fmri(subj);

%EXTRACT DATA FROM fMRI 

clc

nodes = {'rBA1', 'rBA2', 'rBA3', 'rBA22', 'rBA41', 'rBA42', 'rLOC', 'rV1', 'rVO'};
% nodes = {'rNet1', 'rNet2', 'rNet3', 'rNet6', 'rNet7', 'rNet8', 'rNet11'};
% nodes = {'rNet1+11'};
% nodes = {'rGlasser'};
% nodes = {'rShen'};
% nodes = {'rLOC', 'rV1', 'rVO'}; %only visual
rerun = 0;
subj = 10:45;
mode = 'ts'; %'ts' for timesseries and 'se' for salient event

% for subjnum = subj
%     for ROI = 1:length(nodes)
%         which_ROI = nodes{ROI};
% %         extract_TS_from_mask(which_ROI , subjnum  , rerun, mode);
% %         extract_ED_and_SD_from_mask_RANDOM(which_ROI , subjnum  , rerun, mode);
% %         copy_task_volumes(which_ROI , subjnum  , rerun, mode);
% %         extract_TS_from_shen(which_ROI , subjnum  , rerun, mode);
%         extract_ED_and_SD_from_mask(which_ROI , subjnum  , rerun, mode);
% %         extract_ED_and_SD_from_mask_dubois(which_ROI, subjnum, rerun, mode);
%     end 
% end
    
analysis = 'signedChanges';
% GET ACCUM CHANGES
% análise hierárquica
regions = {'visual', 'auditory', 'somatosensory'};
layers = {1, 2, 3};

for subjnum = subj
    for region = 1:length(regions) 
        which_region = regions{region};
        for layer = 1:length(layers)
            which_layer = layers{layer};
            clc
            disp(which_region)
            get_accumulated_changes(subjnum, which_region , which_layer, analysis);
        end
    end
end

% análise de redes
% for subjnum = subj
%     for region = 1:length(nodes) %alterar pra cada caso (hierarquico/rede)
%         which_region = nodes{region};
%         clc
%         get_accumulated_changes_net(subjnum , which_region , analysis);
%     end
% end