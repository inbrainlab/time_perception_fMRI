
%% --------------------------------------------------
%  get_accumulated_changes( participantNumber , region , layer, analysis)
%
%  This is a function to accumulate salient changes from extracted summed voxel
%  data.
%
%  Inputs:
%          participantNumber: numeric
%          region:            'visual', 'auditory' or 'somatosensory'
%          layer:              1, 2 or 3
%          analysis:           'euclideanDistance' (confirmatory model)
%                              'signedDifference' (exploratory model)
%
%  Outputs:
%  Accumulation of Salient Events Predicts Subjective Time
%  Sherman et al. (2019) BiorXiv
%
%  Maxine Sherman
%  m.sherman@sussex.ac.uk
%  Adapted in 2024 (desouza.erickalmeida@gmail.com)


%  --------------------------------------------------
function [d,y] = get_accumulated_changes(subjnum, region, layer, analysis)
rng(11);
addpath('Results')

%% ========================================================================
%  Interpret inputs + get data
%  ========================================================================

% does the file exist? if not, make the mask
path_fMRI      = 'D:/ERICK/MD/';
path_savedat   = fullfile(path_fMRI, 'Results', int2str(subjnum), 'Robustness', region);
if ~isfolder(path_savedat)
    mkdir(path_savedat); % this is the folder where you'll save the data
end
subjpath = fullfile(path_fMRI, 'extracted_voxels', int2str(subjnum));
subjpath
mfile = [subjpath '/' region '_corr.mat'];
disp(mfile)

if isfile(mfile) 
   disp('A pasta existe');
else
    disp('A pasta NÃO existe')
   return
end

subjdata = load(mfile);
fieldnames(subjdata)
disp('mask file found...')

%% ========================================================================
%  Parameters & Setup: c = @(t) slope*exp( - decay_rate * t ) + baseline +
%  noise
%  ========================================================================

% tmax         = 1.5;
% tmin         = -1;
% decay_rate   = 1; % smaller values => steeper
% noise_sd     = 0.05;
% baseline = tmin;
% 
% slope    = tmax + baseline;

% parameters variation
tmax_list         = linspace(0,3,50);
baseline_list         = linspace(-3,0,50);
decay_rate   = 1; % smaller values => steeper
noise_sd     = 0.05;

%% ========================================================================
%  Loop through data + estimate
%  ========================================================================

for tmax = tmax_list
    for baseline = baseline_list
        y = struct('cVal',[],'isSalient',[],'zDiff',[], 'TR', [], 'isNan', []);
        d = struct('nSalient',[],'report',[],'isCity',[],'bias',[],'duration',[]);
        slope = tmax + baseline;
        for irun = 1:numel(subjdata.subjdata)
        
            df{irun} = subjdata.subjdata{irun};
            df{irun}
            
            % get all the data for the run
            run_data = struct('z',[],'t',[],'tr',[],'newTrial',[]);
            
            for itrial = 1:height(df{irun})
                   
                z                 = df{irun}.voxdiff(itrial);
                z = z{1};
                run_data.z        = [run_data.z ; z ];
                run_data.t        = [run_data.t ; repmat(itrial,numel(z),1)];
                run_data.tr       = [run_data.tr; [1:numel(z)]'];
                run_data.newTrial = [run_data.newTrial; [1; zeros(numel(z)-1,1)] ];
            end
            
            % zscore the run data
            run_data.z
            run_data.z = zscore(run_data.z); %run_data.z-mean(run_data.z);%
            
            for itrial = 1:height(df{irun})
        
                % get info
                idx = find(run_data.t == itrial);
                z   = run_data.z(idx);
              
                
                % get the relevant data and zscore so they're comparable across
                % trials/subjects
                nSalient = [];
                
                % initialise the criterion
                reset_criterion        = 0;
                t                      = 0;
        %         c = tmax;
                
                % loop through the TRs
                for itr = 1:numel(z)
                    
                    if z(itr) < 2.5
                        % if the distance on the prev trial was classed as 'salient', reset
                        % the decaying criterion
                        if reset_criterion; t = 1; else; t = t+1; end
                        
                        % get noise for this trial
                        noise  = normrnd( 0 , noise_sd );
                        
                        % get prev criterion
                        if t == 1; c = tmax;
                        else
                            c = y.cVal(end) - slope*exp( -t * decay_rate ) + baseline;
                        end
                       
                        y.cVal = [ y.cVal; c + noise ];
                        
                        % get the data for this trial
                        y.zDiff = [y.zDiff; z(itr)];
                        
                        % has the criterion been passed?
                        y.isSalient = [ y.isSalient; y.zDiff(end) >= y.cVal(end) ];
                        
                        % if it's been surpassed, reset
                        if y.isSalient(end)
                            reset_criterion = true;
                            nSalient        = [nSalient;1];
                        else
                            reset_criterion = false;
                            nSalient        = [nSalient;0];
                        end
                        nSalient
        
                        y.TR    = [ y.TR; itr]; % for plotting
                        y.isNan = [y.isNan; 0];
                        
                    else
                        y.isNan  = [y.isNan; 1];                   
                    end
                end
                
                % load stuff into the d structure
                d.nSalient = [d.nSalient; sum(nSalient)];
                d.report   = [d.report; df{irun}.humanReport(itrial)];
                d.duration = [d.duration; df{irun}.veridicalDuration(itrial)];
                d.isCity   = [d.isCity; df{irun}.isCity(itrial)];
                d.bias     = [d.bias; 100*df{irun}.humanBias(itrial)];
                   
            end
        end

        filename = ['accum_changes_' region '_tmax=' num2str(tmax) '_tmin=' num2str(baseline) '.mat'];
        dfile = [path_savedat '/' filename];
        save(dfile, "d")
    end
end
end
