
%% --------------------------------------------------
%  get_accumulated_changes( participantNumber , region , layer, analysis)
%  --------------------------------------------------
function [d,y] = get_accumulated_changes_RANDOM(subjnum, region , layer, analysis)
rng(11);
addpath('Results')

%% ========================================================================
%  Interpret inputs + get data
%  ========================================================================

switch region
    case 'visual'
        mask = { 'rV1' , 'rLOC' , 'rVO' };
        mask = mask{ layer };
    case 'auditory'
        mask = { 'rBA41' , 'rBA42' , 'rBA22' };
        mask = mask{ layer };
    case 'somatosensory';
        mask = { 'rBA3' , 'rBA1' , 'rBA2' };
        mask = mask{ layer };
end

% does the file exist? if not, make the mask
path_fMRI      = 'D:/ERICK/MD/';
path_savedat   = fullfile(path_fMRI, 'Results', int2str(subjnum));
if ~isfolder(path_savedat)
    mkdir(path_savedat); % this is the folder where you'll save the data
end
subjpath = fullfile(path_fMRI, 'extracted_voxels', int2str(subjnum), 'random_changes');
mfile = [subjpath '/' mask '_corr_RANDOM.mat'];
disp(mfile)
subjdata = load(mfile);
fieldnames(subjdata)
disp('mask file found...')

%% ========================================================================
%  Parameters & Setup: c = @(t) slope*exp( - decay_rate * t ) + baseline +
%  noise
%  ========================================================================

switch layer
    
    case 1 % most changes here
        tmax         = 1.5;
        tmin         = -1;
        decay_rate   = 1; % smaller values => steeper
        noise_sd     = 0.05;
        
    case 2 % middling changes
        tmax         = 1.5;
        tmin         = -0.5;
        decay_rate   = 1; % smaller values => steeper
        noise_sd     = 0.05;
        
    case 3 % fewest changes
        tmax         = 1.5;
        tmin         = 0;
        decay_rate   = 1; % smaller values => steeper
        noise_sd     = 0.05;
        
end

baseline = tmin;
slope    = tmax + baseline;


%% ========================================================================
%  Initialise structures
%  ========================================================================

y = struct('cVal',[],'isSalient',[],'zDiff',[], 'TR', [], 'isNan', []);
d = struct('nSalient',[],'report',[],'isCity',[],'bias',[],'duration',[]);

%% ========================================================================
%  Loop through data + estimate
%  ========================================================================

for irun = 1:numel(subjdata.subjdata)

    df{irun} = subjdata.subjdata{irun};
    df{irun}
    
    % get all the data for the run
    run_data = struct('z',[],'t',[],'tr',[],'newTrial',[]);
    
    for itrial = 1:height(df{irun})

        switch analysis
            case 'euclideanDistance'
                z                 = subjdata{irun,itrial}.euclideanDistance;
            case 'signedChanges'
                z                 = df{irun}.voxdiff(itrial);
                z = z{1};
        end
       
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
        
        % loop through the TRs
        for itr = 1:numel(z)
            
            if z(itr) < 2.5
                % if the distance on the prev trial was classed as 'salient', reset
                % the decaying criterion
                if reset_criterion; t = 1; else; t = t+1; end
                
                % get noise for this trial
                noise  = normrnd( 0 , noise_sd );
                noise
                
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

    filename = ['accum_changes_' region '_' int2str(layer) '_RANDOM.mat'];
    dfile = [path_savedat '/' filename];
    save(dfile, "d")
    
end
end
