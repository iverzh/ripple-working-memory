%% Ripple and Unit Firing Analysis During Working Memory Task
% This script analyzes unit firing rates and co-ripple rates in human 
% intracranial recordings during a Sternberg working memory task.
% Analysis includes within and across region comparisons.

close all 
clc
clear

%% Add required paths
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))

%% Define directories and load file lists
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';

% Create directories if they don't exist
if ~isfolder(exportDir); mkdir(exportDir); end
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

% Get subject list
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;

% Load all necessary file lists
unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '../task/*task*'));
taskfiles = {taskfiles.name}';

%% Analysis parameters
recordingState = 'wake';
location = 'NC';
sfreq = 1e3;  % Sampling frequency
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};  % Brain regions of interest

% Time windows for analysis (in ms)
win = [-9 25] * 1e3;  % Window around probe
winResp = [-10.5 1] * 1e3;  % Window around response
winLength = 12905;  % Total trial length in ms
taskMarkers = [500 2600 4700 6800 9400] + 500;  % Task event markers


%% Initialize variables for all subjects
% Preallocate arrays (using estimated max 6500 trials)
maxTrials = 6.5e3;

% Overall ripple and unit activity arrays
recog_trials_all = nan(maxTrials, length(win(1):win(2)));
recog_trials_units_all = nan(maxTrials, length(win(1):win(2)));
recog_trials_resp_all = nan(maxTrials, length(winResp(1):winResp(2)));
recog_trials_units_resp_all = nan(maxTrials, length(winResp(1):winResp(2)));
whole_trial_all = nan(maxTrials, winLength);
whole_trial_units_all = nan(maxTrials, winLength);

% Region-specific arrays
whole_trials_region_all = nan(maxTrials, winLength, length(regions));
resp_region_all = nan(maxTrials, length(winResp(1):winResp(2)), length(regions));
resp_units_region_all = nan(maxTrials, length(winResp(1):winResp(2)), length(regions));
whole_trials_region_all_shuff = nan(maxTrials, winLength, length(regions));
whole_trials_units_sm_region_all = nan(maxTrials, winLength, length(regions));
whole_trials_units_region_all = nan(maxTrials, winLength, length(regions));

% Cross-region arrays
whole_trials_cross_region_all = cell(length(regions));
whole_trials_cross_region_all_shuff = cell(length(regions));
probe_trials_units_all = cell(3, length(regions));  % 3 load conditions
probe_trials_all = cell(3, length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);

% Initialize cross-region matrices
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        whole_trials_cross_region_all{iRa, iRb} = nan(maxTrials, winLength);       
        whole_trials_cross_region_all_shuff{iRa, iRb} = nan(maxTrials, winLength);       
        for ii = 1:size(whole_trials_unit_cross_region_all,3)
            whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(maxTrials, winLength);
        end
    end
end

% Trial metadata
subjID = nan(maxTrials, 1);
correct_trials_all = nan(maxTrials, 1);
load_trials_all = nan(maxTrials, 1);
probe_in_out = nan(maxTrials, 1);
respLatency = nan(maxTrials, 1);


spikePerBundleRip = nan(maxTrials, 1);
nChRegions = zeros(maxTrials, length(regions));

% Counters and statistics collectors
countChPerRegion = zeros(1, length(regions));
countUnitPerRegion = zeros(1, length(regions));
densityAll = [];
freqAll = [];
durAll = [];
ampAll = [];
uChanAll = [];

%% Main processing loop
tic
cTrial = 1;  % Current trial counter
countBundU = 1;  % Bundle unit counter

for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    fprintf('Processing %s ... \n', subject)
    
    % Load micro LFP data
    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    LFPtimeMicr = micrObj.times;
    
    % Process channel labels
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@isempty, locations)) = [];
    
    % Standardize location names
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
    
    % Add hemisphere labels
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locations(hemi) = cellfun(@(X) ['L' X], locations(hemi), 'UniformOutput', false);
    hemi = contains(hem, 'right');
    locations(hemi) = cellfun(@(X) ['R' X], locations(hemi), 'UniformOutput', false);
    
    % Load task data
    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, '../task/',taskfiles{f}));
    
    % Load ripple statistics and event times
    modifier = '1kHz_template_z25';
    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    
    % Update ripple statistics with processed information
    LFPchanNum = cellfun(@str2double, rippleStats.chanLabels);
    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations)
        error('Error formatting channel labels.');
    end
    
    rippleStats.recordingType = repmat({'micro'}, [1 length(microObj.rippleStats.locs)]);
    rippleStats.locs = cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false);
    microObj.rippleStats.window(cellfun(@isempty, microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', ...
                                 microObj.rippleStats.window, 'UniformOutput', false);
    rippleStats.microTimes = LFPtimeMicr;
    
    % Create ripple masks
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask,1)
        if rippleStats.density{chRipp} > 1
            
            
            % Collect ripple statistics
            densityAll = [densityAll rippleStats.density{chRipp}];
            freqAll = [freqAll mean(rippleStats.oscFreq{chRipp})];
            durAll = [durAll mean(rippleStats.duration{chRipp})];
            ampAll = [ampAll mean(rippleStats.rippleAmp{chRipp})];
            
            % Create ripple mask
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
            
            for ii = 1:length(iE)
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
    end
    
    % Calculate overall ripple activity
    rippAll = sum(rippMask);
    
    % Pad arrays to avoid edge effects
    padSize = 3e3;
    rippMask(:, end:end+padSize) = nan;
    rippAll(end:end+padSize) = nan;
    
    % Load and pad LFP data
    data = micrObj.lfp_data;
    data(:, end:end+padSize) = nan;
    
    % Load spike data
    units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    
    % Filter for pyramidal, interneuron, and multi-unit activity
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    
    % Process unit locations
    uLocations = units(U,end-2);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
    
    % Add hemisphere labels to unit locations
    hemi = contains(uLocations, 'left');
    uLocations(hemi) = cellfun(@(X) ['L' X], uLocations(hemi), 'UniformOutput', false);
    hemi = contains(uLocations, 'right');
    uLocations(hemi) = cellfun(@(X) ['R' X], uLocations(hemi), 'UniformOutput', false);
    uLocations = strrep(uLocations, '_left', '');
    uLocations = strrep(uLocations, '_right', '');
    
    % Create spike masks
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001;  % 1 ms bins
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    
    for iU = 1:nNeuron 
        X = units{unitsAll(iU),2};
        [N,~] = histcounts(X, binEdges);
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian', 100);  % 100 ms smoothing
    end
    
    % Calculate overall unit activity
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian', 100);
    else
        unitAll = smoothdata(spikeMask, 'gaussian', 100);
    end
    unitAll = zscore(unitAll);
    
    % Create baseline mask (1 sec before first image)
    baselineMask = false(1, length(rippMask));
    for iT = 1:length(trials.start_time)
        image1Time = round(trials.timestamps_Encoding1(iT)*1e3);
        baselineMask(image1Time-1000:image1Time-1) = true;
    end
    
    % Initialize region-specific data containers
    unit_trials_region = nan(length(trials.start_time), winLength, length(regions));
    dat_region = cell(1,length(regions));
    erp_dat_region = cell(1,length(regions));
    u_dat_region_sm = cell(1,length(regions));
    u_dat_region = cell(1,length(regions));
    resp_per_chan = nan(size(rippMask,1), 140);
    
    % Process each trial
    for iT = 1:length(trials.start_time)
        % Get trial timing information
        respTime = round(trials.timestamps_Response(iT)*1e3);
        probeTime = round(trials.timestamps_Probe(iT)*1e3);
        respLatency(cTrial) = respTime - probeTime;
        
        % Extract overall co-rippling around probe
        trialData = nan(1, win(2)-win(1)+1);
        dat = rippAll(probeTime+win(1):respTime);
        if length(dat) > length(trialData)
            trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;
        end
        recog_trials_all(cTrial, :) = trialData;
        
        % Extract unit activity around response
        dat = sum(spikeMask_sm(:,respTime+winResp(1):respTime+winResp(2)),1);
        recog_trials_units_resp_all(cTrial, :) = dat;
        
        if sum(dat) == 0; continue; end
        
        % Get trial event times
        trialLoad = trials.loads(iT);
        image1Time = round(trials.timestamps_Encoding1(iT)*1e3);
        image2Time = round(trials.timestamps_Encoding2(iT)*1e3);
        image3Time = round(trials.timestamps_Encoding3(iT)*1e3);
        maintenanceTime = round(trials.timestamps_Maintenance(iT)*1e3);
        probeTime = round(trials.timestamps_Probe(iT)*1e3);
        
        % Extract response period data
        dat = rippAll(respTime+winResp(1):respTime+winResp(2));
        recog_trials_resp_all(cTrial, :) = dat;
        
        % Calculate response relative to baseline for each channel
        dat = rippMask(:, respTime-2000:respTime);
        dat = mean(dat,2) ./ mean(rippMask(:, baselineMask),2);
        resp_per_chan(:,iT) = dat;
        
        % Build full trial data based on load condition
        if trialLoad == 3  % High load (3 items)
            % Extract segments for each task phase
            im1seg = rippAll(image1Time-1000:image1Time+2e3);
            im2seg = rippAll(image2Time:image2Time+2e3);
            im3seg = rippAll(image3Time:image3Time+2e3);
            maintenanceSeg = rippAll(maintenanceTime:maintenanceTime+1250);
            probeSeg = rippAll(probeTime-1250:probeTime+3000);
            
            trialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                        im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
            
            % Channel-wise ripple data
            nCh = size(rippMask,1);
            im1seg = rippMask(:, image1Time-1000:image1Time+2e3);
            im2seg = rippMask(:, image2Time:image2Time+2e3);
            im3seg = rippMask(:, image3Time:image3Time+2e3);
            maintenanceSeg = rippMask(:, maintenanceTime:maintenanceTime+1250);
            probeSeg = rippMask(:, probeTime-1250:probeTime+3000);
            
            ripTrialFullChan = [im1seg, nan(nCh,100), im2seg, nan(nCh,100), ...
                               im3seg, nan(nCh,100), maintenanceSeg, nan(nCh,100), probeSeg];
            
            % Unit-wise spike data
            nUnit = size(spikeMask,1);
            im1seg = spikeMask_sm(:, image1Time-1000:image1Time+2e3);
            im2seg = spikeMask_sm(:, image2Time:image2Time+2e3);
            im3seg = spikeMask_sm(:, image3Time:image3Time+2e3);
            maintenanceSeg = spikeMask_sm(:, maintenanceTime:maintenanceTime+1250);
            probeSeg = spikeMask_sm(:, probeTime-1250:probeTime+3000);
            
            unitTrialFullChan = [im1seg, nan(nUnit,100), im2seg, nan(nUnit,100), ...
                                im3seg, nan(nUnit,100), maintenanceSeg, nan(nUnit,100), probeSeg];
            
        elseif trialLoad == 1  % Low load (1 item)
            % Only first image is shown
            im1seg = rippAll(image1Time-1000:image1Time+2e3);
            im2seg = nan(1, 2e3+1);
            im3seg = nan(1, 2e3+1);
            maintenanceSeg = rippAll(maintenanceTime:maintenanceTime+1250);
            probeSeg = rippAll(probeTime-1250:probeTime+3000);
            
            trialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                        im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
            
            % Channel-wise ripple data
            nCh = size(rippMask,1);
            im1seg = rippMask(:, image1Time-1000:image1Time+2e3);
            im2seg = nan(nCh, 2e3+1);
            im3seg = nan(nCh, 2e3+1);
            maintenanceSeg = rippMask(:, maintenanceTime:maintenanceTime+1250);
            probeSeg = rippMask(:, probeTime-1250:probeTime+3000);
            
            ripTrialFullChan = [im1seg, nan(nCh,100), im2seg, nan(nCh,100), ...
                               im3seg, nan(nCh,100), maintenanceSeg, nan(nCh,100), probeSeg];
            
            % Unit-wise spike data
            nUnit = size(spikeMask,1);
            im1seg = spikeMask_sm(:, image1Time-1000:image1Time+2e3);
            im2seg = nan(nUnit, 2e3+1);
            im3seg = nan(nUnit, 2e3+1);
            maintenanceSeg = spikeMask_sm(:, maintenanceTime:maintenanceTime+1250);
            probeSeg = spikeMask_sm(:, probeTime-1250:probeTime+3000);
            
            unitTrialFullChan = [im1seg, nan(nUnit,100), im2seg, nan(nUnit,100), ...
                                im3seg, nan(nUnit,100), maintenanceSeg, nan(nUnit,100), probeSeg];
        end
        
        whole_trial_all(cTrial,:) = trialFull; 
        whole_trial_units_all(cTrial,:) = sum(unitTrialFullChan, 1) / size(unitTrialFullChan,1);
        
        % Process data by brain region
        for iR = 1:length(regions)
            chRegion = contains(rippleStats.chanLabels, regions(iR));
            uRegion = contains(uLocations, regions(iR));
            nChRegions(cTrial, iR) = sum(chRegion);
            
            % Initialize region-specific data on first encounter
            if isempty(dat_region{iR}) && sum(chRegion) > 0
                dat_region{iR} = sum(rippMask(chRegion,:),1);
                erp_dat_region{iR} = mean(data(chRegion,:),1);
                
                % Calculate spike-ripple coupling for each bundle
                bundles = unique(rippleStats.chanLabels(chRegion));
                for iBund = 1:length(bundles)
                    chBundle = contains(rippleStats.chanLabels, bundles(iBund));
                    uBundle = find(contains(uLocations, bundles(iBund)));
                    rippBundleMask = sum(rippMask(chBundle,:)) > 0;
                    rippBnds = mask2bounds(rippBundleMask);
                    
                    for iUnit = 1:length(uBundle)
                        uSpike = spikeMask(uBundle(iUnit),:);
                        spikePerBundleRip(countBundU) = sum(uSpike(rippBundleMask))/size(rippBnds,1);
                        countBundU = countBundU + 1;
                    end
                end
            elseif sum(chRegion) == 0
                dat_region{iR} = nan(1, length(rippMask));
                dat_region_shuff{iR} = nan(1, length(rippMask));
                erp_dat_region{iR} = nan(1, length(rippMask));
            end
            
            % Initialize unit data for region
            if isempty(u_dat_region_sm{iR}) && sum(uRegion) > 0
                u_dat_region_sm{iR} = sum(spikeMask_sm(uRegion,:),1) / sum(uRegion);
                u_dat_region{iR} = sum(spikeMask(uRegion,:),1) / sum(uRegion);
            elseif sum(uRegion) == 0
                u_dat_region_sm{iR} = nan(1, length(rippMask));
                u_dat_region{iR} = nan(1, length(rippMask));
            end
            
            % Build region-specific trial data based on load condition
            if trialLoad == 3
                % Extract segments for each task phase
                im1seg = dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = dat_region{iR}(image2Time:image2Time+2e3);
                im3seg = dat_region{iR}(image3Time:image3Time+2e3);
                maintenanceSeg = dat_region{iR}(maintenanceTime:maintenanceTime+1250);
                probeSeg = dat_region{iR}(probeTime-1250:probeTime+3000);
                
                trialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                            im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                
                
               
                
                % Smoothed unit data
                im1seg = u_dat_region_sm{iR}(image1Time-1000:image1Time+2e3);
                im2seg = u_dat_region_sm{iR}(image2Time:image2Time+2e3);
                im3seg = u_dat_region_sm{iR}(image3Time:image3Time+2e3);
                maintenanceSeg = u_dat_region_sm{iR}(maintenanceTime:maintenanceTime+1250);
                probeSeg = u_dat_region_sm{iR}(probeTime-1250:probeTime+3000);
                
                uTrialFullsm = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                               im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                
                % Raw unit data
                im1seg = u_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = u_dat_region{iR}(image2Time:image2Time+2e3);
                im3seg = u_dat_region{iR}(image3Time:image3Time+2e3);
                maintenanceSeg = u_dat_region{iR}(maintenanceTime:maintenanceTime+1250);
                probeSeg = u_dat_region{iR}(probeTime-1250:probeTime+3000);
                
                uTrialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                             im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                         
            elseif trialLoad == 1
                % Extract segments for each task phase
                im1seg = dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = dat_region{iR}(probeTime-1250:probeTime+3000);
                
                trialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                            im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                
                
               
                
                % Smoothed unit data
                im1seg = u_dat_region_sm{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = u_dat_region_sm{iR}(maintenanceTime:maintenanceTime+1250);
                probeSeg = u_dat_region_sm{iR}(probeTime-1250:probeTime+3000);
                
                uTrialFullsm = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                               im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                
                % Raw unit data
                im1seg = u_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = u_dat_region{iR}(maintenanceTime:maintenanceTime+1250);
                probeSeg = u_dat_region{iR}(probeTime-1250:probeTime+3000);
                
                uTrialFull = [im1seg, nan(1,100), im2seg, nan(1,100), ...
                             im3seg, nan(1,100), maintenanceSeg, nan(1,100), probeSeg];
                
            end
            
            % Store trial data for this region
            unit_trials_region(iT,:,iR) = uTrialFull;           
            whole_trials_region_all(cTrial, :, iR) = trialFull;
            resp_region_all(cTrial, :, iR) = dat_region{iR}(respTime+winResp(1):respTime+winResp(2));
            resp_units_region_all(cTrial, :, iR) = u_dat_region_sm{iR}(respTime+winResp(1):respTime+winResp(2));
            whole_trials_units_sm_region_all(cTrial, :, iR) = uTrialFullsm;
            whole_trials_units_region_all(cTrial, :, iR) = uTrialFull;
            
            % Calculate unit and ripple rates for baseline, maintenance, and probe periods
            uRegion = find(uRegion);  % Convert to indices
            for ii = 1:length(uRegion)
                baslnUnits = mean(spikeMask(uRegion(ii),image1Time-1e3:image1Time));
                maintUnits = mean(spikeMask(uRegion(ii),maintenanceTime:maintenanceTime+2e3));
                probeUnits = mean(spikeMask(uRegion(ii),probeTime:probeTime+1250));
                probe_trials_units_all{trialLoad,iR} = [probe_trials_units_all{trialLoad,iR}; ...
                    [countUnitPerRegion(iR)+ii baslnUnits maintUnits probeUnits]];
            end
            
            chRegion = find(chRegion);  % Convert to indices
            for ii = 1:length(chRegion)
                baslnRipples = mean(rippMask(chRegion(ii),image1Time-1e3:image1Time));
                maintRipples = mean(rippMask(chRegion(ii),maintenanceTime:maintenanceTime+2e3));
                probeRipples = mean(rippMask(chRegion(ii),probeTime:probeTime+1250));
                probe_trials_all{trialLoad,iR} = [probe_trials_all{trialLoad,iR}; ...
                    [countChPerRegion(iR)+ii baslnRipples maintRipples probeRipples]];
            end
        end
        
        % Calculate cross-region co-rippling
        for iRa = 1:length(regions)
            for iRb = iRa:length(regions)
                datArip = whole_trials_region_all(cTrial,:,iRa);
                datBrip = whole_trials_region_all(cTrial,:,iRb);
                
                datAripShuff = whole_trials_region_all_shuff(cTrial,:,iRa);
                datBripShuff = whole_trials_region_all_shuff(cTrial,:,iRb);
                
                if iRa == iRb
                    % Within-region case
                    datAB = datArip;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                    whole_trials_cross_region_all_shuff{iRa, iRb}(cTrial, :) = datAripShuff;
                    whole_trials_unit_cross_region_all{iRa, iRb, 3}(cTrial, :) = unit_trials_region(iT,:,iRa);
                else
                    % Cross-region case: only count when both regions have ripples
                    datAB = datArip + datBrip;
                    datAB(datArip <= 0 | datBrip <= 0) = 0;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                    
                    datAB = datAripShuff + datBripShuff;
                    datAB(datAripShuff <= 0 | datBripShuff <= 0) = 0;
                    whole_trials_cross_region_all_shuff{iRa, iRb}(cTrial, :) = datAB;
                end
                
                % Find overlapping units and channels between regions
                uRegionA = uChan(contains(uLocations, regions(iRa)));
                uRegionB = uChan(contains(uLocations, regions(iRb)));
                
                ovA = nan(1, length(uRegionA));
                ovAcell = arrayfun(@(X) find(LFPchanNum == X), uRegionA', 'UniformOutput', false);
                ovA(cellfun(@(X) ~isempty(X), ovAcell)) = cell2mat(ovAcell);
                ovB = nan(1, length(uRegionB));
                ovBcell = arrayfun(@(X) find(LFPchanNum == X), uRegionB', 'UniformOutput', false);
                ovB(cellfun(@(X) ~isempty(X), ovBcell)) = cell2mat(ovBcell);
                
                if sum(~isnan(ovA)) < 1 || sum(~isnan(ovB)) < 1; continue; end
            end
        end
        
        % Extract unit activity around probe
        trialData = nan(1, win(2)-win(1)+1);
        dat = unitAll(probeTime+win(1):respTime);
        if length(dat) > length(trialData)
            trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;
        end
        recog_trials_units_all(cTrial,:) = trialData;
        
        % Store trial metadata
        correct_trials_all(cTrial) = trials.response_accuracy(iT);
        load_trials_all(cTrial) = trials.loads(iT);
        probe_in_out(cTrial) = trials.probe_in_out(iT);
        subjID(cTrial) = subj;
        
        cTrial = cTrial + 1;
    end
    
    % Update counters for channels and units per region
    for iR = 1:length(regions)
        chRegion = contains(rippleStats.chanLabels, regions(iR));
        countChPerRegion(iR) = countChPerRegion(iR) + sum(chRegion);
        
        uRegion = contains(uLocations, regions(iR));
        countUnitPerRegion(iR) = countUnitPerRegion(iR) + sum(uRegion);
    end
    
    fprintf('Number of total trials: %i, number of subject trials: %i \n\n', cTrial - 1, length(trials.start_time))
end

% Trim unused rows
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) = [];
probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('Done processing.\n')
toc

%% Analysis: Ripple spike firing  during task vs baseline
close all

% Calculate spike rates per ripple duration for task vs baseline
baseSpikePerBund = [];
baseRipDurPerBund = [];
taskSpikePerBund = [];
taskRipDurPerBund = [];

for iS = 1:max(subjID)
    subjectTrials = find(subjID == iS);
    
    for iR = 1:length(regions)
        baseSpike = 0;
        baseRipDur = 0;
        taskSpike = 0;
        taskRipDur = 0;
        
        for iT = 1:length(subjectTrials)
            ripTrial = whole_trials_region_all(subjectTrials(iT), :, iR);
            uTrial = whole_trials_unit_cross_region_all{iR, iR, 3}(subjectTrials(iT), :);
            
            if sum(uTrial, 'omitnan') == 0; continue; end
            
            % Baseline period (first 1 sec)
            baselineRipp = ripTrial(1:1e3) > 0;
            baselineUnit = uTrial(1:1e3);
            baseSpike = baseSpike + sum(baselineUnit(baselineRipp), 'omitnan');
            baseRipDur = baseRipDur + sum(baselineRipp) / 1e3;
            
            % Task period (after baseline)
            taskRipp = ripTrial(1e3+1:end) > 0;
            taskUnit = uTrial(1e3+1:end);
            taskSpike = taskSpike + sum(taskUnit(taskRipp), 'omitnan');
            taskRipDur = taskRipDur + sum(taskRipp) / 1e3;
        end
        
        if taskSpike == 0 && baseSpike == 0; continue; end
        
        baseSpikePerBund = [baseSpikePerBund baseSpike];
        baseRipDurPerBund = [baseRipDurPerBund baseRipDur];
        taskSpikePerBund = [taskSpikePerBund taskSpike];
        taskRipDurPerBund = [taskRipDurPerBund taskRipDur];
    end
end

% Test if spike rate during ripples increases during task vs baseline
[P, H, STATS] = signrank((taskSpikePerBund ./ taskRipDurPerBund), (baseSpikePerBund ./ baseRipDurPerBund), ...
                         'tail', 'right');
fprintf('Spike-ripple coupling task vs baseline: p = %.4f\n', P);

%% Statistical analysis: Compare ripple and firing rates across task phases
regionColors = brewermap(12, 'Dark2');
plot_condition = load_trials_all == 3;  % High load trials

% Initialize arrays for ripple rates (coR) and firing rates (FR)
coR = nan(length(plot_condition), length(regions)+1, 6);  % 6 task phases
FR = nan(length(plot_condition), length(regions)+1, 6);
pLoadCoR = nan(length(regions)+1, 6);  % p-values for load effect
pLoadFR = nan(length(regions)+1, 6);
pBslnCoR = nan(length(regions)+1, 6, 2);  % p-values vs baseline
pBslnFR = nan(length(regions)+1, 6, 2);

% Calculate mean rates for each task phase
for iR = 1:length(regions)+1
    for iT = 1:length(plot_condition)
        if iR == 6  % Overall (all regions combined)
            if plot_condition(iT)  % High load
                coR(iT, iR, 1) = std(whole_trial_all(iT, 1:taskMarkers(1)), 'omitnan');  % Baseline
                coR(iT, iR, 2) = mean(whole_trial_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');  % Encoding 1
                coR(iT, iR, 3) = mean(whole_trial_all(iT, taskMarkers(2):taskMarkers(3)), 'omitnan');  % Encoding 2
                coR(iT, iR, 4) = mean(whole_trial_all(iT, taskMarkers(3):taskMarkers(4)), 'omitnan');  % Encoding 3
                coR(iT, iR, 5) = mean(whole_trial_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');  % Maintenance
                coR(iT, iR, 6) = mean(whole_trial_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');  % Probe
                
                FR(iT, iR, 1) = mean(whole_trial_units_all(iT, 1:taskMarkers(1)), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trial_units_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                FR(iT, iR, 3) = mean(whole_trial_units_all(iT, taskMarkers(2):taskMarkers(3)), 'omitnan');
                FR(iT, iR, 4) = mean(whole_trial_units_all(iT, taskMarkers(3):taskMarkers(4)), 'omitnan');
                FR(iT, iR, 5) = mean(whole_trial_units_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trial_units_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');
            else  % Low load
                coR(iT, iR, 1) = std(whole_trial_all(iT, 1:taskMarkers(1)), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trial_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                coR(iT, iR, 3:4) = nan;  % No encoding 2 or 3
                coR(iT, iR, 5) = mean(whole_trial_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trial_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');
                
                FR(iT, iR, 1) = mean(whole_trial_units_all(iT, 1:taskMarkers(1)), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trial_units_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                FR(iT, iR, 3:4) = nan;
                FR(iT, iR, 5) = mean(whole_trial_units_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trial_units_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');
            end
        else  % Individual regions
            if plot_condition(iT)  % High load
                coR(iT, iR, 1) = mean(whole_trials_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trials_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                coR(iT, iR, 3) = mean(whole_trials_region_all(iT, taskMarkers(2):taskMarkers(3), iR), 'omitnan');
                coR(iT, iR, 4) = mean(whole_trials_region_all(iT, taskMarkers(3):taskMarkers(4), iR), 'omitnan');
                coR(iT, iR, 5) = mean(whole_trials_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trials_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');
                
                FR(iT, iR, 1) = mean(whole_trials_units_sm_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                FR(iT, iR, 3) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(2):taskMarkers(3), iR), 'omitnan');
                FR(iT, iR, 4) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(3):taskMarkers(4), iR), 'omitnan');
                FR(iT, iR, 5) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');
            else  % Low load
                coR(iT, iR, 1) = mean(whole_trials_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trials_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                coR(iT, iR, 3:4) = nan;
                coR(iT, iR, 5) = mean(whole_trials_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trials_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');
                
                FR(iT, iR, 1) = mean(whole_trials_units_sm_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                FR(iT, iR, 3:4) = nan;
                FR(iT, iR, 5) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');
            end
        end
    end
    
    % Statistical testing for each task phase
    for iTm = [1,2,5,6]  % Test baseline, encoding 1, maintenance, probe
        % Test load effect (high vs low) using linear mixed effects
        dat = coR(:, iR, iTm);
        tab = table(subjID(~isnan(dat)), plot_condition(~isnan(dat)), dat(~isnan(dat)), ...
                   'VariableNames', {'subject', 'condition', 'rates'});
        lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
        pLoadCoR(iR, iTm) = lme.Coefficients(2,6);
        
        dat = FR(:, iR, iTm);
        tab = table(subjID(~isnan(dat)), plot_condition(~isnan(dat)), dat(~isnan(dat)), ...
                   'VariableNames', {'subject', 'condition', 'rates'});
        lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
        pLoadFR(iR, iTm) = lme.Coefficients(2,6);
        
        % Test vs baseline using Wilcoxon rank sum
        pBslnCoR(iR, iTm, 2) = ranksum(coR(plot_condition, iR, iTm), coR(:, iR, 1));  % High load vs baseline
        pBslnCoR(iR, iTm, 1) = ranksum(coR(~plot_condition, iR, iTm), coR(:, iR, 1));  % Low load vs baseline
        
        pBslnFR(iR, iTm, 2) = ranksum(FR(plot_condition, iR, iTm), FR(:, iR, 1));
        pBslnFR(iR, iTm, 1) = ranksum(FR(~plot_condition, iR, iTm), FR(:, iR, 1));
        
        % Also test encoding 2 and 3 for high load condition
        if iTm == 2
            pBslnCoR(iR, iTm+1, 2) = ranksum(coR(plot_condition, iR, iTm+1), coR(:, iR, 1));
            pBslnCoR(iR, iTm+2, 2) = ranksum(coR(plot_condition, iR, iTm+2), coR(:, iR, 1));
            
            pBslnFR(iR, iTm+1, 2) = ranksum(FR(plot_condition, iR, iTm+1), FR(:, iR, 1));
            pBslnFR(iR, iTm+2, 2) = ranksum(FR(plot_condition, iR, iTm+2), FR(:, iR, 1));
        end
    end    
end

% Convert firing rates to Hz
FR = FR * 1e3;

% Apply FDR correction for multiple comparisons
pLoadCoR_adj = nan(size(pLoadCoR));
[~, ~, ~, pLoadCoR_adj(~isnan(pLoadCoR))] = fdr_bh(pLoadCoR(~isnan(pLoadCoR)), 0.05, 'pdep', 'yes');

pLoadFR_adj = nan(size(pLoadFR));
[~, ~, ~, pLoadFR_adj(~isnan(pLoadFR))] = fdr_bh(pLoadFR(~isnan(pLoadFR)), 0.05, 'pdep', 'yes');

pBslnCoR_adj = nan(size(pBslnCoR));
[~, ~, ~, pBslnCoR_adj(~isnan(pBslnCoR))] = fdr_bh(pBslnCoR(~isnan(pBslnCoR)), 0.05, 'pdep', 'yes');

pBslnFR_adj = nan(size(pBslnFR));
[~, ~, ~, pBslnFR_adj(~isnan(pBslnFR))] = fdr_bh(pBslnFR(~isnan(pBslnFR)), 0.05, 'pdep', 'yes');

%% Figure: Co-ripple rates across task phases by region
figure('Position', [669 805 285 round(729*(6/5))]);

for iR = 1:length(regions)+1
    if iR < 6
        subplot(6,1,iR+1)
        cBox = [0 0 0; regionColors(iR,:)];
    else
        subplot(6,1,1)
        cBox = [0 0 0; 1 0 0];
    end
    
    load1 = squeeze(coR(~plot_condition, iR, :))';
    load3 = squeeze(coR(plot_condition, iR, :))';
    
    group_num = [ones(1, length(load1)), 2*ones(1, length(load3))];
    dat = [load1, load3]';
    h = daboxplot(dat, 'groups', group_num, 'mean', 1, 'outlier', 0, 'color', cBox, ...
                  'xshift', 0, 'boxalpha', 0.2, 'boxwidth', 1.3, 'meansize', 4, 'mediansize', 2);
    hold on;
    
    ax = gca;
    ax.XTick = 1:6;
    ax.XTickLabel = {'B','E1','E2','E3','M','P'};  % Baseline, Encoding 1-3, Maintenance, Probe
    
    % Add baseline reference line
    hl = hline(median(coR(:,iR,1), 'all', 'omitnan'));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    % Add significance markers for baseline comparisons
    pVals = find(pBslnCoR_adj(iR,:,1) < 0.05);  % Low load
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(h.gpos(1, iTm), y_max, repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnCoR_adj(iR,iTm,1), 'maxNum',1)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
    
    pVals = find(pBslnCoR_adj(iR,:,2) < 0.05);  % High load
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(mean(h.gpos(2, iTm)), y_max, repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnCoR_adj(iR,iTm,2), 'maxNum',1)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
    
    % Add significance markers for load comparisons
    pVals = find(pLoadCoR_adj(iR,:) < 0.05);
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.3/1.27 * y_max;
        line([h.gpos(1, iTm) h.gpos(2, iTm)], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);
        text(mean(h.gpos(:, iTm)), y_max + 0.05/1.27 * y_max, ...
             repmat('*', 1, countLeadingZerosAfterDecimal(pLoadCoR_adj(iR,iTm), 'maxNum',3)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegionsBoxPlot_%s.pdf', tag)))

%% Figure: Firing rates across task phases by region
figure('Position', [669 805 285 round(729*(6/5))]);

for iR = 1:length(regions)+1
    if iR < 6
        subplot(6,1,iR+1)
        cBox = [0 0 0; regionColors(iR,:)];
    else
        subplot(6,1,1)
        cBox = [0 0 0; 1 0 0];
    end
    
    load1 = squeeze(FR(~plot_condition, iR, :))';
    load3 = squeeze(FR(plot_condition, iR, :))';
    
    group_num = [ones(1, length(load1)), 2*ones(1, length(load3))];
    dat = [load1, load3]';
    h = daboxplot(dat, 'groups', group_num, 'mean', 1, 'outlier', 0, 'color', cBox, ...
                  'xshift', 0, 'boxalpha', 0.2, 'boxwidth', 1.3, 'meansize', 4, 'mediansize', 2);
    hold on;
    
    ax = gca;
    ax.XTick = 1:6;
    ax.XTickLabel = {'B','E1','E2','E3','M','P'};
    
    % Add baseline reference line
    hl = hline(median(FR(:,iR,1), 'all', 'omitnan'));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    % Add significance markers for baseline comparisons
    pVals = find(pBslnFR_adj(iR,:,1) < 0.05);  % Low load
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(h.gpos(1, iTm), y_max, repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnFR_adj(iR,iTm,1), 'maxNum',1)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
    
    pVals = find(pBslnFR_adj(iR,:,2) < 0.05);  % High load
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(mean(h.gpos(2, iTm)), y_max, repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnFR_adj(iR,iTm,2), 'maxNum',1)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
    
    % Add significance markers for load comparisons
    pVals = find(pLoadFR_adj(iR,:) < 0.05);
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]) + 0.30;
        y_max = y_max + 0.3/1.27 * y_max;
        line([h.gpos(1, iTm) h.gpos(2, iTm)], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5);
        text(mean(h.gpos(:, iTm)), y_max + 0.05/1.27 * y_max, ...
             repmat('*', 1, countLeadingZerosAfterDecimal(pLoadFR_adj(iR,iTm), 'maxNum',3)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
    end
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegionsBoxPlot_%s.pdf', tag)))



close all

%% Initialize visualization parameters
taskMarkers = [500 2600 4700 6800 9400] + 500;  % Task event markers
respLatency = respLatency(~isnan(respLatency));  % Clean response latencies
sm = 100;  % Smoothing window
win = [-3e3 500];  % Analysis window around response
regionColors = brewermap(12, 'Dark2');  % Color scheme for regions

%% Figure: Response-locked ripple and unit activity comparison
% Compare probe-in vs probe-out trials
figure('Position', [904 414 round(285/2) 522]);

% Define comparison conditions
plot_condition = probe_in_out == 1;    % Probe-in trials
plot_condition_2 = probe_in_out == 0;  % Probe-out trials

% Clean infinite values
recog_trials_resp_all(isinf(recog_trials_resp_all)) = nan;

% Subplot 1: Response-locked ripple activity
subplot(3,1,1)
plot_data = smoothdata(recog_trials_resp_all(plot_condition_2,:), 2, 'gaussian', sm);
[bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                       std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2)), ...
                       'k', 'nan', 'gap');
hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(recog_trials_resp_all(plot_condition,:), 2, 'gaussian', sm);
[bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                       std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), ...
                       'r', 'nan', 'gap');
bf.FaceAlpha = 0.7;
vline(0)
xlim([win(1) win(2)])
ylabel('Ripple rate')
title('Response-locked activity')

% Subplot 2: Response time distributions
subplot(3,1,2)
histogram(-respLatency(plot_condition), win(1):100:win(2), 'Normalization', 'probability');
hold on;
histogram(-respLatency(plot_condition_2), win(1):100:win(2), 'Normalization', 'probability')
ylabel('Probability')
xlabel('Time from response (ms)')

% Subplot 3: Response-locked unit activity
subplot(3,1,3)
plot_data = smoothdata(recog_trials_units_resp_all(plot_condition,:), 2, 'gaussian', sm);
[bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                       std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), ...
                       'r', 'nan', 'gap');
hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(recog_trials_units_resp_all(plot_condition_2,:), 2, 'gaussian', sm);
[bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                       std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2)), ...
                       'k', 'nan', 'gap');
bf.FaceAlpha = 0.7;
vline(0)
xlim([win(1) win(2)])
ylabel('Unit firing rate')
xlabel('Time from response (ms)')

%% Figure 2: Response time histogram
figure('Position', [893 434 175 196]);
h1 = histogram(-respLatency(plot_condition), win(1):100:win(2), 'Normalization', 'probability');
h1.FaceColor = 'r';
hold on;
h2 = histogram(-respLatency(plot_condition_2), win(1):100:win(2), 'Normalization', 'probability');
h2.FaceColor = 'k';
fig = gcf;
fig.Color = 'w';
ylabel('Proportion')
xlabel('Time from key press (ms)')
savepdf(gcf, fullfile(exportDirFigs, sprintf('respTime_load_%s.pdf', tag)))


%% Figure: Task-evoked co-ripple activity 
figure('Position', [669 805 285 round(729/5)]);

% Compare load 1 vs load 2
plot_condition = load_trials_all == 3;
plot_condition_2 = load_trials_all == 1;

% Find example trials for NaN masking
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');

% Plot slow response trials
plot_data = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian', 100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:))) = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap');
hold on;
bf.FaceAlpha = 0.7;

% Plot fast response trials
plot_data = smoothdata(whole_trial_all(plot_condition_2,:), 2, 'gaussian', 100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2));
plot_data_mu(isnan(whole_trial_all(condTrial2,:))) = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;

[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap');
bf.FaceAlpha = 0.7;

% Add task markers and baseline reference
vl = vline([500 2600 4700 6800 9400]-500);
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_all(:,1:1e3), 'all', 'omitnan'));
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-1e3 max(taskMarkers) + 2e3])
fig = gcf;
fig.Color = 'w';
ylabel('Co-ripples')
xlabel('Time from trial start (ms)')
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTask_allChannels_load_%s.pdf', tag)))

%% Figure: Task-evoked unit firing (fast vs slow responses)
figure('Position', [669 805 285 round(729/5)]);

% Load 3 trials
plot_data = smoothdata(whole_trial_units_all(plot_condition,:), 2, 'gaussian', 100) * 1e3;
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_units_all(condTrial1,:))) = nan;
plot_data_sem(isnan(whole_trial_units_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap');
hold on;
bf.FaceAlpha = 0.7;

% Load 1 trials
plot_data = smoothdata(whole_trial_units_all(~plot_condition,:), 2, 'gaussian', 100) * 1e3;
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_units_all(condTrial2,:))) = nan;
plot_data_sem(isnan(whole_trial_units_all(condTrial2,:))) = nan;

[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap');
bf.FaceAlpha = 0.7;

% Add task markers and baseline
vl = vline([500 2600 4700 6800 9400]-500);
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_units_all(:,1:1e3), 'all', 'omitnan')*1e3);
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-1e3 max(taskMarkers) + 2e3])
fig = gcf;
fig.Color = 'w';
ylabel('Unit firing (Hz)')
xlabel('Time from trial start (ms)')
savepdf(gcf, fullfile(exportDirFigs, sprintf('Firing_allChannels_load_%s.pdf', tag)))

%% Figure: Probe in/out comparison
figure('Position', [669 805 285 round(729/5)]);

plot_condition = probe_in_out == 1;  % Probe-in trials
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');

% Plot probe-in trials
plot_data = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian', 200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:))) = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap');
hold on;
bf.FaceAlpha = 0.7;

% Plot probe-out trials
plot_data = smoothdata(whole_trial_all(~plot_condition,:), 2, 'gaussian', 200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial2,:))) = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;

[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap');
bl2.Color = [247, 133, 177]./255;
bf.FaceColor = bl2.Color;
bf.FaceAlpha = 0.25;

% Add markers
vl = vline([500 2600 4700 6800 9400]-500);
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_all(:,1:1e3), 'all', 'omitnan'));
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-1e3 max(taskMarkers) + 2e3])
fig = gcf;
fig.Color = 'w';
ylabel('Co-ripples')
xlabel('Time from trial start (ms)')
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTask_allChannels_probe_%s.pdf', tag)))



%% Figure: Regional ripple dynamics (load comparison)
figure('Position', [669 805 285 729]);

plot_condition = load_trials_all == 3;
plot_condition_2 = load_trials_all == 1;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');

for iR = 1:length(regions)
    subplot(5,1,iR)
    
    % Fast response trials
    plot_data_1 = whole_trials_region_all(plot_condition_2,:,iR);
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian', 100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian', 100);
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:))) = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    
    [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap');
    hold on;
    bf.FaceAlpha = 0.7;
    
    % Slow response trials
    plot_data_3 = whole_trials_region_all(plot_condition,:,iR);
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian', 100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian', 100);
    plot_data_mu_3(isnan(whole_trial_all(condTrial1,:))) = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;
    
    [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap');
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    % Add baseline reference
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    % Statistical testing: Response period
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2);
    probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2);
    probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'}, 'paired', 'off', ...
                            'method', 'perm', 'naccu', 10000, 'verbose', 'off');
    fprintf('%s ripples response p = %.2f\n', regions{iR}, p_perm)
    
    % Statistical testing: Maintenance period
    probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2);
    probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2);
    probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'}, 'paired', 'off', ...
                            'method', 'perm', 'naccu', 10000, 'verbose', 'off');
    fprintf('%s ripples maintenance p = %.2f\n', regions{iR}, p_perm)
    
    % Format axes
    xlim([-1e3 max(taskMarkers) + 2e3])
    vl = vline([500 2600 4700 6800 9400]-500);
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
    ylabel('Co-ripples')
    if iR == length(regions); xlabel('Time from trial start (ms)'); end
    
    ax = gca;
    ax.FontSize = 10;
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegions_%s.pdf', tag)))
%% Figure: Response-locked regional ripple activity
figure('Position', [669 805 round(285/2) 729]);

plot_condition = load_trials_all == 3;

for iR = 1:length(regions)
    subplot(5,1,iR)
    
    % Calculate baseline for normalization (per subject)
    baseline = nan(6500,1);
    for iS = 1:max(subjID)
        ii = find(subjID == iS);
        baseline(ii) = mean(whole_trials_region_all(ii,1:1e3,iR), 'all', 'omitnan');
    end
    
    % Low load trials
    plot_data_1 = resp_region_all(~plot_condition,:,iR);
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian', 100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian', 100);
    
    [bl2, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap');
    hold on;
    bf.FaceAlpha = 0.7;
    
    % High load trials
    plot_data_3 = resp_region_all(plot_condition,:,iR);
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian', sm);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian', sm);
    
    [bl1, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap');
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    % Add baseline reference
    hl = hline(mean([plot_data_mu_1(end-1e3:end), plot_data_mu_3(end-1e3:end)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    xlim([-3.0e3 500])
    vl = vline(0);
    vl.LineWidth = 1.0;
    
    if iR == length(regions); xlabel('Time from key press (ms)'); end
    
    ax = gca;
    ax.FontSize = 9;
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegionsResp_%s.pdf', tag)))
%% Figure: Response-locked regional unit activity
figure('Position', [669 805 round(285/2) 729]);

for iR = 1:length(regions)
    subplot(5,1,iR)
    
    % Calculate baseline per subject
    baseline = nan(6500,1);
    for iS = 1:max(subjID)
        ii = find(subjID == iS);
        baseline(ii) = mean(whole_trials_units_sm_region_all(ii,1:1e3,iR), 'all', 'omitnan');
    end
    
    % Fast response trials
    plot_data_1 = resp_units_region_all(~plot_condition,:,iR) * 1e3;  % Convert to Hz
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    
    [bl2, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap');
    hold on;
    bf.FaceAlpha = 0.7;
    
    % Slow response trials
    plot_data_3 = resp_units_region_all(plot_condition,:,iR) * 1e3;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    
    [bl1, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap');
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    % Add baseline reference
    hl = hline(mean([plot_data_mu_1(end-1e3:end), plot_data_mu_3(end-1e3:end)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    xlim([-3e3 500])
    vl = vline(0);
    vl.LineWidth = 1.0;
    
    if iR == length(regions); xlabel('Time from key press (ms)'); end
    
    ax = gca;
    ax.FontSize = 9;
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegionsResp_%s.pdf', tag)))
%% Figure: Regional unit firing dynamics
figure('Position', [669 805 285 729]);

condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');

for iR = 1:length(regions)
    subplot(5,1,iR)
    
    % Fast response trials
    plot_data_1 = whole_trials_units_sm_region_all(plot_condition_2,:,iR) * 1e3;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:))) = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    
    [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap');
    hold on;
    bf.FaceAlpha = 0.7;
    
    % Slow response trials
    plot_data_3 = whole_trials_units_sm_region_all(plot_condition,:,iR) * 1e3;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3(isnan(whole_trial_all(condTrial1,:))) = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;
    
    [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap');
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    % Add baseline reference
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 1.0;
    
    % Statistical testing for probe period
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2);
    probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2);
    probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'}, 'paired', 'off', ...
                            'method', 'perm', 'naccu', 10000, 'verbose', 'off');
    fprintf('%s firing p = %.2f\n', regions{iR}, p_perm)
    
    % Format axes
    xlim([-1e3 max(taskMarkers) + 2e3])
    vl = vline([500 2600 4700 6800 9400]-500);
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
    ylabel('Unit firing (Hz)')
    if iR == length(regions); xlabel('Time from trial start (ms)'); end
    
    ax = gca;
    ax.FontSize = 10;
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegions_%s.pdf', tag)))

%% Figure: Cross-regional co-rippling analysis
figure('Position', [1281 840 772 733]);

% Initialize storage for statistics
tabCoRip = cell(5);  % Tables for linear mixed effects models
lmeCoRip = cell(5);  % LME model objects
lmePval = nan(5);    % P-values from LME

for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        % Calculate subplot position
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)
        
        % Low load trials
        plot_data_1 = whole_trials_cross_region_all{iRa,iRb}(~plot_condition,:);
        plot_data_mu_1 = mean(plot_data_1, 'omitnan');
        plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~plot_condition));
        plot_data_mu_1 = smoothdata(plot_data_mu_1, 'gaussian', 100);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian', 100);
        plot_data_mu_1(isnan(whole_trial_all(condTrial2,:))) = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
        
        [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap');
        hold on;
        bf.FaceAlpha = 0.7;
        bl2.LineWidth = 0.5;
        
        % High load trials
        plot_data_3 = whole_trials_cross_region_all{iRa,iRb}(plot_condition,:);
        plot_data_mu_3 = mean(plot_data_3, 'omitnan');
        plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(plot_condition));
        plot_data_mu_3 = smoothdata(plot_data_mu_3, 'gaussian', 100);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian', 100);
        plot_data_mu_3(isnan(whole_trial_all(condTrial1,:))) = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;
        
        [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap');
        hold on;
        bf.FaceAlpha = 0.7;
        bl1.Color = mean([regionColors(iRa,:); regionColors(iRb,:)]);
        bf.FaceColor = mean([regionColors(iRa,:); regionColors(iRb,:)]);
        bl1.MarkerFaceColor = bl1.Color;
        bf.FaceAlpha = 0.3;
        bl1.LineWidth = 0.5;
        
        % Add task markers
        vl = vline([500 2600 4700 6800 9400]-500);
        for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
        
        % Add baseline reference
        hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
        hl.Color = 'k';
        hl.LineStyle = '-';
        hl.LineWidth = 0.5;
        
        % Calculate co-ripple rates for different task periods
        coRresp = mean(whole_trials_cross_region_all{iRa,iRb}(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2, 'omitnan');
        coRmaintenance = mean(whole_trials_cross_region_all{iRa,iRb}(:,taskMarkers(end-1):taskMarkers(end)), 2, 'omitnan');
        coRbaseline = mean(whole_trials_cross_region_all{iRa,iRb}(:,1:taskMarkers(1)), 2, 'omitnan');
        
        % Remove excess trials
        coRresp(cTrial:end) = [];
        coRmaintenance(cTrial:end) = [];
        coRbaseline(cTrial:end) = [];
        
        % Create table for LME analysis
        tabCoRip{iRa,iRb} = table(subjID(~isnan(coRresp)), ...
                                  plot_condition(~isnan(coRresp)), ...
                                  coRresp(~isnan(coRresp)), ...
                                  coRmaintenance(~isnan(coRresp)), ...
                                  coRbaseline(~isnan(coRresp)), ...
                                  'VariableNames', {'subject', 'condition', 'rates_resp', ...
                                                   'rates_maintenance', 'rates_baseline'});
        
        % Fit linear mixed effects model for maintenance period
        lmeCoRip{iRa,iRb} = fitlme(tabCoRip{iRa,iRb}, 'rates_maintenance ~ condition + (1|subject)');
        lmePval(iRa,iRb) = lmeCoRip{iRa,iRb}.Coefficients(2,6);
        
        % Calculate and report percent change
        pct = (mean(coRmaintenance(plot_condition), 'omitnan') - ...
               mean(coRmaintenance(~plot_condition), 'omitnan')) / ...
               mean(coRmaintenance(~plot_condition), 'omitnan');
        fprintf('%s <--> %s ripples maintenance %.2f%%\n', regions{iRa}, regions{iRb}, pct*100)
        
        % Format subplot
        if iRa == 1; ylabel('Co-ripples'); end
        if iRb == length(regions); xlabel('Time (ms)'); end
        
        plotMarkers = [500 2600 4700 6800 9400]-500;
        xlim([plotMarkers(end-1)-200 plotMarkers(end)+2.0e3])
    end
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskAcrossRegions_%s.pdf', tag)))

%% Figure: Statistical significance matrix for cross-regional comparisons
% Apply FDR correction to p-values
adj_p = nan(size(lmePval));
[~, ~, ~, adj_p(~isnan(lmePval))] = fdr_bh(lmePval(~isnan(lmePval)), 0.05, 'pdep', 'yes');

figure('Position', [1738 1327 246 182]);
imagesc(adj_p', [min(adj_p(:))-1e-2 1]);
hold on;
colorbar;

% Add grid lines
hl = hline((1/2):size(adj_p,1)+(1/2), 'k-');
vl = vline((1/2):size(adj_p,1)+(1/2), 'k-');

% Create custom colormap (white to magma)
cmap = flipud(slanCM('magma'));
cmap = [1 1 1; cmap];
colormap(cmap)

% Mark significant comparisons
[yy, xx] = find(adj_p' < 0.05);
plot(xx, yy, 'r*');

% Format axes
ax = gca;
ax.YTick = 1:size(adj_p,1);
ax.XTick = 1:size(adj_p,1);
ax.XTickLabel = regions;
ax.YTickLabel = regions;
ax.Title.String = 'Cross-regional co-rippling (p-values)';

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskAcrossRegions_stats_maintenance_%s.pdf', tag)))
%% Figure: Summary bar plot of significant cross-regional effects
figure('Position', [1209 1255 254 294]);

% Create bar plot with counts of significant effects
% Note: These values are  calculated from the actual statistical results
% from oscillations in all gamma bands
br = bar([5 12; 0 1; 2 7]');  
hold on;
br(1).FaceColor = [230, 159, 0]/255;  % Orange
br(2).FaceColor = [0, 158, 115]/255;  % Green
br(3).FaceColor = [213, 94, 0]/255;   % Red
ylim([0 15]);

ax = gca;
ax.XTickLabel = {'Maintenance', 'Probe'};
ax.YLabel.String = 'Number of significant pairs';

legend('Co-ripples', 'Co-LG', 'Co-vHG', 'location', 'northwest')
box off

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('AcrossRegionsBar_stats_%s.pdf', tag)))

















