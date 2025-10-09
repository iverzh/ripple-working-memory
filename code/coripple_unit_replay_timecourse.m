

%% Replay Timecourse Analysis During Working Memory Task
% This script analyzes the temporal relationship between neural replay events
% and hippocampal ripples during encoding and probe periods of a Sternberg task.
% Analysis focuses on co-firing patterns within ripple windows.

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
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';

% Create directories if needed
if ~isfolder(exportDir); mkdir(exportDir); end
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

% Get subject list
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;

% Load all file lists
unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory,'../task', '*task*'));
taskfiles = {taskfiles.name}';

%% Analysis parameters
recordingState = 'wake';
location = 'NC';
modifier = '1kHz_template_z25';
tag = [recordingState,'_',location,'_',modifier];

% Define regions with hemispheres
regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};  % For visualization

% Task timing parameters
taskMarkers = [500 2600 4700 6800 9400] + 500;  % Task event markers in ms
winLength = 12905;  % Total trial length in ms
win = [-9 25] * 1e3;  % Analysis window around events
winResp = [-1.5 1] * 1e3;  % Response window

% Co-firing window parameter (in ms)
coFireWindow = round(25/2);  % 25ms total window = ?12.5ms

% Load channel curation information
channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end,3:end-1));
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2,3:end-1));
bundleNotes = table2cell(bundleNotes(3:end,3:end-1));

%% Initialize analysis variables
maxTrials = 6.5e3;  % Maximum expected trials

% Trial-level data storage
recog_trials_all = nan(maxTrials, length(win(1):win(2)));
recog_trials_units_all = nan(maxTrials, length(win(1):win(2)));
recog_trials_resp_all = nan(maxTrials, length(winResp(1):winResp(2)));
whole_trial_all = nan(maxTrials, winLength);
whole_trials_region_all = nan(maxTrials, winLength, length(regions));
whole_trials_units_region_all = nan(maxTrials, winLength, length(regions));

% Cross-region storage
whole_trials_cross_region_all = cell(length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);

% Replay analysis storage - co-ripple periods
coR_replay_PRTHe = nan(5e5, length(regions), length(regions));  % Encoding PSTH
coR_replay_PRTHp = nan(5e5, length(regions), length(regions));  % Probe PSTH
coR_replay_PRTHe_shuff = nan(5e5, length(regions), length(regions));  % Shuffled encoding
coR_replay_PRTHp_shuff = nan(5e5, length(regions), length(regions));  % Shuffled probe

% Replay analysis storage - single region ripples
rA_replay_PRTHe = nan(5e5, length(regions), length(regions));  % Region A ripples encoding
rA_replay_PRTHp = nan(5e5, length(regions), length(regions));  % Region A ripples probe

% Counters for normalization
countPRTHp = ones(length(regions), length(regions));  % Probe PSTH counter
countPRTHe = ones(length(regions), length(regions));  % Encoding PSTH counter
countCoRp = ones(length(regions), length(regions));   % Co-ripple probe counter
countCoRe = ones(length(regions), length(regions));   % Co-ripple encoding counter

countPRTHp_shuff = ones(length(regions), length(regions));
countPRTHe_shuff = ones(length(regions), length(regions));
countCoRp_shuff = ones(length(regions), length(regions));
countCoRe_shuff = ones(length(regions), length(regions));

countPRTHpA = ones(length(regions), length(regions));
countPRTHeA = ones(length(regions), length(regions));
countRAe = ones(length(regions), length(regions));
countRAp = ones(length(regions), length(regions));

% Trial metadata
subjID = nan(maxTrials, 1);
correct_trials_all = nan(maxTrials, 1);
load_trials_all = nan(maxTrials, 1);
probe_in_out = nan(maxTrials, 1);
respLatency = nan(maxTrials, 1);

% Statistics collectors
densityAll = [];
uChanAll = [];

%% Main processing loop
tic
cTrial = 1;  % Trial counter

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
    trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
    
    % Load ripple statistics
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    LFPchanNum = cellfun(@str2double, rippleStats.chanLabels);
    
    rippleStats.chanLabels = locations;
    if length(rippleStats.locs) ~= length(locations)
        error('Error formatting channel labels.');
    end
    
    % Update ripple statistics
    rippleStats.recordingType = repmat({'micro'}, [1 length(microObj.rippleStats.locs)]);
    rippleStats.locs = cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false);
    microObj.rippleStats.window(cellfun(@isempty, microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', ...
                                 microObj.rippleStats.window, 'UniformOutput', false);
    rippleStats.microTimes = LFPtimeMicr;
    
    % Create ripple mask
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask,1)
        if rippleStats.density{chRipp} > 1
            % Skip macro channels and SPE channels
            if strcmp(rippleStats.recordingType{chRipp}, 'macro'); continue; end
            if contains(rippleStats.chanLabels{chRipp}, 'SPE')
                rippMask(chRipp,:) = nan(1,length(rippMask));
                continue
            end
            
            % Collect ripple statistics
            densityAll = [densityAll rippleStats.density{chRipp}];
            
            % Create ripple mask for this channel
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
            
            for ii = 1:length(iE)
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
    end
    
    % Pad arrays to avoid edge effects
    padSize = 3e3;
    rippMask(:, end:end+padSize) = nan;
    micrObj.lfp_data(:, end:end+padSize) = nan;
    
    rippAll = sum(rippMask);
    
    % Load spike data
    units = LoadSpikeTimes(subject, 'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    
    % Filter for valid unit types
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
    
    % Add hemisphere labels
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
        [N, ~] = histcounts(X, binEdges);
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian', 500);
    end
    
    % Calculate overall unit activity
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian', 500);
    else
        unitAll = smoothdata(spikeMask, 'gaussian', 500);
    end
    unitAll = zscore(unitAll);
    
    %% Process cross-regional replay during co-ripples
    for iRa = 1:length(regions)
        for iRb = iRa+1:length(regions)
            
            % Get unit indices for each region
            uRegionA = find(contains(uLocations, regions(iRa)));
            uRegionB = find(contains(uLocations, regions(iRb)));
            
            % Get ripple masks for each region
            datArip = sum(rippMask(contains(locations, regions(iRa)),:), 'omitnan') > 0;
            datBrip = sum(rippMask(contains(locations, regions(iRb)),:), 'omitnan') > 0;
            
            if sum(datArip) == 0 || sum(datBrip) == 0; continue; end
            
            % Process each trial
            for iT = 1:length(trials.start_time)
                
                % Skip probe-out trials for main analysis
                if ~trials.probe_in_out(iT); continue; end
                
                trialLoad = trials.loads(iT);
                
                % Get trial timing
                imageTime(1) = round(trials.timestamps_Encoding1(iT)*1e3);
                imageTime(2) = round(trials.timestamps_Encoding2(iT)*1e3);
                imageTime(3) = round(trials.timestamps_Encoding3(iT)*1e3);
                imageTime(4) = round(trials.timestamps_Maintenance(iT)*1e3);
                probeTime = round(trials.timestamps_Probe(iT)*1e3);
                respTime = round(trials.timestamps_Response(iT)*1e3);
                
                % Get image IDs
                probeIm = trials.PicIDs_Probe(iT);
                encodIm(1) = trials.PicIDs_Encoding1(iT);
                encodIm(2) = trials.PicIDs_Encoding2(iT);
                encodIm(3) = trials.PicIDs_Encoding3(iT);
                
                % Process unit pairs between regions
                for uA = 1:length(uRegionA)
                    for uB = 1:length(uRegionB)
                        % Get spike data for probe period (1 sec after probe)
                        spikeMaskP = spikeMask(:, probeTime+1:probeTime+1e3);
                        
                        % Identify co-ripple and single-ripple periods
                        datABcoripP = datArip(probeTime+1:probeTime+1e3) > 0 & ...
                                     datBrip(probeTime+1:probeTime+1e3) > 0;
                        datAripP = datArip(probeTime+1:probeTime+1e3) > 0 & ...
                                  datBrip(probeTime+1:probeTime+1e3) == 0;
                        
                        % Get encoding period data based on trial load
                        if trialLoad == 3
                            % Find which encoding image matches probe
                            iE = find(encodIm(:) == probeIm);
                            if isempty(iE); continue; end
                            
                            % Get spike and ripple data for encoding period (2 sec after image)
                            spikeMaskE = spikeMask(:, imageTime(iE):imageTime(iE)+2e3);
                            datABcoripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & ...
                                         datBrip(imageTime(iE):imageTime(iE)+2e3) > 0;
                            datAripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & ...
                                      datBrip(imageTime(iE):imageTime(iE)+2e3) == 0;
                        else
                            % For load 1, only use first image
                            spikeMaskE = spikeMask(:, imageTime(1):imageTime(1)+2e3);
                            datABcoripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & ...
                                         datBrip(imageTime(1):imageTime(1)+2e3) > 0;
                            datAripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & ...
                                      datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                        end
                        
                        %% Calculate co-firing during co-ripples
                        % Get spike times for each unit
                        datAP = find(spikeMaskP(uRegionA(uA),:));
                        datBP = find(spikeMaskP(uRegionB(uB),:));
                        datAE = find(spikeMaskE(uRegionA(uA),:));
                        datBE = find(spikeMaskE(uRegionB(uB),:));
                        
                        if isempty(datAP) || isempty(datBP) || isempty(datAE) || isempty(datBE)
                            continue;
                        end
                        
                        % Calculate co-firing within time window for probe
                        sP = datAP - datBP';
                        [coFpB, coFpA] = find(sP >= -coFireWindow*2 & sP <= coFireWindow*2);
                        
                        % Calculate co-firing within time window for encoding
                        sE = datAE - datBE';
                        [coFeB, coFeA] = find(sE >= -coFireWindow*2 & sE <= coFireWindow*2);
                        
                        if ~isempty(coFeA) && ~isempty(coFpA)
                            % Get co-ripple boundaries
                            bCoRipE = mean(mask2bounds(datABcoripE), 2)';
                            bCoRipP = mean(mask2bounds(datABcoripP), 2)';
                            bRipAE = mean(mask2bounds(datAripE), 2)';
                            bRipAP = mean(mask2bounds(datAripP), 2)';
                            
                            % Calculate mean co-firing times
                            mat = nan(2, length(coFpA));
                            mat(1,:) = datAP(coFpA);
                            mat(2,:) = datBP(coFpB);
                            coFp = mean(mat, 1);
                            
                            mat = nan(2, length(coFeA));
                            mat(1,:) = datAE(coFeA);
                            mat(2,:) = datBE(coFeB);
                            coFe = mean(mat, 1);
                            
                            % Calculate PSTH: co-firing time relative to co-ripple center
                            cE = bCoRipE - coFe';  % Encoding
                            cP = bCoRipP - coFp';  % Probe
                            
                            % Store co-ripple PSTH data
                            coR_replay_PRTHe(countPRTHe(iRa,iRb):countPRTHe(iRa,iRb)+length(cE(:))-1, iRa, iRb) = cE(:);
                            countPRTHe(iRa,iRb) = countPRTHe(iRa,iRb) + length(cE(:));
                            countCoRe(iRa,iRb) = countCoRe(iRa,iRb) + length(bCoRipE);
                            
                            coR_replay_PRTHp(countPRTHp(iRa,iRb):countPRTHp(iRa,iRb)+length(cP(:))-1, iRa, iRb) = cP(:);
                            countPRTHp(iRa,iRb) = countPRTHp(iRa,iRb) + length(cP(:));
                            countCoRp(iRa,iRb) = countCoRp(iRa,iRb) + length(bCoRipP);
                            
                            % Calculate PSTH for single region A ripples
                            cE = bRipAE - coFe';
                            cP = bRipAP - coFp';
                            
                            rA_replay_PRTHe(countPRTHeA(iRa,iRb):countPRTHeA(iRa,iRb)+length(cE(:))-1, iRa, iRb) = cE(:);
                            countPRTHeA(iRa,iRb) = countPRTHeA(iRa,iRb) + length(cE(:));
                            countRAe(iRa,iRb) = countRAe(iRa,iRb) + length(bRipAE);
                            
                            rA_replay_PRTHp(countPRTHpA(iRa,iRb):countPRTHpA(iRa,iRb)+length(cP(:))-1, iRa, iRb) = cP(:);
                            countPRTHpA(iRa,iRb) = countPRTHpA(iRa,iRb) + length(cP(:));
                            countRAp(iRa,iRb) = countRAp(iRa,iRb) + length(bRipAP);
                        end
                        
                        %% Generate shuffled control
                        % Randomize spike times within window
                        datAP_shuff = randi(length(spikeMaskP), [1, length(datAP)]);
                        datBP_shuff = randi(length(spikeMaskP), [1, length(datBP)]);
                        datAE_shuff = randi(length(spikeMaskE), [1, length(datAE)]);
                        datBE_shuff = randi(length(spikeMaskE), [1, length(datBE)]);
                        
                        % Calculate shuffled co-firing
                        sP = datAP_shuff - datBP_shuff';
                        [coFpB, coFpA] = find(sP >= -coFireWindow*2 & sP <= coFireWindow*2);
                        
                        sE = datAE_shuff - datBE_shuff';
                        [coFeB, coFeA] = find(sE >= -coFireWindow*2 & sE <= coFireWindow*2);
                        
                        if ~isempty(coFeA) && ~isempty(coFpA)
                            % Get co-ripple boundaries (same as real data)
                            bCoRipE = mean(mask2bounds(datABcoripE), 2)';
                            bCoRipP = mean(mask2bounds(datABcoripP), 2)';
                            
                            % Calculate mean shuffled co-firing times
                            mat = nan(2, length(coFpA));
                            mat(1,:) = datAP_shuff(coFpA);
                            mat(2,:) = datBP_shuff(coFpB);
                            coFp = mean(mat, 1);
                            
                            mat = nan(2, length(coFeA));
                            mat(1,:) = datAE_shuff(coFeA);
                            mat(2,:) = datBE_shuff(coFeB);
                            coFe = mean(mat, 1);
                            
                            % Calculate shuffled PSTH
                            cE = bCoRipE - coFe';
                            cP = bCoRipP - coFp';
                            
                            % Store shuffled data
                            coR_replay_PRTHe_shuff(countPRTHe_shuff(iRa,iRb):countPRTHe_shuff(iRa,iRb)+length(cE(:))-1, iRa, iRb) = cE(:);
                            countPRTHe_shuff(iRa,iRb) = countPRTHe_shuff(iRa,iRb) + length(cE(:));
                            countCoRe_shuff(iRa,iRb) = countCoRe_shuff(iRa,iRb) + length(bCoRipE);
                            
                            coR_replay_PRTHp_shuff(countPRTHp_shuff(iRa,iRb):countPRTHp_shuff(iRa,iRb)+length(cP(:))-1, iRa, iRb) = cP(:);
                            countPRTHp_shuff(iRa,iRb) = countPRTHp_shuff(iRa,iRb) + length(cP(:));
                            countCoRp_shuff(iRa,iRb) = countCoRp_shuff(iRa,iRb) + length(bCoRipP);
                        end
                    end % uB
                end % uA
            end % trials
        end % iRb
    end % iRa
end % subjects

% Trim unused rows
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) = [];
probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('Done processing.\n')
toc
%% Visualization: PSTH of replay relative to ripple centers
regionColors = brewermap(12, 'Dark2');
contraColor = brewermap(12, 'Accent');
close all

% Load region location data if available
load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')

% PSTH parameters
psthWindow = 130;  % Window for PSTH calculation (ms)
plotWindow = 100;  % Window for plotting (ms)
smoothWindow = 11;  % Smoothing window
binSize = 1;       % Bin size (ms)

% Initialize figures for encoding and probe PSTHs
figure('Position', [1324 421 364 136]);  % Encoding figure
figure('Position', [1324 421 364 136]);  % Probe figure

%% Calculate and plot PSTH for ipsilateral pairs
iRapl = find(contains(regions, regionsPlot(1:5)));
iRbpl = find(contains(regions, regionsPlot(1:5)));

% Initialize storage for averaged PSTHs
Np = [];  % Probe co-ripple
Ne = [];  % Encoding co-ripple
NpA = []; % Probe single ripple
NeA = []; % Encoding single ripple
NpShuff = [];  % Probe shuffled
NeShuff = [];  % Encoding shuffled

% Process ipsilateral pairs (same hemisphere)
for Aloop = iRapl
    for Bloop = iRbpl
        iiA = min([Aloop Bloop]);
        iiB = max([Aloop Bloop]);
        
        % Skip if no data or different hemispheres
        if countPRTHp(iiA,iiB) == 1; continue; end
        if ~strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
        
        % Calculate encoding PSTH
        times = coR_replay_PRTHe_shuff(:, iiA, iiB);
        tempShuff = smoothdata(histcounts(times(~isnan(times)), ...
                              -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRe_shuff(iiA,iiB), ...
                              'gaussian', smoothWindow);
        
        % Z-score normalize together
        z = zscore([temp, tempA, tempShuff]);
        Ne = [Ne; z(1:length(z)/3)];
        NeA = [NeA; z((length(z)/3)+1:(2*length(z)/3))];
        NeShuff = [NeShuff; z((2*length(z)/3)+1:end)];
        
        % Calculate probe PSTH
        times = coR_replay_PRTHp(:, iiA, iiB);
        temp = smoothdata(histcounts(times(~isnan(times)), ...
                         -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRp(iiA,iiB), ...
                         'gaussian', smoothWindow);
        
        times = rA_replay_PRTHp(:, iiA, iiB);
        tempA = smoothdata(histcounts(times(~isnan(times)), ...
                          -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countRAp(iiA,iiB), ...
                          'gaussian', smoothWindow);
        
        times = coR_replay_PRTHp_shuff(:, iiA, iiB);
        tempShuff = smoothdata(histcounts(times(~isnan(times)), ...
                              -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRp_shuff(iiA,iiB), ...
                              'gaussian', smoothWindow);
        
        % Z-score normalize together
        z = zscore([temp, tempA, tempShuff]);
        Np = [Np; z(1:length(z)/3)];
        NpA = [NpA; z((length(z)/3)+1:(2*length(z)/3))];
        NpShuff = [NpShuff; z((2*length(z)/3)+1:end)];
    end
end

%% Plot ipsilateral encoding PSTH
figure(1)
subplot(1,2,1)

% Add shaded region for co-firing window
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255);
hold on;
f.LineStyle = '--';
f.LineWidth = 1;

% Plot co-ripple PSTH
[bl1, bf1] = boundedline(-psthWindow:binSize:psthWindow, mean(Ne), ...
                        std(Ne)/sqrt(size(Ne,1)), 'b');
bl1.Color = regionColors(6,:);
bl1.LineWidth = 1.5;
bf1.FaceColor = regionColors(6,:);
bf1.FaceAlpha = 0.3;

% Plot single ripple PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NeA), ...
                        std(NeA)/sqrt(size(NeA,1)));
bl3.Color = 0.8*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

% Plot shuffled PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NeShuff), ...
                        std(NeShuff)/sqrt(size(NeShuff,1)));
bl3.Color = 0.5*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

xlim([-plotWindow plotWindow])
vl = vline(0);
vl.LineWidth = 1.5;
box off
xlabel('Time from ripple center (ms)')
ylabel('Z-scored co-firing rate')
title('Encoding - Ipsilateral')

% Statistical testing
times = -psthWindow:binSize:psthWindow;
xx = [];
xx(:,1) = mean(Ne(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NeA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NeShuff(:, times >= -25 & times <= 25), 2);
p = anova1(xx, [], 'off');
fprintf('Encoding ipsilateral ANOVA p: %.4f\n', p)

ax = gca;
ax.LineWidth = 1;

%% Plot ipsilateral probe PSTH
figure(2)
subplot(1,2,1)

% Add shaded region
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255);
hold on;
f.LineStyle = '--';
f.LineWidth = 1;

% Plot co-ripple PSTH
[bl2, bf2] = boundedline(-psthWindow:binSize:psthWindow, mean(Np), ...
                        std(Np)/sqrt(size(Np,1)), 'r');
bl2.Color = regionColors(6,:);
bl2.LineWidth = 1.5;
bf2.FaceColor = regionColors(6,:);
bf2.FaceAlpha = 0.3;

% Plot single ripple PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NpA), ...
                        std(NpA)/sqrt(size(NpA,1)));
bl3.Color = 0.8*bl2.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

% Plot shuffled PSTH
[bl4, bf4] = boundedline(-psthWindow:binSize:psthWindow, mean(NpShuff), ...
                        std(NpShuff)/sqrt(size(NpShuff,1)));
bl4.Color = 0.5*bl2.Color;
bl4.LineWidth = 1.5;
bf4.FaceColor = bl3.Color;
bf4.FaceAlpha = 0.3;

xlim([-plotWindow plotWindow])
vl = vline(0);
vl.LineWidth = 1.5;
box off
xlabel('Time from ripple center (ms)')
ylabel('Z-scored co-firing rate')
title('Probe - Ipsilateral')

% Statistical testing
xx = [];
xx(:,1) = mean(Np(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NpA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NpShuff(:, times >= -25 & times <= 25), 2);
p = anova1(xx, [], 'off');
fprintf('Probe ipsilateral ANOVA p: %.4f\n', p)

ax = gca;
ax.LineWidth = 1;

%% Calculate and plot PSTH for contralateral pairs
% Reset storage
Np = [];
Ne = [];
NpA = [];
NeA = [];
NpShuff = [];
NeShuff = [];

% Process contralateral pairs (different hemispheres)
for Aloop = iRapl
    for Bloop = iRbpl
        iiA = min([Aloop Bloop]);
        iiB = max([Aloop Bloop]);
        
        % Skip if no data or same hemisphere
        if countPRTHp(iiA,iiB) == 1; continue; end
        if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
        
        % Calculate encoding PSTH
        times = coR_replay_PRTHe(:, iiA, iiB);
        temp = smoothdata(histcounts(times(~isnan(times)), ...
                         -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRe(iiA,iiB), ...
                         'gaussian', smoothWindow);
        
        times = rA_replay_PRTHe(:, iiA, iiB);
        tempA = smoothdata(histcounts(times(~isnan(times)), ...
                          -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countRAe(iiA,iiB), ...
                          'gaussian', smoothWindow);
        
        times = coR_replay_PRTHe_shuff(:, iiA, iiB);
        tempShuff = smoothdata(histcounts(times(~isnan(times)), ...
                              -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRe_shuff(iiA,iiB), ...
                              'gaussian', smoothWindow);
        
        % Z-score normalize
        z = zscore([temp, tempA, tempShuff]);
        Ne = [Ne; z(1:length(z)/3)];
        NeA = [NeA; z((length(z)/3)+1:(2*length(z)/3))];
        NeShuff = [NeShuff; z((2*length(z)/3)+1:end)];
        
        % Calculate probe PSTH
        times = coR_replay_PRTHp(:, iiA, iiB);
        temp = smoothdata(histcounts(times(~isnan(times)), ...
                         -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRp(iiA,iiB), ...
                         'gaussian', smoothWindow);
        
        times = rA_replay_PRTHp(:, iiA, iiB);
        tempA = smoothdata(histcounts(times(~isnan(times)), ...
                          -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countRAp(iiA,iiB), ...
                          'gaussian', smoothWindow);
        
        times = coR_replay_PRTHp_shuff(:, iiA, iiB);
        tempShuff = smoothdata(histcounts(times(~isnan(times)), ...
                              -psthWindow-(binSize/2):binSize:psthWindow+(binSize/2)) / countCoRp_shuff(iiA,iiB), ...
                              'gaussian', smoothWindow);
        
        % Z-score normalize
        z = zscore([temp, tempA, tempShuff]);
        Np = [Np; z(1:length(z)/3)];
        NpA = [NpA; z((length(z)/3)+1:(2*length(z)/3))];
        NpShuff = [NpShuff; z((2*length(z)/3)+1:end)];
    end
end

fprintf('\n')

%% Plot contralateral encoding PSTH
figure(1)
subplot(1,2,2)

% Add shaded region
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255);
hold on;
f.LineStyle = '--';
f.LineWidth = 1;

% Plot co-ripple PSTH
[bl1, bf1] = boundedline(-psthWindow:binSize:psthWindow, mean(Ne), ...
                        std(Ne)/sqrt(size(Ne,1)), 'b');
bl1.Color = contraColor(7,:);
bl1.LineWidth = 1.5;
bf1.FaceColor = contraColor(7,:);
bf1.FaceAlpha = 0.3;

% Plot single ripple PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NeA), ...
                        std(NeA)/sqrt(size(NeA,1)));
bl3.Color = 0.8*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

% Plot shuffled PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NeShuff), ...
                        std(NeShuff)/sqrt(size(NeShuff,1)));
bl3.Color = 0.5*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

xlim([-plotWindow plotWindow])
vl = vline(0);
vl.LineWidth = 1.5;
box off
xlabel('Time from ripple center (ms)')
ylabel('Z-scored co-firing rate')
title('Encoding - Contralateral')

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplayPRTH_%s_E.pdf', tag)))

% Statistical testing
times = -psthWindow:binSize:psthWindow;
xx = [];
xx(:,1) = mean(Ne(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NeA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NeShuff(:, times >= -25 & times <= 25), 2);
p = anova1(xx, [], 'off');
fprintf('Encoding contralateral ANOVA p: %.4f\n', p)

ax = gca;
ax.LineWidth = 1;

%% Plot contralateral probe PSTH
figure(2)
subplot(1,2,2)

% Add shaded region
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255);
hold on;
f.LineStyle = '--';
f.LineWidth = 1;

% Plot co-ripple PSTH
[bl2, bf2] = boundedline(-psthWindow:binSize:psthWindow, mean(Np), ...
                        std(Np)/sqrt(size(Np,1)), 'r');
bl2.Color = contraColor(7,:);
bl2.LineWidth = 1.5;
bf2.FaceColor = contraColor(7,:);
bf2.FaceAlpha = 0.3;

% Plot single ripple PSTH
[bl3, bf3] = boundedline(-psthWindow:binSize:psthWindow, mean(NpA), ...
                        std(NpA)/sqrt(size(NpA,1)));
bl3.Color = 0.8*bl2.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

% Plot shuffled PSTH
[bl4, bf4] = boundedline(-psthWindow:binSize:psthWindow, mean(NpShuff), ...
                        std(NpShuff)/sqrt(size(NpShuff,1)));
bl4.Color = 0.5*bl2.Color;
bl4.LineWidth = 1.5;
bf4.FaceColor = bl3.Color;
bf4.FaceAlpha = 0.3;

xlim([-plotWindow plotWindow])
vl = vline(0);
vl.LineWidth = 1.5;
box off
xlabel('Time from ripple center (ms)')
ylabel('Z-scored co-firing rate')
title('Probe - Contralateral')

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplayPRTH_%s_P.pdf', tag)))

% Statistical testing
xx = [];
xx(:,1) = mean(Np(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NpA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NpShuff(:, times >= -25 & times <= 25), 2);
p = anova1(xx, [], 'off');
fprintf('Probe contralateral ANOVA p: %.4f\n', p)
















