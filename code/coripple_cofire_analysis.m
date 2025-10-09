close all 
clc
clear

%% ========================================================================
%  CO-RIPPLE CO-FIRING ANALYSIS FOR WORKING MEMORY TASK
%  ========================================================================
%  This script analyzes co-firing  between brain regions during
%  ripple events in a Sternberg working memory task. It computes co-firing
%  rates during ripples, non-ripples, and overall, comparing across task
%  phases and memory loads.
%


%% Add Required Paths
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))

%% Define Directories
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';

% Create output directories if they don't exist
if ~isfolder(exportDir); mkdir(exportDir); end
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

%% Get Subject and File Lists
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '../task/*task*'));
taskfiles = {taskfiles.name}';

%% Analysis Parameters
recordingState = 'wake';
location = 'NC';
sfreq = 1e3;  % Sampling frequency in Hz

% Brain regions of interest
regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};

% Co-firing window (ms) - spikes within ?coFireWindow are considered co-firing
coFireWindow = 25;
coFireWindow = round(coFireWindow/2);

% Time windows for analysis (ms)
win = [-9 25] * 1e3;           % Full trial window
winResp = [-1.5 1] * 1e3;      % Response window
winLength = 12905;              % Total trial length (ms)

% Subsampling factor for co-firing analysis
coFireSubsample = 1;

% Task timing markers (ms from trial start)
taskMarkers = [500 2600 4700 6800 9400] + 500;

%% Preallocate Data Structures
% Trial-level data
nMaxTrials = 6.5e3;
subjID = nan(nMaxTrials, 1);
correct_trials_all = nan(nMaxTrials, 1);
load_trials_all = nan(nMaxTrials, 1);
probe_in_out = nan(nMaxTrials, 1);
respLatency = nan(nMaxTrials, 1);

recog_trials_all = nan(nMaxTrials, length(win(1):win(2)));
recog_trials_units_all = nan(nMaxTrials, length(win(1):win(2)));
recog_trials_resp_all = nan(nMaxTrials, length(winResp(1):winResp(2)));

whole_trial_all = nan(nMaxTrials, winLength);
whole_trials_region_all = nan(nMaxTrials, winLength, length(regions));
whole_trials_units_region_all = nan(nMaxTrials, winLength, length(regions));

% Region pair data structures
nMaxPairs = 4.5e5;
nMaxTrialsPair = 3.5e5;
nRegions = length(regions);

whole_trials_cross_region_all = cell(nRegions);
whole_trials_unit_cross_region_all = cell(nRegions, nRegions, 9);
whole_trials_cross_region_dur = cell(nRegions, nRegions, 9);
load_trials_all_unitPair = cell(nRegions);
unitPairID = cell(nRegions);
testCoFire = cell(nRegions);
UnitPairType = cell(nRegions);
UnitPairFR = cell(nRegions);

trialsProbe = zeros(nRegions, nRegions, 1e3);

% Initialize cell arrays for each region pair
for iRa = 1:nRegions
    for iRb = iRa:nRegions
        whole_trials_cross_region_all{iRa, iRb} = nan(nMaxTrials, winLength);       
        testCoFire{iRa, iRb} = zeros(1, ceil(winLength/coFireSubsample));       
        load_trials_all_unitPair{iRa, iRb} = nan(1, nMaxTrialsPair);    
        unitPairID{iRa, iRb} = nan(nMaxPairs, 3);
        UnitPairType{iRa, iRb} = cell(2, nMaxPairs);
        UnitPairFR{iRa, iRb} = nan(2, nMaxPairs);
        
        for ii = 1:9
            whole_trials_unit_cross_region_all{iRa, iRb, ii} = nan(4.5e4, 4);
            whole_trials_cross_region_dur{iRa, iRb, ii} = nan(4.5e4, 4);
        end
    end
end

% Ripple statistics
rippleFR = nan(length(subj_list_full), 100);
CoRippleRatesRegionAll = nan(nRegions, nRegions, length(subj_list_full));
densityAll = [];
uChanAll = [];

%%  MAIN ANALYSIS LOOP - Process Each Subject
tic
cTrial = 1;
cUnitPair = ones(nRegions);

% Ripple detection modifier
modifier = '1kHz_template_z25';
tag = [recordingState, '_', location, '_', modifier];

for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    fprintf('Processing %s...\n', subject)
    
    %% Load LFP Data
    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    LFPtimeMicr = micrObj.times;
    
    % Parse channel labels and assign anatomical locations
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@(X) isempty(X), locations)) = [];
    
    % Standardize region names
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
    
    % Add hemisphere prefix
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locations(hemi) = cellfun(@(X) ['L' X], locations(hemi), 'UniformOutput', false);
    hemi = contains(hem, 'right');
    locations(hemi) = cellfun(@(X) ['R' X], locations(hemi), 'UniformOutput', false);
    
    %% Load Task Data
    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
    
    %% Load Ripple Detection Results
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    LFPchanNum = cellfun(@(X) str2double(X), rippleStats.chanLabels);
    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations)
        error('Error formatting channel labels.')
    end
    
    % Convert ripple timing to microsecond times
    rippleStats.recordingType = repmat({'micro'}, [1 length(microObj.rippleStats.locs)]);
    rippleStats.locs = cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false);
    microObj.rippleStats.window(cellfun(@(X) isempty(X), microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', ...
                                  microObj.rippleStats.window, 'UniformOutput', false);
    rippleStats.microTimes = LFPtimeMicr;
    
    %% Create Ripple Mask
    % Binary mask indicating ripple periods for each channel
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask, 1)
        if rippleStats.density{chRipp} <= 1; continue; end
        if ~strcmp(rippleStats.recordingType{chRipp}, 'micro'); continue; end
        if contains(rippleStats.chanLabels{chRipp}, 'SPE')
            rippMask(chRipp, :) = nan(1, length(rippMask));
            continue
        end
        
        times = rippleStats.microTimes;
        densityAll = [densityAll rippleStats.density{chRipp}];
        
        % Mark ripple periods in mask
        iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
        iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
        
        for ii = 1:length(iE)
            if any([iS(ii) iE(ii)] <= 0); continue; end
            rippMask(chRipp, iS(ii):iE(ii)) = 1;
        end
    end
    
    % Pad the ending to avoid edge effects
    rippMask(:, end:end+3e3) = nan;
    micrObj.lfp_data(:, end:end+3e3) = nan;
    
    rippAll = sum(rippMask);
    
    %% Load and Process Unit Data
    units = LoadSpikeTimes(subject, 'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    
    % Filter for pyramidal and interneuron units
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    
    % Assign anatomical locations to units
    uLocations = units(U, end-2);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
    
    % Add hemisphere prefix
    hemi = contains(uLocations, 'left');
    uLocations(hemi) = cellfun(@(X) ['L' X], uLocations(hemi), 'UniformOutput', false);
    hemi = contains(uLocations, 'right');
    uLocations(hemi) = cellfun(@(X) ['R' X], uLocations(hemi), 'UniformOutput', false);
    uLocations = strrep(uLocations, '_left', '');
    uLocations = strrep(uLocations, '_right', '');
    
    %% Create Spike Masks and Compute Firing Rates
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001; % 1ms bins
    
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    FR = zeros(nNeuron, 1);
    
    for iU = 1:nNeuron
        X = units{unitsAll(iU), 2};
        [N, EDGES] = histcounts(X, binEdges);
        FR(iU) = sum(N) / rippleStats.recordingLength * 1000;
        spikeMask(iU, :) = N;
        spikeMask_sm(iU, :) = smoothdata(N, 'gaussian', 500);
    end
    
    % Compute population activity
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian', 500);
    else
        unitAll = smoothdata(spikeMask, 'gaussian', 500);
    end
    unitAll = zscore(unitAll);
    
    %% Create Task Phase Masks
    % Create binary masks for different task phases
    taskMask = false(4, length(rippMask));
    taskMaskLoad3 = false(4, length(rippMask));
    taskMaskLoad1 = false(4, length(rippMask));
    
    for iT = 1:length(trials.start_time)
        trialLoad = trials.loads(iT);
        
        % Get timing of task events
        image1Time = round(trials.timestamps_Encoding1(iT) * 1e3);
        maintenanceTime = round(trials.timestamps_Maintenance(iT) * 1e3);
        probeTime = round(trials.timestamps_Probe(iT) * 1e3);
        
        % Define task phase masks:
        % 1: Baseline (1s before first image)
        % 2: Encoding (first image presentation)
        % 3: Maintenance (delay period)
        % 4: Probe (probe presentation)
        taskMask(1, image1Time-1000:image1Time-1) = true;
        taskMask(2, image1Time:image1Time+1e3) = true;
        taskMask(3, maintenanceTime:probeTime) = true;
        taskMask(4, probeTime+1:probeTime+1e3) = true;
        
        % Separate masks by memory load
        if trialLoad == 3
            taskMaskLoad3(1, image1Time-1000:image1Time-1) = true;
            taskMaskLoad3(2, image1Time:image1Time+1e3) = true;
            taskMaskLoad3(3, maintenanceTime:probeTime) = true;
            taskMaskLoad3(4, probeTime+1:probeTime+1e3) = true;
        else
            taskMaskLoad1(1, image1Time-1000:image1Time-1) = true;
            taskMaskLoad1(2, image1Time:image1Time+1e3) = true;
            taskMaskLoad1(3, maintenanceTime:probeTime) = true;
            taskMaskLoad1(4, probeTime+1:probeTime+1e3) = true;
        end
    end
    
    %  Cross-Region Co-Firing Analysis
    % Analyze co-firing between all region pairs
    for iRa = 1:nRegions
        for iRb = iRa:nRegions
            % Get units in each region
            uLFPregionA = uChan(contains(uLocations, regions(iRa)));
            uLFPregionB = uChan(contains(uLocations, regions(iRb)));
            uRegionA = find(contains(uLocations, regions(iRa)));
            uRegionB = find(contains(uLocations, regions(iRb)));
            
            % Get ripple activity for each region
            datArip = sum(rippMask(contains(locations, regions(iRa)), :), 'omitnan') > 0;
            datBrip = sum(rippMask(contains(locations, regions(iRb)), :), 'omitnan') > 0;
            
            if sum(datArip) == 0 || sum(datBrip) == 0; continue; end
            
            % Analyze all unit pairs between regions
            for uA = 1:length(uRegionA)
                for uB = 1:length(uRegionB)
                    if uLFPregionA(uA) == uLFPregionB(uB); continue; end
                    
                    typeA = units{uRegionA(uA), 3};
                    typeB = units{uRegionB(uB), 3};
                    
                    % Analyze each task phase
                    for ii = 1:4
                        % Define co-ripple and no-ripple periods
                        datABcorip = datArip > 0 & datBrip > 0 & taskMask(ii, :);
                        datABnorip = datArip == 0 & datBrip == 0 & taskMask(ii, :);
                        
                        %% Co-Firing During Co-Ripples
                        datA = find(spikeMask(uRegionA(uA), datABcorip));
                        datB = find(spikeMask(uRegionB(uB), datABcorip));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 1}(cUnitPair(iRa, iRb), ii) = sum(datABcorip);
                        
                        % Load 3 specific
                        datABcoripLoad3 = datABcorip & taskMaskLoad3(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABcoripLoad3));
                        datB = find(spikeMask(uRegionB(uB), datABcoripLoad3));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 4}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 4}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 3);
                        
                        % Load 1 specific
                        datABcoripLoad1 = datABcorip & taskMaskLoad1(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABcoripLoad1));
                        datB = find(spikeMask(uRegionB(uB), datABcoripLoad1));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 5}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 5}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 1);
                        
                        %% Co-Firing During No-Ripples
                        datA = find(spikeMask(uRegionA(uA), datABnorip));
                        datB = find(spikeMask(uRegionB(uB), datABnorip));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 2}(cUnitPair(iRa, iRb), ii) = sum(datABnorip);
                        
                        % Load 3 specific
                        datABnoripLoad3 = datABnorip & taskMaskLoad3(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABnoripLoad3));
                        datB = find(spikeMask(uRegionB(uB), datABnoripLoad3));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 6}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 6}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 3);
                        
                        % Load 1 specific
                        datABnoripLoad1 = datABnorip & taskMaskLoad1(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABnoripLoad1));
                        datB = find(spikeMask(uRegionB(uB), datABnoripLoad1));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 7}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 7}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 1);
                        
                        %% All Co-Firing (Regardless of Ripples)
                        datA = find(spikeMask(uRegionA(uA), taskMask(ii, :)));
                        datB = find(spikeMask(uRegionB(uB), taskMask(ii, :)));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 3}(cUnitPair(iRa, iRb), ii) = sum(taskMask(ii, :));
                        
                        % Load 3 specific
                        datABLoad3 = taskMaskLoad3(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABLoad3));
                        datB = find(spikeMask(uRegionB(uB), datABLoad3));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 8}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 8}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 3);
                        
                        % Load 1 specific
                        datABLoad1 = taskMaskLoad1(ii, :);
                        datA = find(spikeMask(uRegionA(uA), datABLoad1));
                        datB = find(spikeMask(uRegionB(uB), datABLoad1));
                        s = datA - datB';
                        whole_trials_unit_cross_region_all{iRa, iRb, 9}(cUnitPair(iRa, iRb), ii) = ...
                            sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                        whole_trials_cross_region_dur{iRa, iRb, 9}(cUnitPair(iRa, iRb), ii) = sum(trials.loads == 1);
                    end
                    
                    % Store unit pair metadata
                    unitPairID{iRa, iRb}(cUnitPair(iRa, iRb), 1) = subj;
                    unitPairID{iRa, iRb}(cUnitPair(iRa, iRb), 2:3) = [uA uB];
                    UnitPairType{iRa, iRb}{1, cUnitPair(iRa, iRb)} = typeA;
                    UnitPairType{iRa, iRb}{2, cUnitPair(iRa, iRb)} = typeB;
                    UnitPairFR{iRa, iRb}(1, cUnitPair(iRa, iRb)) = FR(uA);
                    UnitPairFR{iRa, iRb}(2, cUnitPair(iRa, iRb)) = FR(uB);
                    cUnitPair(iRa, iRb) = cUnitPair(iRa, iRb) + 1;
                end
            end
        end
    end
end

% Clean up unused preallocated space
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) = [];
probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('Done processing.\n')
toc

%%  VISUALIZATION AND STATISTICAL ANALYSIS

%% Setup for Plotting
regionColors = brewermap(12, 'Dark2');
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};
thresh = 0; % Threshold for minimum co-firing rate

% Initialize results matrices
coRegionColors = nan(10, 3);
coRall = [];
noRall = [];
baseAll = [];
Npair = nan(5, 5);
pBaseline = nan(length(regionsPlot), length(regionsPlot), 3, 3);
modBaseline = nan(length(regionsPlot), length(regionsPlot), 3, 3);
keepAll = cell(length(regionsPlot), length(regionsPlot));

%% Figure: Co-Firing Rates Across Task Phases
figure('Position', [6 896 771 677])

for iRa = 1:length(regionsPlot)
    for iRb = iRa:length(regionsPlot)
        % Get indices for this region pair
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        % Aggregate data across hemispheres
        coR = [];
        noR = [];
        base = [];
        fireA = [];
        fireB = [];
        
        for iiA = iRapl
            for iiB = iRbpl
                if iiA == iiB; continue; end
                if isempty(whole_trials_unit_cross_region_all{iiA, iiB, 3}); continue; end
                
                types = UnitPairType{iiA, iiB};
                iiType = strcmp(types(1, :), 'pyr') | strcmp(types(1, :), 'int');
                
                tempBase = [];
                tempCoRcount = [];
                tempNoR = [];
                
                % Calculate rates for each task phase
                for ii = 1:4
                    if ii == 1
                        % Baseline uses overall co-firing as reference
                        tempBase(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 3}(iiType, ii) ./ ...
                                         whole_trials_cross_region_dur{iiA, iiB, 3}(iiType, ii) * 1e3;
                        tempCoRcount(ii, :) = tempBase(ii, :);
                        tempNoR(ii, :) = tempBase(ii, :);
                    else
                        tempBase(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 3}(iiType, ii) ./ ...
                                         whole_trials_cross_region_dur{iiA, iiB, 3}(iiType, ii) * 1e3;
                        tempCoRcount(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 1}(iiType, ii) ./ ...
                                             whole_trials_cross_region_dur{iiA, iiB, 1}(iiType, ii) * 1e3;
                        tempNoR(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 2}(iiType, ii) ./ ...
                                        whole_trials_cross_region_dur{iiA, iiB, 2}(iiType, ii) * 1e3;
                    end
                end
                
                coR = [coR tempCoRcount];
                noR = [noR tempNoR];
                base = [base tempBase];
                fireA = [fireA UnitPairFR{iiA, iiB}(1, iiType)];
                fireB = [fireB UnitPairFR{iiA, iiB}(2, iiType)];
            end
        end
        
        % Create subplot for this region pair
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c, :) = mean([regionColors(iRa, :); regionColors(iRb, :)]);
        subplot(5, 5, c)
        
        % Filter pairs with sufficient activity
        keep = coR(2, :) > thresh | coR(3, :) > thresh | coR(4, :) > thresh;
        keepAll{iRa, iRb} = keep;
        Npair(iRa, iRb) = sum(keep);
        
        % Statistical testing
        for iT = 1:3
            [pBaseline(iRa, iRb, iT, 1), ~, ~] = pairedPermutationTest(base(1, keep)', coR(iT+1, keep)', 10000);
            [pBaseline(iRa, iRb, iT, 3), ~, ~] = pairedPermutationTest(base(1, keep)', noR(iT+1, keep)', 10000);
            [pBaseline(iRa, iRb, iT, 2), ~, ~] = pairedPermutationTest(base(1, keep)', base(iT+1, keep)', 10000);
            
            modBaseline(iRa, iRb, iT, 1) = (median(coR(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
            modBaseline(iRa, iRb, iT, 3) = (median(noR(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
            modBaseline(iRa, iRb, iT, 2) = (median(base(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
        end
        
        % Aggregate across all region pairs
        coRall = [coRall coR(:, keep)];
        noRall = [noRall noR(:, keep)];
        baseAll = [baseAll base(:, keep)];
        
        % Plot data
        xlim([0.5 4.5])
        hl = hline(median(base(1, keep), 'omitnan'), 'k-');
        hl.LineWidth = 0.5;
        
        errobarPatch(1, base(1, keep)', 'k', 0.2);
        errobarPatch([2:4]-0.2, base(2:end, keep)', coRegionColors(c, :)*0.6, 0.1);
        errobarPatch([2:4]+0.2-0.2, coR(2:end, keep)', coRegionColors(c, :), 0.1);
        errobarPatch([2:4]-0.2-0.2, noR(2:end, keep)', coRegionColors(c, :)*0.3, 0.1);
        
        xlim([0.5 4.7])
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'B', 'E1', 'M', 'P'};
        
        if iRa == 1; ylabel('co-fire rate [Hz]'); end
        
        % Print modulation during maintenance
        pct = [median(coR(3, keep), 'omitnan') - median(base(1, keep), 'omitnan')] / ...
              median(base(1, keep), 'omitnan');
        fprintf('%s <--> %s ripples maintenance %.2f\n', regionsPlot{iRa}, regionsPlot{iRb}, pct)
        
        axRng = range(ax.YLim);
        ax.YLim(1) = ax.YLim(1) - 0.2*axRng;
    end
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegions_%s.pdf', tag)))

% Kruskal-Wallis test across conditions
[p, tbl, stats] = kruskalwallis([baseAll(4, :); coRall(3, :); noRall(4, :)]');

%% Statistical Heatmap with FDR Correction
adj_p = nan(size(pBaseline));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pBaseline))] = fdr_bh(pBaseline(~isnan(pBaseline)), 0.05, 'pdep', 'yes');

figure('Position', [1219 1378 600 600]);
ha = tight_subplot(3, 3, [0.05, 0.05], [.05 .05], [.07 .07]);

pArray = [1 0.8 0.5 0.2 0.05 1e-3 1e-4 1e-5];
pSize = [2 3 6 10 15 20 25];

% Create color map for modulation values
cmap = slanCM('RdBu');
cmapRGB = [cmap(1, :); [0.9 0.9 0.9]; cmap(end, :)];
xq = linspace(0, 1, 256*3);
colorPositions = [0, 1/3, 1];
cmap = interp1(colorPositions, cmapRGB, xq, 'linear');

for ii = 1:3
    for jj = 1:3
        adj_pFL = fliplr(adj_p(:, :, jj, ii));
        modBaselineFL = fliplr(modBaseline(:, :, jj, ii));
        c = (ii - 1) * 3 + jj;
        axes(ha(c))
        
        [xx, yy] = find(~isnan(adj_pFL));
        for p = 1:length(xx)
            pl = plot(xx(p), yy(p), 'rs'); hold on;
            
            pVal = adj_pFL(xx(p), yy(p));
            iiP = arrayfun(@(X) pArray(X) > pVal & pArray(X+1) <= pVal, 1:length(pArray)-1);
            if pVal == 1
                iiP(1) = true;
            elseif pVal < 1e-5
                iiP(end) = true;
            end
            pl.MarkerSize = pSize(iiP);
            
            pl.MarkerFaceColor = interp1(linspace(-0.25, 0.5, 256*3), cmap, ...
                                         modBaselineFL(xx(p), yy(p)), 'nearest', 'extrap');
            pl.MarkerEdgeColor = interp1(linspace(-0.25, 0.5, 256*3), cmap, ...
                                         modBaselineFL(xx(p), yy(p)), 'nearest', 'extrap');
            
            % Highlight diagonal
            if xx(p) + yy(p) == 6
                pl = plot(xx(p), yy(p), 'rs'); hold on;
                pl.MarkerEdgeColor = 'g';
                pl.MarkerFaceColor = 'none';
                pl.MarkerSize = max(pSize) + 3;
                pl.LineWidth = 1.5;
            end
        end
        
        cmap = [1 1 1; cmap];
        
        [xx, yy] = find(adj_pFL < 0.05);
        if ~isempty(xx)
            pl = plot(xx, yy, '*'); hold on;
            pl.MarkerSize = 4;
            pl.Color = [0.8 0 0];
        end
        
        ylim([0.5 5.5])
        xlim([0.5 5.5])
        box off
        ax = gca;
        ax.YTickLabel = fliplr(regionsPlot);
        ax.XTickLabel = regionsPlot;
        
        fig = gcf;
        fig.Color = 'w';
    end
end

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_%s.pdf', tag)))

% Create colorbars and legends
figure('Position', [1214 1349 172 212]);
imagesc(zeros(5, 5), [-0.25 0.5]); hold on;
h = colorbar;
caxis([-0.25 0.5]*100);
h.Ticks = [-0.25:0.25:0.5]*100;
colormap(cmap)
fig = gcf;
fig.Color = 'w';
h.Location = 'southoutside';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_modlegend_%s.pdf', tag)))

figure('Position', [1219 1378 600 202]);
for p = 1:length(pSize)
    pl = plot(1+p*(1.3/3), 1, 'rs'); hold on;
    pl.MarkerSize = pSize(p);
    pl.MarkerFaceColor = cmap(end, :);
    pl.MarkerEdgeColor = cmap(end, :);
end
axis off
xlim([0 p])
fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_plegend_%s.pdf', tag)))

%% Figure: Aggregate Box Plot Across All Pairs
clr = brewermap(10, 'Paired');
cBox = [clr(8, :); clr(8, :)*0.6; clr(1, :)];

figure('Position', [370 1318 362 209])
group_num = [ones(1, length(coRall)), 2*ones(1, length(noRall)), 3*ones(1, length(baseAll))];
dat = [coRall(2:end, :), baseAll(2:end, :), noRall(2:end, :)]';
h = daboxplot(dat, 'groups', group_num, 'mean', 1, 'outlier', 0, 'color', cBox, 'xshift', 1, 'boxalpha', 0.5); hold on;
hB = daboxplot(log10(baseAll(1, :)+1)', 'mean', 1, 'outlier', 0, 'color', 'k', 'boxalpha', 0.5, 'boxwidth', 0.6); hold on;

ax = gca;
ax.XTick = [1:4];
ax.XTickLabel = {'B', 'E1', 'M', 'P'};
ax.LineWidth = 1.0;
fig = gcf;
fig.Color = 'w';
xlim([0.5 4.5])
ylabel('co-fire rate [Hz]');
hl = hline(mean(baseAll(1, :)));
hl.LineStyle = '-';

%% Figure: Percent Change Analysis
figure('Position', [682 1313 425 202])

dat = nan(size(baseAll, 1), size(baseAll, 2), 3);
pVal = nan(3, 3);
effectsize = nan(3, 3);

for t = 2:4
    [pVal(1, t-1), ~, effectsize(1, t-1)] = pairedPermutationTest(baseAll(1, :)', baseAll(t, :)', 2000);
    [pVal(2, t-1), ~, effectsize(2, t-1)] = pairedPermutationTest(baseAll(1, :)', coRall(t, :)', 2000);
    [pVal(3, t-1), ~, effectsize(3, t-1)] = pairedPermutationTest(baseAll(1, :)', noRall(t, :)', 2000);
    
    bs = baseAll(1, :);
    xx = coRall(t, :);
    dat(t-1, :, 1) = (xx-bs) / median(bs) * 100;
    effectsize(1, t-1) = median(dat(t-1, :, 1));
    
    xx = baseAll(t, :);
    dat(t-1, :, 2) = (xx-bs) / median(bs) * 100;
    effectsize(2, t-1) = median(dat(t-1, :, 2));
    
    xx = noRall(t, :);
    dat(t-1, :, 3) = (xx-bs) / median(bs) * 100;
    effectsize(3, t-1) = median(dat(t-1, :, 3));
end

% FDR correction
p_adj = nan(size(pVal));
[~, crit_p, adj_ci_cvrg, p_adj(~isnan(pVal))] = fdr_bh(pVal(~isnan(pVal)), 0.05, 'pdep', 'yes');

xlim([1.5 4.5])
hl = hline(0); hold on;
hl.LineStyle = '-';
hl.Color = 'k';

cBox = [clr(8, :); mean([clr(1, :); clr(8, :)]); clr(1, :)];

for ii = 1:3
    bf = bar(h.gpos(ii, 1), median(dat(1, :, ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii, :);
    
    bf = bar(h.gpos(ii, 2), median(dat(2, :, ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii, :);
    
    bf = bar(h.gpos(ii, 3), median(dat(3, :, ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii, :);
    
    quantile(dat(1, :, ii), [0.25 0.75]);
end

ax = gca;
ax.XTick = [1:4];
ax.XTickLabel = {'B', 'E1', 'M', 'P'};
ax.LineWidth = 1.0;
fig = gcf;
fig.Color = 'w';
xlim([1.5 4.5])

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsAll_%s.pdf', tag)))

%%  ANATOMICAL DISTANCE ANALYSIS

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));

fs = 1e3;
ctxParc = {'ACC', 'SMA', 'OFC', 'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};

figure('Position', [1224 1316 328 257]);

% Aggregate co-firing and distance data
x = [];
xMu = [];
yCo = [];
yNo = [];
iRaAll = [];
iRbAll = [];
hemA = [];
hemB = [];
parcA = [];
parcB = [];
withinBundle = [];

for iRa = 1:nRegions
    for iRb = iRa:nRegions
        parcelA = regions{iRa}(2:end);
        parcelB = regions{iRb}(2:end);
        
        types = UnitPairType{iRa, iRb};
        iiType = strcmp(types(1, :), 'pyr') | strcmp(types(1, :), 'int');
        fireA = UnitPairFR{iRa, iRb}(1, iiType);
        fireB = UnitPairFR{iRa, iRb}(2, iiType);
        
        tempCoRcount = [];
        tempCoRdur = [];
        tempNoRcount = [];
        tempNoRdur = [];
        tempBaseCount = [];
        
        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb, 1}(iiType, :), 2);
        tempCoRcount = [tempCoRcount temp'];
        temp = sum(whole_trials_cross_region_dur{iRa, iRb, 1}(iiType, :), 2);
        tempCoRdur = [tempCoRdur temp'];
        
        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb, 2}(iiType, :), 2);
        tempNoRcount = [tempNoRcount temp'];
        temp = sum(whole_trials_cross_region_dur{iRa, iRb, 2}(iiType, :), 2);
        tempNoRdur = [tempNoRdur temp'];
        
        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb, 3}(iiType, :), 2);
        tempBaseCount = [tempBaseCount temp'];
        
        tCo = ((tempCoRcount ./ tempCoRdur) * 1e3);
        tNo = ((tempNoRcount ./ tempNoRdur) * 1e3);
        
        z = [(tCo(tempBaseCount >= 0)) (tNo(tempBaseCount >= 0))];
        tCo = z(1:length(z)/2);
        tNo = z((length(z)/2)+1:end);
        
        hemA = [hemA repmat({regions{iRa}(1)}, [1 length(tCo)])];
        hemB = [hemB repmat({regions{iRb}(1)}, [1 length(tCo)])];
        parcA = [parcA repmat({parcelA}, [1 length(tCo)])];
        parcB = [parcB repmat({parcelB}, [1 length(tCo)])];
        samebund = iRa == iRb;
        withinBundle = [withinBundle repmat(samebund, [1 length(tCo)])];
        
        broadmanA = broadman{strcmp(ctxParc, parcelA)};
        HCPa = sprintf('%s_%s', hemA{end}, broadmanA);
        broadmanB = broadman{strcmp(ctxParc, parcelB)};
        HCPb = sprintf('%s_%s', hemB{end}, broadmanB);
        
        d = tractLenStruct.tractLengths(strcmp(tractLenStruct.parcelIDs, HCPa), ...
                                        strcmp(tractLenStruct.parcelIDs, HCPb));
        
        yCo = [yCo tCo];
        yNo = [yNo tNo];
        iRaAll = [iRaAll repmat(iRa, [1 length(tCo)])];
        iRbAll = [iRbAll repmat(iRb, [1 length(tCo)])];
        xMu = [xMu d];
        x = [x repmat(d, [1 length(tCo)])];
        
        yCoMu = [yCoMu mean(tCo)];
        yNoMu = [yNoMu mean(tNo)];
    end
end

% Bin by distance and compute means
yCoMu = [];
yNoMu = [];
yCoSEM = [];
yNoSEM = [];
xMu = [];
xBin = quantile(x, 6);
xBin = [0 xBin];
xBin = [xBin max(x)];

for iB = 2:length(xBin)
    ii = x >= xBin(iB-1) & x <= xBin(iB);
    pctDelta = yCo(ii);
    yCoMu = [yCoMu mean(pctDelta, 'omitnan')];
    yNoMu = [yNoMu mean(yNo(ii))];
    yCoSEM = [yCoSEM std(pctDelta, 'omitnan')/sqrt(sum(ii))];
    yNoSEM = [yNoSEM std(yNo(ii))/sqrt(sum(ii))];
    xMu = [xMu mean([xBin(iB-1) xBin(iB)])];
end

clr = brewermap(10, 'Paired');

[bl2, bf] = boundedline(xMu, yCoMu, yCoSEM, 'rs-'); hold on;
bl2.Color = clr(8, :);
bl2.MarkerFaceColor = clr(8, :);
bf.FaceColor = clr(8, :);
bf.FaceAlpha = 0.3;

[bl2, bf] = boundedline(xMu, yNoMu, yNoSEM, 'bs-'); hold on;
bl2.Color = clr(1, :);
bl2.MarkerFaceColor = clr(1, :);
bf.FaceColor = clr(1, :);
bf.FaceAlpha = 0.3;

ax = gca;
ax.FontSize = 11;
ylh = ylabel('co-fire / $\sqrt{FR_a FR_b}$');
ylh.Interpreter = 'latex';
xlabel('distance [mm]')

fig = gcf;
fig.Color = 'w';

xlim([0 250])

% Correlation analysis
[RHO, PVAL] = corr(x', yNo', 'type', 'Kendall');
[RHO, PVAL] = corr(x', yCo', 'type', 'Kendall');
[RHO, PVAL] = corr(x', [yCo-yNo]', 'type', 'Kendall');

[H, P] = signrank(yCo, yNo, 'tail', 'both');
q = quantile((yCo-yNo)./yNo, [0.25 0.5 0.75]);

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireDistanceGeoMean_%s.pdf', tag)))

%% Figure: Co-Firing by Anatomical Connection Type
figure('Position', [100 368 round(535*(2/3)) round(232*(2/3))]);

num_variables = 6;
num_conditions = 2;

% Create custom x positions for grouped boxplots
pos = zeros(num_variables * num_conditions, 1);
gap_within = 0.3;
gap_between = 0.8;

for var = 1:num_variables
    base_pos = (var-1) * (gap_between);
    pos((var-1)*num_conditions + 1) = base_pos;
    pos((var-1)*num_conditions + 2) = base_pos + gap_within;
end

all_data = [];
all_groups = [];

% Within bundle
datCond1 = yCo(withinBundle == 1);
all_data = [all_data, datCond1];
idx = (1-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 1);
all_data = [all_data, datCond2];
idx = (1-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

% AMY to frontal
datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'AMY'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond1];
idx = (2-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'AMY'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond2];
idx = (2-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H, P] = signrank(datCond1, datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H)

% HIP to frontal
datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'HIP'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond1];
idx = (3-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'HIP'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond2];
idx = (3-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H, P] = signrank(datCond1, datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H)

% HIP-AMY
datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'HIP', 'AMY'}) & ismember(parcB, {'HIP', 'AMY'}));
all_data = [all_data, datCond1];
idx = (4-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'HIP', 'AMY'}) & ismember(parcB, {'HIP', 'AMY'}));
all_data = [all_data, datCond2];
idx = (4-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H, P] = signrank(datCond1, datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H)

% Frontal ipsilateral
datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ...
               ismember(parcB, {'OFC', 'ACC', 'SMA'}) & strcmp(hemA, hemB));
all_data = [all_data, datCond1];
idx = (5-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ...
               ismember(parcB, {'OFC', 'ACC', 'SMA'}) & strcmp(hemA, hemB));
all_data = [all_data, datCond2];
idx = (5-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H, P] = signrank(datCond1, datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H)

% Frontal contralateral
datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ...
               ismember(parcB, {'OFC', 'ACC', 'SMA'}) & ~strcmp(hemA, hemB));
all_data = [all_data, datCond1];
idx = (6-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ...
               ismember(parcB, {'OFC', 'ACC', 'SMA'}) & ~strcmp(hemA, hemB));
all_data = [all_data, datCond2];
idx = (6-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H, P] = signrank(datCond1, datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H)

% Create boxplot
bx = boxplot(all_data, all_groups, 'positions', pos, 'width', 0.2, ...
             'Colors', [0 0 0], 'Symbol', '', 'Whisker', 0);

tick_pos = zeros(num_variables, 1);
tick_labels = cell(num_variables, 1);

for var = 1:num_variables
    tick_pos(var) = pos((var-1)*num_conditions + 1) + gap_within/2;
end

set(gca, 'XTick', tick_pos);
set(gca, 'XTickLabel', tick_labels);
set(gca, 'XLim', [min(pos)-0.5, max(pos)+0.5]);
set(gca, 'LineWidth', 0.75);
set(gca, 'FontSize', 7);
set(bx, 'LineWidth', 0.75);

% Customize box colors
for var = 1:num_variables
    h = findobj(gca, 'Tag', 'Box');
    for cond = 1:num_conditions
        idx = (var-1)*num_conditions + (num_conditions-cond+1);
        if cond == 1
            patch(get(h(idx), 'XData'), get(h(idx), 'YData'), clr(8, :), 'FaceAlpha', 0.5, 'LineWidth', 0.75);
        else
            patch(get(h(idx), 'XData'), get(h(idx), 'YData'), clr(1, :), 'FaceAlpha', 0.5, 'LineWidth', 0.75);
        end
    end
end

ylabel('co-firing rate [Hz]', 'FontSize', 7);
box off

fig = gcf;
fig.Color = 'w';
ylim([-0.1 1])

qWithin = quantile(yCo(withinBundle == 1), [0.25 0.5 0.75]);
qAcross = quantile(yCo(withinBundle == 0), [0.25 0.5 0.75]);

[H, P] = ranksum(yCo(withinBundle == 1), yCo(withinBundle == 0), 'tail', 'both');
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireSummary_%s.pdf', tag)))

%%  LOAD-DEPENDENT MODULATION ANALYSIS

EdgeColors = [clr(8, :); mean([clr(1, :); clr(8, :)]); clr(1, :)];

% Initialize results
Pval = nan(length(regionsPlot), length(regionsPlot), 3, 3);
modulAll = nan(length(regionsPlot), length(regionsPlot), 3, 3);
modulAllPairs = cell(1, 2);

spacing = 0.4;
strng = {'M', 'P'};

% Analyze modulation for maintenance and probe phases
for t = 3:4
    figure('Position', [1 847 410 296])
    xx_jitt = [];
    YY = [];
    
    for iRa = 1:length(regionsPlot)
        for iRb = iRa:length(regionsPlot)
            iRapl = find(contains(regions, regionsPlot{iRa}));
            iRbpl = find(contains(regions, regionsPlot{iRb}));
            
            coR = [];
            noR = [];
            base = [];
            
            for iiA = iRapl
                for iiB = iRbpl
                    if iiA == iiB; continue; end
                    if isempty(whole_trials_unit_cross_region_all{iiA, iiB, 3}); continue; end
                    
                    types = UnitPairType{iiA, iiB};
                    iiType = strcmp(types(1, :), 'pyr') | strcmp(types(1, :), 'int');
                    
                    fprintf('%s %s\n', regions{iiA}, regions{iiB})
                    
                    tempBase = [];
                    tempCoRcount = [];
                    tempNoR = [];
                    
                    for ii = 1:2
                        tempBase(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 7+ii}(iiType, t) ./ ...
                                         whole_trials_cross_region_dur{iiA, iiB, 7+ii}(iiType, t);
                        tempCoRcount(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 3+ii}(iiType, t) ./ ...
                                             whole_trials_cross_region_dur{iiA, iiB, 3+ii}(iiType, t);
                        tempNoR(ii, :) = whole_trials_unit_cross_region_all{iiA, iiB, 5+ii}(iiType, t) ./ ...
                                        whole_trials_cross_region_dur{iiA, iiB, 5+ii}(iiType, t);
                    end
                    
                    coR = [coR tempCoRcount];
                    noR = [noR tempNoR];
                    base = [base tempBase];
                end
            end
            
            c = (iRb - 1) * 5 + iRa;
            coRegionColors(c, :) = mean([regionColors(iRa, :); regionColors(iRb, :)]);
            
            keep = keepAll{iRa, iRb};
            
            % Calculate load modulation
            baseModul = [base(1, keep)-base(2, keep)]./mean([base(1, keep) base(2, keep)]) * 100;
            coRmodul = [coR(1, keep)-coR(2, keep)]./mean([coR(1, keep) coR(2, keep)]) * 100;
            noRmodul = [noR(1, keep)-noR(2, keep)]./mean([noR(1, keep) noR(2, keep)]) * 100;
            modul = [coRmodul; baseModul; noRmodul];
            
            modulAllPairs{t-2} = [modulAllPairs{t-2}, modul];
            
            YY(:, iRa, iRb) = mean(modul, 2, 'omitnan');
            
            xx = [(t-2)*3+1:(t-2)*3+3]+((t-2)*spacing);
            
            xx_jitt(:, iRa, iRb) = (1:3) + ((rand(1) - 0.5)*0.25);
            pl = plot(xx_jitt(:, iRa, iRb), YY(:, iRa, iRb), '-'); hold on;
            pl.Color = coRegionColors(c, :);
            pl.LineWidth = 1.5;
            
            xlim([0 3.5])
            fig = gcf;
            fig.Color = 'w';
            hl = hline(0);
            hl.LineStyle = '-';
            hl.Color = 'k';
            hl.LineWidth = 1;
            uistack(hl, 'bottom')
            ax = gca;
            ax.XTick = [(3-2)*3+2+((3-2)*spacing) (4-2)*3+2+((4-2)*spacing)];
            box off
            
            % Statistical tests
            [H, P, CI, STATS] = ttest(base(1, keep), base(2, keep), 'tail', 'right');
            Pval(iRa, iRb, 2, t-1) = P;
            [H, P, CI, STATS] = ttest(noR(1, keep), noR(2, keep), 'tail', 'right');
            Pval(iRa, iRb, 3, t-1) = P;
            [H, P, CI, STATS] = ttest(coR(1, keep), coR(2, keep), 'tail', 'right');
            Pval(iRa, iRb, 1, t-1) = P;
            
            modulAll(iRa, iRb, 1, t-1) = YY(1);
            modulAll(iRa, iRb, 3, t-1) = YY(3);
            modulAll(iRa, iRb, 2, t-1) = YY(2);
            
            fprintf('\n')
        end
    end
    
    % Add markers for significant results
    pValueColorMap = [
        0.8, 0.8, 0.8;    % p > 0.05
        0.1, 0.7, 0.9;    % p < 0.05
        0.4, 0.0, 0.7;    % p < 0.01
        0.8, 0.0, 0.0     % p < 0.001
    ];
    
    adj_p = nan(size(Pval));
    [h, crit_p, adj_ci_cvrg, adj_p(~isnan(Pval))] = fdr_bh(Pval(~isnan(Pval)), 0.05, 'pdep', 'yes');
    
    for iRa = 1:length(regionsPlot)
        for iRb = iRa:length(regionsPlot)
            c = (iRb - 1) * 5 + iRa;
            
            pReg = adj_p(iRa, iRb, :, t-1);
            for ii = 1:3
                pl = plot(xx_jitt(ii, iRa, iRb), YY(ii, iRa, iRb), 'o'); hold on;
                pl.MarkerSize = 7;
                pl.LineWidth = 1.5;
                pl.Color = coRegionColors(c, :);
                
                if pReg(ii) < 0.05 && pReg(ii) >= 0.01
                    pl.MarkerFaceColor = pValueColorMap(2, :);
                elseif pReg(ii) < 0.01 && pReg(ii) >= 0.001
                    pl.MarkerFaceColor = pValueColorMap(3, :);
                elseif pReg(ii) < 0.001
                    pl.MarkerFaceColor = pValueColorMap(4, :);
                else
                    pl.MarkerFaceColor = pValueColorMap(1, :);
                end
            end
        end
    end
    
    ax = gca;
    ax.XTick = 1:3;
    ax.XTickLabel = {'coR', 'all', 'noR'};
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
    
    savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskLoadAcrossRegionsSingle%s_%s.pdf', strng{t-2}, tag)))
end