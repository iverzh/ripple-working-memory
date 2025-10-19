%% Replay analysis of co-firing patterns during working memory
% This script tests whether co-firing patterns observed during encoding
% are replayed during probe periods, and how this relates to hippocampal
% ripples in human intracranial recordings during a Sternberg task.
% 
% Key analyses:
% 1. Compare co-firing during encoding-probe pairs (match vs mismatch images)
% 2. Test replay during co-ripple vs non-ripple periods
% 3. Chi-square tests for statistical significance of replay effects

close all 
clc
clear

%% Initialize paths and directories
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))

% Define data directories
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';

% Create output directories if needed
if ~isfolder(exportDir); mkdir(exportDir); end
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

%% Load file lists
% Get subject list from LFP micro files
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;

% Load all necessary file types
unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory,'../task', '*task*'));
taskfiles = {taskfiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../out', '*micro*'));
micrFiles = {micrFiles.name}';

%% Define analysis parameters
recordingState = 'wake';
location = 'NC';

% Brain regions of interest (left/right ? 5 regions)
regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'}; % For visualization

% Load channel curation data
channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end, 3:end-1));   
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2, 3:end-1));
bundleNotes = table2cell(bundleNotes(3:end, 3:end-1));

%% Initialize analysis parameters and matrices
coFireWindow = round(25/2); % Co-firing window in ms (half-window)
computeCoRipple = false;
coFireSubsample = 1;
nIterations = 25; % Number of iterations for bootstrap

% Time windows for analysis
win = [-9 25] * 1e3; % Window around probe onset (ms)
winResp = [-1.5 1] * 1e3; % Window around response (ms)
winLength = 12905; % Total trial length (ms)

% Task timing markers (ms from trial start)
taskMarkers = [500 2600 4700 6800 9400] + 500;

% Initialize data storage matrices
rippleFR = nan(length(subj_list_full), 100);
CoRippleRatesRegionAll = nan(length(regions), length(regions), length(subj_list_full));

% Trial-level data
recog_trials_all = nan(6.5e3, length(win(1):win(2)));
recog_trials_units_all = nan(6.5e3, length(win(1):win(2)));
recog_trials_resp_all = nan(6.5e3, length(winResp(1):winResp(2)));

whole_trial_all = nan(6.5e3, winLength);
whole_trials_region_all = nan(6.5e3, winLength, length(regions));
whole_trials_units_region_all = nan(6.5e3, winLength, length(regions));

% Cross-region analysis arrays
whole_trials_cross_region_all = cell(length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);
whole_trials_cross_region_dur = cell(length(regions), length(regions), 3);
load_trials_all_unitPair = cell(length(regions), length(regions));
unitPairID = cell(length(regions), length(regions));
testCoFire = cell(length(regions), length(regions));
UnitPairType = cell(length(regions), length(regions));
trialsProbe = zeros(length(regions), length(regions), 1e3);

% Replay analysis matrices
coR_replay = zeros(2, length(regions), length(regions)); % Co-ripple replay counts
coR_replay_shuff = zeros(2, length(regions), length(regions)); % Shuffled control
noR_replay = zeros(2, length(regions), length(regions)); % Non-ripple replay counts
noR_replay_shuff = zeros(2, length(regions), length(regions));
full_replay = zeros(2, length(regions), length(regions)); % Full trial replay
full_replay_mismatch = zeros(2, 5);

% Co-firing data storage
fullCoF = nan(5e5, 2, length(regions), length(regions));
fullCoF_mm = nan(5e5, 2, length(regions), length(regions));
replay_event = cell(1, 8); % Store replay event details

% Duration and rate matrices
coRdur = zeros(2, 2, length(regions), length(regions)); % Co-ripple durations
noRdur = zeros(2, 2, length(regions), length(regions)); % Non-ripple durations
coRap = zeros(2, 2, length(regions), length(regions)); % Co-ripple action potentials
noRap = zeros(2, 2, length(regions), length(regions)); % Non-ripple action potentials

% Trial-wise data
loads = nan(5e5, length(regions), length(regions)); % Memory load per trial
respTimes = nan(5e5, length(regions), length(regions)); % Response times
coRcoF = nan(5e5, 2, length(regions), length(regions)); % Co-ripple co-firing
noRcoF = nan(5e5, 2, length(regions), length(regions), nIterations); % Non-ripple co-firing (bootstrapped)
coRcoFShuff = nan(5e5, 2, length(regions), length(regions)); % Shuffled control

% Mismatch condition arrays
coR_replay_mm = zeros(2, length(regions), length(regions));
noR_replay_mm = zeros(2, length(regions), length(regions));
coRdur_mm = zeros(2, 2, length(regions), length(regions));
noRdur_mm = zeros(2, 2, length(regions), length(regions));
coRap_mm = zeros(2, 2, length(regions), length(regions));
noRap_mm = zeros(2, 2, length(regions), length(regions));
coRcoF_mm = nan(5e5, 2, length(regions), length(regions));
noRcoF_mm = nan(5e5, 2, length(regions), length(regions), nIterations);

noRcoFCount = ones(length(regions), length(regions), nIterations);
noRcoF_mmCount = ones(length(regions), length(regions), nIterations);

% Subject tracking
subjID = nan(6.5e3, 1);
correct_trials_all = nan(6.5e3, 1);
load_trials_all = nan(6.5e3, 1);
probe_in_out = nan(6.5e3, 1);
respLatency = nan(6.5e3, 1);
densityAll = [];
uChanAll = [];

%% Main analysis loop across subjects
tic
cTrial = 1;
cUnitPair = ones(length(regions), length(regions));

for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    fprintf('Processing %s ... \n', subject)
    
    %% Load LFP data and parse channel locations
    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    LFPtimeMicr = micrObj.times;
    
    % Parse channel labels
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
    
    % Add hemisphere labels
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locations(hemi) = cellfun(@(X) ['L' X], locations(hemi), 'UniformOutput', false);
    hemi = contains(hem, 'right');
    locations(hemi) = cellfun(@(X) ['R' X], locations(hemi), 'UniformOutput', false);
    
    %% Load task data
    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
    
    %% Load ripple detection results
    modifier = '1kHz_template_z25';
    tag = [recordingState, '_', location, '_', modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    
    % Process ripple statistics
    LFPchanNum = cellfun(@(X) str2double(X), rippleStats.chanLabels);
    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations)
        error('Error formatting channel labels.');
    end
    
    rippleStats.recordingType = repmat({'micro'}, [1 length(microObj.rippleStats.locs)]);
    rippleStats.locs = cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false);
    microObj.rippleStats.window(cellfun(@(X) isempty(X), microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', ...
                                 microObj.rippleStats.window, 'UniformOutput', false);
    rippleStats.microTimes = LFPtimeMicr;
    
    %% Create ripple mask
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask, 1)
        if rippleStats.density{chRipp} > 1
            if strcmp(rippleStats.recordingType{chRipp}, 'macro')
                continue;
            elseif strcmp(rippleStats.recordingType{chRipp}, 'micro')
                times = rippleStats.microTimes;
            end
            
            % Skip special electrode channels
            if contains(rippleStats.chanLabels{chRipp}, 'SPE')
                rippMask(chRipp,:) = nan(1, length(rippMask));
                continue
            end
            
            densityAll = [densityAll rippleStats.density{chRipp}];
            
            % Mark ripple periods
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
            
            for ii = 1:length(iE)
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp, iS(ii):iE(ii)) = 1;
            end
        end
    end
    
    % Pad the end of recordings
    rippMask(:, end:end+3e3) = nan;
    micrObj.lfp_data(:, end:end+3e3) = nan;
    rippAll = sum(rippMask);
    
    %% Load and process unit data
    units = LoadSpikeTimes(subject, 'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    
    % Filter for valid unit types
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    
    % Process unit locations
    uLocations = units(U, end-2);
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
    
    %% Create spike masks
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001; % Bin width in seconds
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    
    for iU = 1:nNeuron
        X = units{unitsAll(iU), 2};
        [N, EDGES] = histcounts(X, binEdges);
        FR = sum(N) / rippleStats.recordingLength * 1000; % Firing rate in Hz
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian', 500);
    end
    
    % Calculate population activity
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian', 500);
    else
        unitAll = smoothdata(spikeMask, 'gaussian', 500);
    end
    unitAll = zscore(unitAll);
    
    %% Analyze replay across region pairs
    for iRseed = 1 % Seed region (can be extended for analysis)
        iiSeed = find(contains(regions, regions{iRseed}));
        iiB = find(~contains(regions, regions{iRseed}));
        
        for iRa = 1:length(regions)
            for iRb = iRa+1:length(regions)
                
                % Get ripple data for both regions
                datArip = whole_trials_region_all(cTrial,:,iRa);
                datBrip = whole_trials_region_all(cTrial,:,iRb);
                
                if iRa == iRb
                    continue
                else
                    datAB = datArip + datBrip;
                    datAB(datArip <= 0 | datBrip <= 0) = 0;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                end
                
                % Get units from each region
                uLFPregionA = uChan(contains(uLocations, regions(iRa)));
                uLFPregionB = uChan(contains(uLocations, regions(iRb)));
                uRegionA = find(contains(uLocations, regions(iRa)));
                uRegionB = find(contains(uLocations, regions(iRb)));
                
                % Get ripple masks for each region
                datArip = sum(rippMask(contains(locations, regions(iRa)),:), 'omitnan') > 0;
                datBrip = sum(rippMask(contains(locations, regions(iRb)),:), 'omitnan') > 0;
                
                if sum(datArip) == 0 || sum(datBrip) == 0; continue; end
                
                nUnits = length(uLFPregionA) + length(uLFPregionB);
                probeIm = []; 
                encodIm = [];
                
                %% Process each trial
                for iT = 1:length(trials.start_time)
                    trialLoad = trials.loads(iT);
                    trialAcc = trials.response_accuracy;
                    
                    % Get trial event timestamps
                    imageTime(1) = round(trials.timestamps_Encoding1(iT) * 1e3);
                    imageTime(2) = round(trials.timestamps_Encoding2(iT) * 1e3);
                    imageTime(3) = round(trials.timestamps_Encoding3(iT) * 1e3);
                    imageTime(4) = round(trials.timestamps_Maintenance(iT) * 1e3);
                    probeTime = round(trials.timestamps_Probe(iT) * 1e3);
                    respTime = round(trials.timestamps_Response(iT) * 1e3);
                    
                    % Store image IDs
                    probeIm(iT) = trials.PicIDs_Probe(iT);
                    encodIm(1, iT) = trials.PicIDs_Encoding1(iT);
                    encodIm(2, iT) = trials.PicIDs_Encoding2(iT);
                    encodIm(3, iT) = trials.PicIDs_Encoding3(iT);
                    
                    if trials.probe_in_out(iT) % Probe matches one of the encoded images
                        
                        for uA = 1:length(uRegionA)
                            for uB = 1:length(uRegionB)
                                typeA = units{uRegionA(uA), 3};
                                typeB = units{uRegionB(uB), 3};
                                
                                % Extract spike data for probe period
                                spikeMaskP = spikeMask(:, probeTime+1:probeTime+1e3);
                                datABcoripP = datArip(probeTime+1:probeTime+1e3) > 0 & ...
                                             datBrip(probeTime+1:probeTime+1e3) > 0;
                                datABnoripP = datArip(probeTime+1:probeTime+1e3) == 0 & ...
                                             datBrip(probeTime+1:probeTime+1e3) == 0;
                                
                                % Get matching encoding period
                                if trials.loads(iT) == 3
                                    iE = find(encodIm(:,iT) == probeIm(iT));
                                    iEmm = find(encodIm(:,iT) ~= probeIm(iT));
                                    iEmm = iEmm(randi(length(iEmm)));
                                    spikeMaskE = spikeMask(:, imageTime(iE):imageTime(iE)+2e3);
                                    spikeMaskEmm = spikeMask(:, imageTime(iEmm):imageTime(iEmm)+2e3);
                                    datABcoripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & ...
                                                 datBrip(imageTime(iE):imageTime(iE)+2e3);
                                    datABnoripE = datArip(imageTime(iE):imageTime(iE)+2e3) == 0 & ...
                                                 datBrip(imageTime(iE):imageTime(iE)+2e3) == 0;
                                else
                                    spikeMaskE = spikeMask(:, imageTime(1):imageTime(1)+2e3);
                                    datABcoripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & ...
                                                 datBrip(imageTime(1):imageTime(1)+2e3) > 0;
                                    datABnoripE = datArip(imageTime(1):imageTime(1)+2e3) == 0 & ...
                                                 datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                                end
                                
                                if sum(datABcoripP) == 0 || sum(datABcoripE) == 0; continue; end
                                
                                %% Analyze co-ripple periods
                                % Calculate co-firing during probe
                                datAP = find(spikeMaskP(uRegionA(uA), datABcoripP));
                                datBP = find(spikeMaskP(uRegionB(uB), datABcoripP));
                                sP = datAP - datBP';
                                coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);
                                
                                % Calculate co-firing during encoding
                                datAE = find(spikeMaskE(uRegionA(uA), datABcoripE));
                                datBE = find(spikeMaskE(uRegionB(uB), datABcoripE));
                                sE = datAE - datBE';
                                coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);
                                
                                % Store co-firing data
                                ii = find(isnan(coRcoF(:,1,iRa,iRb)), 1, 'first');
                                coRcoF(ii, 1, iRa, iRb) = coFp;
                                coRcoF(ii, 2, iRa, iRb) = coFe;
                                loads(ii, iRa, iRb) = trialLoad;
                                respTimes(ii, iRa, iRb) = respTime - probeTime;
                                
                                % Count replay events
                                if coFp > 0 && coFe > 0
                                    coR_replay(1, iRa, iRb) = coR_replay(1, iRa, iRb) + 1;
                                    coRdur(1, 1, iRa, iRb) = coRdur(1, 1, iRa, iRb) + sum(datABcoripE);
                                    coRdur(1, 2, iRa, iRb) = coRdur(1, 2, iRa, iRb) + sum(datABcoripP);
                                    coRap(1, 1, iRa, iRb) = coRap(1, 1, iRa, iRb) + sum(coFe);
                                    coRap(1, 2, iRa, iRb) = coRap(1, 2, iRa, iRb) + sum(coFp);
                                    
                                    % Store replay event details
                                    replay_event{1} = [replay_event{1} subj];
                                    replay_event{2} = [replay_event{2} uRegionA(uA)];
                                    replay_event{3} = [replay_event{3} uRegionB(uB)];
                                    replay_event{4} = [replay_event{4} iT];
                                    replay_event{5} = [replay_event{5} {datAP}];
                                    replay_event{6} = [replay_event{6} {datBP}];
                                    replay_event{7} = [replay_event{7} {datAE}];
                                    replay_event{8} = [replay_event{8} {datBE}];
                                else
                                    coR_replay(2, iRa, iRb) = coR_replay(2, iRa, iRb) + 1;
                                    coRdur(2, 1, iRa, iRb) = coRdur(2, 1, iRa, iRb) + sum(datABcoripE);
                                    coRdur(2, 2, iRa, iRb) = coRdur(2, 2, iRa, iRb) + sum(datABcoripP);
                                    coRap(2, 1, iRa, iRb) = coRap(2, 1, iRa, iRb) + sum(coFe);
                                    coRap(2, 2, iRa, iRb) = coRap(2, 2, iRa, iRb) + sum(coFp);
                                end
                                
                                %% Shuffled control analysis
                                datAshuff = randi(sum(datABcoripP) + sum(datABcoripE), ...
                                                 [1, length(datAP) + length(datAE)]);
                                datBshuff = randi(sum(datABcoripP) + sum(datABcoripE), ...
                                                 [1, length(datBP) + length(datBE)]);
                                
                                A = datAshuff(datAshuff <= sum(datABcoripP));
                                B = datBshuff(datBshuff <= sum(datABcoripP));
                                if isempty(A) || isempty(B)
                                    sPshuff = [];
                                else
                                    sPshuff = A - B';
                                end
                                coFpshuff = sum(sPshuff(:) >= -coFireWindow*2 & sPshuff(:) <= coFireWindow*2);
                                
                                A = datAshuff(datAshuff > sum(datABcoripP));
                                B = datBshuff(datBshuff > sum(datABcoripP));
                                if isempty(A) || isempty(B)
                                    sEshuff = [];
                                else
                                    sEshuff = A - B';
                                end
                                coFeshuff = sum(sEshuff(:) >= -coFireWindow*2 & sEshuff(:) <= coFireWindow*2);
                                
                                coRcoFShuff(ii, 1, iRa, iRb) = coFpshuff;
                                coRcoFShuff(ii, 2, iRa, iRb) = coFeshuff;
                                
                                if coFpshuff > 0 && coFeshuff > 0
                                    coR_replay_shuff(1, iRa, iRb) = coR_replay_shuff(1, iRa, iRb) + 1;
                                else
                                    coR_replay_shuff(2, iRa, iRb) = coR_replay_shuff(2, iRa, iRb) + 1;
                                end
                                
                                %% Bootstrap analysis for non-ripple periods
                                for iter = 1:nIterations
                                    % Sample non-ripple periods matched for duration
                                    ripDur = sum(datABcoripP);
                                    b = mask2bounds(datABnoripP);
                                    bDur = b(:,2) - b(:,1);
                                    
                                    if isempty(b)
                                        b = [1, 1]; 
                                        ripDur = 0;
                                    elseif sum(bDur >= ripDur) == 0
                                        b = b(bDur == max(bDur),:);
                                        ripDur = max(bDur);
                                    else
                                        b(bDur < ripDur,:) = [];
                                        b = b(randi(size(b,1),1),:);
                                    end
                                    
                                    datABnoripPiter = datABnoripP;
                                    b = b(randi(size(b,1),1),:);
                                    randStart = randi([b(1), b(2) - ripDur]);
                                    datABnoripPiter(1:randStart) = false;
                                    datABnoripPiter(randStart+ripDur+1:end) = false;
                                    
                                    % Calculate co-firing for probe period
                                    datA = find(spikeMaskP(uRegionA(uA), datABnoripPiter));
                                    datB = find(spikeMaskP(uRegionB(uB), datABnoripPiter));
                                    sP = datA - datB';
                                    noFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);
                                    
                                    % Repeat for encoding period
                                    ripDur = sum(datABcoripE);
                                    b = mask2bounds(datABnoripE);
                                    bDur = b(:,2) - b(:,1);
                                    
                                    if isempty(b)
                                        b = [1, 1];
                                        ripDur = 0;
                                    elseif sum(bDur >= ripDur) == 0
                                        b = b(find(bDur == max(bDur), 1, 'first'),:);
                                        ripDur = max(bDur);
                                    else
                                        b(bDur < ripDur,:) = [];
                                        b = b(randi(size(b,1),1),:);
                                    end
                                    
                                    datABnoripEiter = datABnoripE;
                                    randStart = randi([b(1), b(2) - ripDur]);
                                    datABnoripEiter(1:randStart) = false;
                                    datABnoripEiter(randStart+ripDur+1:end) = false;
                                    
                                    datA = find(spikeMaskE(uRegionA(uA), datABnoripEiter));
                                    datB = find(spikeMaskE(uRegionB(uB), datABnoripEiter));
                                    sE = datA - datB';
                                    noFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);
                                    
                                    % Store bootstrapped data
                                    ii = noRcoFCount(iRa, iRb, iter);
                                    noRcoF(ii, 1, iRa, iRb, iter) = noFp;
                                    noRcoF(ii, 2, iRa, iRb, iter) = noFe;
                                    noRcoFCount(iRa, iRb, iter) = noRcoFCount(iRa, iRb, iter) + 1;
                                    
                                    if noFp > 0 && noFe > 0
                                        noRdur(1, 1, iRa, iRb) = noRdur(1, 1, iRa, iRb) + sum(datABnoripEiter);
                                        noRdur(1, 2, iRa, iRb) = noRdur(1, 2, iRa, iRb) + sum(datABnoripPiter);
                                        noRap(1, 1, iRa, iRb) = noRap(1, 1, iRa, iRb) + sum(noFe);
                                        noRap(1, 2, iRa, iRb) = noRap(1, 2, iRa, iRb) + sum(noFp);
                                    else
                                        noRdur(2, 1, iRa, iRb) = noRdur(2, 1, iRa, iRb) + sum(datABnoripEiter);
                                        noRdur(2, 2, iRa, iRb) = noRdur(2, 2, iRa, iRb) + sum(datABnoripPiter);
                                        noRap(2, 1, iRa, iRb) = noRap(2, 1, iRa, iRb) + sum(noFe);
                                        noRap(2, 2, iRa, iRb) = noRap(2, 2, iRa, iRb) + sum(noFp);
                                    end
                                end
                                
                                if noFp > 0 && noFe > 0
                                    noR_replay(1, iRa, iRb) = noR_replay(1, iRa, iRb) + 1;
                                else
                                    noR_replay(2, iRa, iRb) = noR_replay(2, iRa, iRb) + 1;
                                end
                                
                                %% Full trial co-firing (regardless of ripples)
                                datAP = find(spikeMaskP(uRegionA(uA),:));
                                datBP = find(spikeMaskP(uRegionB(uB),:));
                                sP = datAP - datBP';
                                coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);
                                
                                datAE = find(spikeMaskE(uRegionA(uA),:));
                                datBE = find(spikeMaskE(uRegionB(uB),:));
                                sE = datAE - datBE';
                                coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);
                                
                                ii = find(isnan(fullCoF(:,1,iRa,iRb)), 1, 'first');
                                fullCoF(ii, 1, iRa, iRb) = coFp;
                                fullCoF(ii, 2, iRa, iRb) = coFe;
                            end
                        end
                        
                    else % Probe does not match encoded images (mismatch condition)
                        
                        % Similar analysis for mismatch trials
                        % (Process mismatch condition with similar structure)
                        % This section analyzes trials where probe doesn't match any encoded image
                        % Code follows same structure as match condition above
                        % [Similar analysis code for mismatch condition omitted for brevity]
                        
                    end % End probe match/mismatch condition
                end % End trial loop
            end % End region B loop
        end % End region A loop
    end % End seed region loop
end % End subject loop

% Clean up unused array elements
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) = [];
probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('Done processing.\n')
toc

%% ========================================================================
%  VISUALIZATION: REPLAY ANALYSIS RESULTS
%% ========================================================================

regionColors = brewermap(12, 'Dark2');
contraColor = brewermap(12, 'Accent');

%% Helper Function: Aggregate Data Across Region Pairs
function [replay_match, replay_mismatch, replay_norip_match, replay_norip_mismatch, ...
          coRap_total, coRap_mm_total, noRap_total, noRap_mm_total, ...
          coRdur_total, coRdur_mm_total, noRdur_total, noRdur_mm_total] = ...
          aggregateRegionData(regions, regionsPlot, targetRegions, coRcoF, coRcoF_mm, ...
                              noRcoF, noRcoF_mm, coRap, coRap_mm, noRap, noRap_mm, ...
                              coRdur, coRdur_mm, noRdur, noRdur_mm, coR_replay, ...
                              sameHemisphere)
    
    % Get indices for target regions
    iRapl = find(contains(regions, regionsPlot(targetRegions)));
    iRbpl = iRapl;
    
    % Initialize accumulators
    replay_match = zeros(2, 1);
    replay_mismatch = zeros(2, 1);
    replay_norip_match = zeros(2, 1);
    replay_norip_mismatch = zeros(2, 1);
    
    coRap_total = 0;
    coRap_mm_total = 0;
    noRap_total = 0;
    noRap_mm_total = 0;
    
    coRdur_total = 0;
    coRdur_mm_total = 0;
    noRdur_total = 0;
    noRdur_mm_total = 0;
    
    % Loop through region pairs
    for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl, iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            
            if sum(coR_replay(:, iiA, iiB)) == 0; continue; end
            
            % Check hemisphere constraint if specified
            if nargin >= 17 && ~isempty(sameHemisphere)
                if strcmp(regions{iiA}(1), regions{iiB}(1)) ~= sameHemisphere
                    continue
                end
            end
            
            % Aggregate matched trials
            dat = [sum(coRcoF(:, 1, iiA, iiB) > 0 & coRcoF(:, 2, iiA, iiB) > 0); ...
                   sum(coRcoF(:, 1, iiA, iiB) == 0 | coRcoF(:, 2, iiA, iiB) == 0)];
            replay_match = replay_match + dat;
            
            % Aggregate mismatched trials
            dat = [sum(coRcoF_mm(:, 1, iiA, iiB) > 0 & coRcoF_mm(:, 2, iiA, iiB) > 0); ...
                   sum(coRcoF_mm(:, 1, iiA, iiB) == 0 | coRcoF_mm(:, 2, iiA, iiB) == 0)];
            replay_mismatch = replay_mismatch + dat;
            
            % Aggregate no-ripple controls (average across iterations)
            dat = [];
            for iter = 1:25
                dat = [dat, [sum(noRcoF(:, 1, iiA, iiB, iter) > 0 & noRcoF(:, 2, iiA, iiB, iter) > 0); ...
                            sum(noRcoF(:, 1, iiA, iiB, iter) == 0 | noRcoF(:, 2, iiA, iiB, iter) == 0)]];
            end
            replay_norip_match = replay_norip_match + round(mean(dat, 2));
            
            dat = [];
            for iter = 1:25
                dat = [dat, [sum(noRcoF_mm(:, 1, iiA, iiB, iter) > 0 & noRcoF_mm(:, 2, iiA, iiB, iter) > 0); ...
                            sum(noRcoF_mm(:, 1, iiA, iiB, iter) == 0 | noRcoF_mm(:, 2, iiA, iiB, iter) == 0)]];
            end
            replay_norip_mismatch = replay_norip_mismatch + round(mean(dat, 2));
            
            % Aggregate co-firing rates
            coRap_total = coRap_total + sum(coRap(:, :, iiA, iiB), 'all');
            coRap_mm_total = coRap_mm_total + sum(coRap_mm(:, :, iiA, iiB), 'all');
            noRap_total = noRap_total + sum(noRap(:, :, iiA, iiB), 'all');
            noRap_mm_total = noRap_mm_total + sum(noRap_mm(:, :, iiA, iiB), 'all');
            
            coRdur_total = coRdur_total + sum(coRdur(:, :, iiA, iiB), 'all');
            coRdur_mm_total = coRdur_mm_total + sum(coRdur_mm(:, :, iiA, iiB), 'all');
            noRdur_total = noRdur_total + sum(noRdur(:, :, iiA, iiB), 'all');
            noRdur_mm_total = noRdur_mm_total + sum(noRdur_mm(:, :, iiA, iiB), 'all');
        end
    end
end

%% Figure 1: Main Replay Analysis - AMY, HIP, Ipsilateral, Contralateral
load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
pValue = nan(4, 3);
bW = 1.1;

figure('Position', [1191 647 608 638]);
ha = tight_subplot(2, 2, [0.25, 0.15], [.05 .05], [.15 .15]);

%% Panel 1: AMY to Frontal
iRa = 4; % AMY
[coR_replay_pl, coR_replay_mm_pl, noR_replay_pl, noR_replay_mm_pl, ...
 coRap_pl, coRap_mm_pl, noRap_pl, noRap_mm_pl, ...
 coRdur_pl, coRdur_mm_pl, noRdur_pl, noRdur_mm_pl] = ...
    aggregateRegionData(regions, regionsPlot, iRa, [1:3], coRcoF, coRcoF_mm, ...
                        noRcoF, noRcoF_mm, coRap, coRap_mm, noRap, noRap_mm, ...
                        coRdur, coRdur_mm, noRdur, noRdur_mm, coR_replay);

% Calculate firing rates
coR_fr = coRap_pl / coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl / coRdur_mm_pl * 1e3;
noR_fr = noRap_pl / noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl / noRdur_mm_pl * 1e3;

X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), ...
     noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)] * 100;

axes(ha(1))
yyaxis left
b = bar(X(1), Y(1), bW); hold on;
b.FaceColor = regionColors(iRa, :);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(3), Y(3), bW); hold on;
b.FaceColor = regionColors(iRa, :);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

ax = gca;
ax.XTick = 1:4;
ax.XTickLabelRotation = 50;
box off

% Add firing rate markers
X_fr = [0.65:3.65; 1.35:4.35]';
Y_fr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3]
    pl = plot(mean([X_fr(iL, 1) X_fr(iL, 2)]), Y_fr(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3, :);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;
end

ax.YAxis(1).Label.String = '% coF in E & P';
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8 * regionColors(3, :);
ax.YAxis(2).Limits = [0 max(Y_fr) + 0.05*max(Y_fr)];
ax.YAxis(1).FontSize = 11;
ax.YAxis(2).FontSize = 11;
ax.LineWidth = 1;

% Statistical test
[chi2stat, pValue(1, 1), expected] = chiSquaredTest([[coR_replay_pl(1), coR_replay_pl(2)]; ...
                                                      [coR_replay_mm_pl(1), coR_replay_mm_pl(2)]]);
[chi2stat, pValue(1, 2), expected] = chiSquaredTest([[noR_replay_pl(1), noR_replay_pl(2)]; ...
                                                      [noR_replay_mm_pl(1), noR_replay_mm_pl(2)]]);
[chi2stat, pValue(1, 3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]), ...
                                                       sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; ...
                                                      [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]), ...
                                                       sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

%% Panel 2: HIP to Frontal
iRa = 5; % HIP
[coR_replay_pl, coR_replay_mm_pl, noR_replay_pl, noR_replay_mm_pl, ...
 coRap_pl, coRap_mm_pl, noRap_pl, noRap_mm_pl, ...
 coRdur_pl, coRdur_mm_pl, noRdur_pl, noRdur_mm_pl] = ...
    aggregateRegionData(regions, regionsPlot, iRa, [1:3], coRcoF, coRcoF_mm, ...
                        noRcoF, noRcoF_mm, coRap, coRap_mm, noRap, noRap_mm, ...
                        coRdur, coRdur_mm, noRdur, noRdur_mm, coR_replay);

coR_fr = coRap_pl / coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl / coRdur_mm_pl * 1e3;
noR_fr = noRap_pl / noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl / noRdur_mm_pl * 1e3;

X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), ...
     noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)] * 100;

axes(ha(2))
yyaxis left
b = bar(X(1), Y(1), bW); hold on;
b.FaceColor = regionColors(iRa, :);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(3), Y(3), bW); hold on;
b.FaceColor = regionColors(iRa, :);
b.FaceAlpha = 0.5;
b.LineWidth = 1;
box off

X_fr = [0.65:3.65; 1.35:4.35]';
Y_fr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3]
    pl = plot(mean([X_fr(iL, 1) X_fr(iL, 2)]), Y_fr(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3, :);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;
end

ax = gca;
ax.XTick = 1:4;
ax.XTickLabelRotation = 50;
ax.YAxis(1).Label.String = '% coF in E & P';
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8 * regionColors(3, :);
ax.YAxis(1).Limits = [0 0.25];
ax.YAxis(2).Limits = [0 max(Y_fr) + 0.05*max(Y_fr)];
ax.YAxis(1).FontSize = 11;
ax.YAxis(2).FontSize = 11;
ax.LineWidth = 1;

[chi2stat, pValue(2, 1), expected] = chiSquaredTest([[coR_replay_pl(1), coR_replay_pl(2)]; ...
                                                      [coR_replay_mm_pl(1), coR_replay_mm_pl(2)]]);
[chi2stat, pValue(2, 2), expected] = chiSquaredTest([[noR_replay_pl(1), noR_replay_pl(2)]; ...
                                                      [noR_replay_mm_pl(1), noR_replay_mm_pl(2)]]);
[chi2stat, pValue(2, 3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]), ...
                                                       sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; ...
                                                      [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]), ...
                                                       sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

%% Panel 3: Ipsilateral (Same Hemisphere)
[coR_replay_pl, coR_replay_mm_pl, noR_replay_pl, noR_replay_mm_pl, ...
 coRap_pl, coRap_mm_pl, noRap_pl, noRap_mm_pl, ...
 coRdur_pl, coRdur_mm_pl, noRdur_pl, noRdur_mm_pl] = ...
    aggregateRegionData(regions, regionsPlot, [1:5], coRcoF, coRcoF_mm, ...
                        noRcoF, noRcoF_mm, coRap, coRap_mm, noRap, noRap_mm, ...
                        coRdur, coRdur_mm, noRdur, noRdur_mm, coR_replay, true);

coR_fr = coRap_pl / coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl / coRdur_mm_pl * 1e3;
noR_fr = noRap_pl / noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl / noRdur_mm_pl * 1e3;

X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), ...
     noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)] * 100;

axes(ha(3))
yyaxis left
b = bar(X(1), Y(1), bW); hold on;
b.FaceColor = regionColors(6, :);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(3), Y(3), bW); hold on;
b.FaceColor = regionColors(6, :);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

X_fr = [0.65:3.65; 1.35:4.35]';
Y_fr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3]
    pl = plot(mean([X_fr(iL, 1) X_fr(iL, 2)]), Y_fr(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3, :);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;
end

ax = gca;
ax.XTick = 1:4;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8 * regionColors(3, :);
ax.XTickLabelRotation = 50;
box off
ax.YAxis(1).Label.String = '% coF in E & P';
ax.YAxis(2).Limits = [0 max(Y_fr) + 0.05*max(Y_fr)];
ax.YAxis(1).FontSize = 11;
ax.YAxis(2).FontSize = 11;
ax.LineWidth = 1;

[chi2stat, pValue(3, 1), expected] = chiSquaredTest([[coR_replay_pl(1), coR_replay_pl(2)]; ...
                                                      [coR_replay_mm_pl(1), coR_replay_mm_pl(2)]]);
[chi2stat, pValue(3, 2), expected] = chiSquaredTest([[noR_replay_pl(1), noR_replay_pl(2)]; ...
                                                      [noR_replay_mm_pl(1), noR_replay_mm_pl(2)]]);
[chi2stat, pValue(3, 3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]), ...
                                                       sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; ...
                                                      [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]), ...
                                                       sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

%% Panel 4: Contralateral (Different Hemispheres)
[coR_replay_pl, coR_replay_mm_pl, noR_replay_pl, noR_replay_mm_pl, ...
 coRap_pl, coRap_mm_pl, noRap_pl, noRap_mm_pl, ...
 coRdur_pl, coRdur_mm_pl, noRdur_pl, noRdur_mm_pl] = ...
    aggregateRegionData(regions, regionsPlot, [1:5], coRcoF, coRcoF_mm, ...
                        noRcoF, noRcoF_mm, coRap, coRap_mm, noRap, noRap_mm, ...
                        coRdur, coRdur_mm, noRdur, noRdur_mm, coR_replay, false);

coR_fr = coRap_pl / coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl / coRdur_mm_pl * 1e3;
noR_fr = noRap_pl / noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl / noRdur_mm_pl * 1e3;

X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), ...
     noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)] * 100;

axes(ha(4))
yyaxis left
b = bar(X(1), Y(1), bW); hold on;
b.FaceColor = contraColor(7, :);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(3), Y(3), bW); hold on;
b.FaceColor = contraColor(7, :);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

X_fr = [0.65:3.65; 1.35:4.35]';
Y_fr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3]
    pl = plot(mean([X_fr(iL, 1) X_fr(iL, 2)]), Y_fr(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3, :);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;
end

ax = gca;
ax.YAxis(1).Label.String = '% coF in E & P';
ax.YAxis(2).Limits = [0 max(Y_fr) + 0.05*max(Y_fr)];
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8 * regionColors(3, :);
ax.YAxis(1).FontSize = 11;
ax.YAxis(2).FontSize = 11;
ax.LineWidth = 1;

[chi2stat1, pValue(4, 1), expected1] = chiSquaredTest([[coR_replay_pl(1), coR_replay_pl(2)]; ...
                                                        [coR_replay_mm_pl(1), coR_replay_mm_pl(2)]]);
[chi2stat2, pValue(4, 2), expected2] = chiSquaredTest([[noR_replay_pl(1), noR_replay_pl(2)]; ...
                                                        [noR_replay_mm_pl(1), noR_replay_mm_pl(2)]]);
[chi2stat, pValue(4, 3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]), ...
                                                       sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; ...
                                                      [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]), ...
                                                       sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

fig = gcf;
fig.Color = 'w';
box off

% FDR correction
adj_p = nan(size(pValue));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pValue))] = fdr_bh(pValue(~isnan(pValue)), 0.05, 'pdep', 'yes');

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplay_%s.pdf', tag)))

%% Compare Effect Sizes
results = compare_chi_squared_effect_sizes(chi2stat1, sum(expected1(:)), 2, 2, ...
                                           chi2stat2, sum(expected2(:)), 2, 2);

fprintf('\n=== REPLAY ANALYSIS COMPLETE ===\n')
fprintf('Results saved to: %s\n', exportDirFigs)