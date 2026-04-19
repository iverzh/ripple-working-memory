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
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/ripple-working-memory/code'))


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
subjID = nan(5e5, length(regions),length(regions)); %subject / session identifier
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
                                subjID(ii,iRa, iRb) = subj;
                                
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



fprintf('Done processing.\n')
toc


%% ========================================================================
%  VISUALIZATION: REPLAY ANALYSIS RESULTS
%% ========================================================================

regionColors = brewermap(12, 'Dark2');
contraColor = brewermap(12, 'Accent');


%% Figure 1: Main Replay Analysis - AMY, HIP, Ipsilateral, Contralateral
load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
pValue = nan(4, 3);
bW = 1.1;

figure('Position', [1191 647 608 638]);
ha = tight_subplot(2, 2, [0.25, 0.15], [.05 .05], [.15 .15]);

% Shared input for aggregateRegionData (panel-specific fields set below)
in_base = struct();
in_base.regions     = regions;     % region labels (cell array, per channel)
in_base.regionsPlot = regionsPlot; % region groups to query against `regions`
in_base.coRcoF      = coRcoF;      % co-ripple co-firing counts, matched trials
in_base.coRcoF_mm   = coRcoF_mm;   % co-ripple co-firing counts, mismatched trials
in_base.noRcoF      = noRcoF;      % no-ripple control co-firing counts, matched
in_base.noRcoF_mm   = noRcoF_mm;   % no-ripple control co-firing counts, mismatched
in_base.coRap       = coRap;       % co-ripple action-potential counts, matched
in_base.coRap_mm    = coRap_mm;    % co-ripple AP counts, mismatched
in_base.noRap       = noRap;       % no-ripple control AP counts, matched
in_base.noRap_mm    = noRap_mm;    % no-ripple control AP counts, mismatched
in_base.coRdur      = coRdur;      % co-ripple durations (ms), matched
in_base.coRdur_mm   = coRdur_mm;   % co-ripple durations, mismatched
in_base.noRdur      = noRdur;      % no-ripple control durations, matched
in_base.noRdur_mm   = noRdur_mm;   % no-ripple control durations, mismatched
in_base.coR_replay  = coR_replay;  % co-ripple replay mask per region pair

% Panel 1: AMY to Frontal
iRa = 4; % AMY
in = in_base;
in.targetRegionsA = iRa;       % AMY
in.targetRegionsB = [1:3];     % Frontal group
in.respTimes      = respTimes; % trial response times for fast/slow split
in.loads          = loads;     % memory load per trial (for load-filtered median)
out = aggregateRegionData(in);

% Calculate firing rates
coR_fr    = out.coRap_total    / out.coRdur_total    * 1e3;
coR_mm_fr = out.coRap_mm_total / out.coRdur_mm_total * 1e3;
noR_fr    = out.noRap_total    / out.noRdur_total    * 1e3;
noR_mm_fr = out.noRap_mm_total / out.noRdur_mm_total * 1e3;

X = [1:4];
Y = [out.replay_match(1)/sum(out.replay_match), out.replay_mismatch(1)/sum(out.replay_mismatch), ...
     out.replay_norip_match(1)/sum(out.replay_norip_match), out.replay_norip_mismatch(1)/sum(out.replay_norip_mismatch)] * 100;

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
[chi2stat, pValue(1, 1), expected] = chiSquaredTest([[out.replay_match(1),       out.replay_match(2)]; ...
                                                      [out.replay_mismatch(1),    out.replay_mismatch(2)]]);
[chi2stat, pValue(1, 2), expected] = chiSquaredTest([[out.replay_norip_match(1), out.replay_norip_match(2)]; ...
                                                      [out.replay_norip_mismatch(1), out.replay_norip_mismatch(2)]]);
[chi2stat, pValue(1, 3), expected] = chiSquaredTest([[sum([out.replay_match(1), out.replay_mismatch(1)]), ...
                                                       sum([out.replay_match(2), out.replay_mismatch(2)])]; ...
                                                      [sum([out.replay_norip_match(1), out.replay_norip_mismatch(1)]), ...
                                                       sum([out.replay_norip_match(2), out.replay_norip_mismatch(2)])]]);

% Panel 2: HIP to Frontal
iRa = 5; % HIP
in = in_base;
in.targetRegionsA = iRa;   % HIP
in.targetRegionsB = [1:3]; % Frontal group
out = aggregateRegionData(in);

coR_fr    = out.coRap_total    / out.coRdur_total    * 1e3;
coR_mm_fr = out.coRap_mm_total / out.coRdur_mm_total * 1e3;
noR_fr    = out.noRap_total    / out.noRdur_total    * 1e3;
noR_mm_fr = out.noRap_mm_total / out.noRdur_mm_total * 1e3;

X = [1:4];
Y = [out.replay_match(1)/sum(out.replay_match), out.replay_mismatch(1)/sum(out.replay_mismatch), ...
     out.replay_norip_match(1)/sum(out.replay_norip_match), out.replay_norip_mismatch(1)/sum(out.replay_norip_mismatch)] * 100;

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

[chi2stat, pValue(2, 1), expected] = chiSquaredTest([[out.replay_match(1),       out.replay_match(2)]; ...
                                                      [out.replay_mismatch(1),    out.replay_mismatch(2)]]);
[chi2stat, pValue(2, 2), expected] = chiSquaredTest([[out.replay_norip_match(1), out.replay_norip_match(2)]; ...
                                                      [out.replay_norip_mismatch(1), out.replay_norip_mismatch(2)]]);
[chi2stat, pValue(2, 3), expected] = chiSquaredTest([[sum([out.replay_match(1), out.replay_mismatch(1)]), ...
                                                       sum([out.replay_match(2), out.replay_mismatch(2)])]; ...
                                                      [sum([out.replay_norip_match(1), out.replay_norip_mismatch(1)]), ...
                                                       sum([out.replay_norip_match(2), out.replay_norip_mismatch(2)])]]);

% Panel 3: Ipsilateral (Same Hemisphere)
in = in_base;
in.targetRegionsA  = [1:5]; % all region groups, A side
in.targetRegionsB  = [1:5]; % all region groups, B side
in.sameHemisphere  = true;  % restrict to pairs in the same hemisphere
out = aggregateRegionData(in);

coR_fr    = out.coRap_total    / out.coRdur_total    * 1e3;
coR_mm_fr = out.coRap_mm_total / out.coRdur_mm_total * 1e3;
noR_fr    = out.noRap_total    / out.noRdur_total    * 1e3;
noR_mm_fr = out.noRap_mm_total / out.noRdur_mm_total * 1e3;

X = [1:4];
Y = [out.replay_match(1)/sum(out.replay_match), out.replay_mismatch(1)/sum(out.replay_mismatch), ...
     out.replay_norip_match(1)/sum(out.replay_norip_match), out.replay_norip_mismatch(1)/sum(out.replay_norip_mismatch)] * 100;

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

[chi2stat, pValue(3, 1), expected] = chiSquaredTest([[out.replay_match(1),       out.replay_match(2)]; ...
                                                      [out.replay_mismatch(1),    out.replay_mismatch(2)]]);
[chi2stat, pValue(3, 2), expected] = chiSquaredTest([[out.replay_norip_match(1), out.replay_norip_match(2)]; ...
                                                      [out.replay_norip_mismatch(1), out.replay_norip_mismatch(2)]]);
[chi2stat, pValue(3, 3), expected] = chiSquaredTest([[sum([out.replay_match(1), out.replay_mismatch(1)]), ...
                                                       sum([out.replay_match(2), out.replay_mismatch(2)])]; ...
                                                      [sum([out.replay_norip_match(1), out.replay_norip_mismatch(1)]), ...
                                                       sum([out.replay_norip_match(2), out.replay_norip_mismatch(2)])]]);

% Panel 4: Contralateral (Different Hemispheres)
in = in_base;
in.targetRegionsA  = [1:5]; % all region groups, A side
in.targetRegionsB  = [1:5]; % all region groups, B side
in.sameHemisphere  = false; % restrict to pairs in opposite hemispheres
out = aggregateRegionData(in);

coR_fr    = out.coRap_total    / out.coRdur_total    * 1e3;
coR_mm_fr = out.coRap_mm_total / out.coRdur_mm_total * 1e3;
noR_fr    = out.noRap_total    / out.noRdur_total    * 1e3;
noR_mm_fr = out.noRap_mm_total / out.noRdur_mm_total * 1e3;

X = [1:4];
Y = [out.replay_match(1)/sum(out.replay_match), out.replay_mismatch(1)/sum(out.replay_mismatch), ...
     out.replay_norip_match(1)/sum(out.replay_norip_match), out.replay_norip_mismatch(1)/sum(out.replay_norip_mismatch)] * 100;

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

[chi2stat1, pValue(4, 1), expected1] = chiSquaredTest([[out.replay_match(1),       out.replay_match(2)]; ...
                                                        [out.replay_mismatch(1),    out.replay_mismatch(2)]]);
[chi2stat2, pValue(4, 2), expected2] = chiSquaredTest([[out.replay_norip_match(1), out.replay_norip_match(2)]; ...
                                                        [out.replay_norip_mismatch(1), out.replay_norip_mismatch(2)]]);
[chi2stat, pValue(4, 3), expected] = chiSquaredTest([[sum([out.replay_match(1), out.replay_mismatch(1)]), ...
                                                       sum([out.replay_match(2), out.replay_mismatch(2)])]; ...
                                                      [sum([out.replay_norip_match(1), out.replay_norip_mismatch(1)]), ...
                                                       sum([out.replay_norip_match(2), out.replay_norip_mismatch(2)])]]);

fig = gcf;
fig.Color = 'w';
box off

% FDR correction
adj_p = nan(size(pValue));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pValue))] = fdr_bh(pValue(~isnan(pValue)), 0.05, 'pdep', 'yes');

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplay_%s.pdf', tag)))

%% Figure 2: Replay by subject
in = struct();
in.regions        = regions;        % region labels (cell array, per channel)
in.regionsPlot    = regionsPlot;    % region groups to query against `regions`
in.targetRegionsA = [1:5];          % all region groups, A side
in.targetRegionsB = [1:5];          % all region groups, B side
in.coRcoF         = coRcoF;         % co-ripple co-firing counts, matched trials
in.noRcoF         = noRcoF;         % no-ripple control co-firing counts, matched
in.coRap          = coRap;          % co-ripple action-potential counts, matched
in.coRap_mm       = coRap_mm;       % co-ripple AP counts, mismatched
in.noRap          = noRap;          % no-ripple control AP counts, matched
in.noRap_mm       = noRap_mm;       % no-ripple control AP counts, mismatched
in.coRdur         = coRdur;         % co-ripple durations (ms), matched
in.coRdur_mm      = coRdur_mm;      % co-ripple durations, mismatched
in.noRdur         = noRdur;         % no-ripple control durations, matched
in.noRdur_mm      = noRdur_mm;      % no-ripple control durations, mismatched
in.subjID         = subjID;         % subject index per trial and region pair
in.respTimes      = respTimes;      % trial response times for fast/slow split
in.loads          = loads;          % memory load per trial

out = aggregateRegionDataSubj(in);

coSubj_fn     = out.coSubj;
noSubj_fn     = out.noSubj;
coFastSubj_fn = out.coFastSubj;
coSlowSubj_fn = out.coSlowSubj;
nSubj         = numel(coSubj_fn);

% Fix NaN edge cases: if one RT bin is NaN but the other isn't, set to 0
coFastSubj_fn(isnan(coFastSubj_fn) & ~isnan(coSlowSubj_fn)) = 0;
coSlowSubj_fn(~isnan(coFastSubj_fn) & isnan(coSlowSubj_fn)) = 0;

% Statistics
p_coNo   = pairedPermutationTest(coSubj_fn(~isnan(coSubj_fn))', noSubj_fn(~isnan(noSubj_fn))', 10e3);
p_rtSplit = pairedPermutationTest(coFastSubj_fn(~isnan(coFastSubj_fn))', coSlowSubj_fn(~isnan(coSlowSubj_fn))', 10e3);
[P_rt, observeddifference_rt, effectsize_rt] = permutationTestMed(coFastSubj_fn(~isnan(coFastSubj_fn))', coSlowSubj_fn(~isnan(coSlowSubj_fn))', 10e3);
[P_rt_sr, H_rt_sr] = signrank(coFastSubj_fn, coSlowSubj_fn, 'tail', 'both');

clr = brewermap(10, 'Paired');

figure('Position', [738 907 201 248]);

% --- Subplot 1: co-ripple vs no-ripple replay rate per subject ---
subplot(1, 2, 1)
for iS = 1:nSubj
    jitterX = 0.5 * (rand - 0.5);

    pl = plot([1+jitterX 2+jitterX], [coSubj_fn(iS), noSubj_fn(iS)], '-'); hold on;
    pl.Color = [0.8 0.8 0.8];
    pl.LineWidth = 0.25;

    pl = scatter(1+jitterX, coSubj_fn(iS), 'o'); hold on;
    pl.SizeData = 20;
    pl.MarkerFaceColor = clr(8, :);
    pl.MarkerEdgeAlpha = 0.0;

    pl = scatter(2+jitterX, noSubj_fn(iS), 'o'); hold on;
    pl.SizeData = 20;
    pl.MarkerFaceColor = clr(1, :);
    pl.MarkerEdgeAlpha = 0.0;
end

boxplot(coSubj_fn, 'positions', 1, 'width', 0.2, 'Colors', [0 0 0], 'Symbol', '', 'Whisker', 1);
boxplot(noSubj_fn, 'positions', 2, 'width', 0.2, 'Colors', [0 0 0], 'Symbol', '', 'Whisker', 1);

xlim([0 3])
ax = gca;
ax.XTick = 1:2;
box off

% --- Subplot 2: fast RT vs slow RT co-ripple replay rate per subject ---
subplot(1, 2, 2)
for iS = 1:nSubj
    jitterX = 0.5 * (rand - 0.5);

    pl = plot([1+jitterX 2+jitterX], [coFastSubj_fn(iS), coSlowSubj_fn(iS)], '-'); hold on;
    pl.Color = [0.8 0.8 0.8];
    pl.LineWidth = 0.25;

    pl = scatter(1+jitterX, coFastSubj_fn(iS), 'o'); hold on;
    pl.SizeData = 20;
    pl.MarkerFaceColor = clr(8, :);
    pl.MarkerEdgeAlpha = 0.0;

    pl = scatter(2+jitterX, coSlowSubj_fn(iS), 'o'); hold on;
    pl.SizeData = 20;
    pl.MarkerFaceColor = 0.6 * clr(8, :);
    pl.MarkerEdgeAlpha = 0.0;
end

boxplot(coFastSubj_fn, 'positions', 1, 'width', 0.2, 'Colors', [0 0 0], 'Symbol', '', 'Whisker', 1);
boxplot(coSlowSubj_fn, 'positions', 2, 'width', 0.2, 'Colors', [0 0 0], 'Symbol', '', 'Whisker', 1);

ax = gca;
ax.XTick = 1:2;
xlim([0 3])
box off

fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirFigs, 'ReplayBySubject_fn.pdf'))


%% Figure 3: Fast v slow v shuffle v no ripple

regionColors = brewermap(12, 'Dark2');
contraColor  = brewermap(12, 'Accent');
bW = 0.8;

% --- Shared base struct (data arrays common to both functions) ---
in_base = struct();
in_base.regions     = regions;       % region labels (cell array, per channel)
in_base.regionsPlot = regionsPlot;   % region groups to query against `regions`
in_base.coRcoF      = coRcoF;        % co-ripple co-firing counts, matched trials
in_base.coRcoF_mm   = coRcoF_mm;     % co-ripple co-firing counts, mismatched trials
in_base.noRcoF      = noRcoF;        % no-ripple control co-firing counts, matched
in_base.noRcoF_mm   = noRcoF_mm;     % no-ripple control co-firing counts, mismatched
in_base.coRap       = coRap;         % co-ripple action-potential counts, matched
in_base.coRap_mm    = coRap_mm;      % co-ripple AP counts, mismatched
in_base.noRap       = noRap;         % no-ripple control AP counts, matched
in_base.noRap_mm    = noRap_mm;      % no-ripple control AP counts, mismatched
in_base.coRdur      = coRdur;        % co-ripple durations (ms), matched
in_base.coRdur_mm   = coRdur_mm;     % co-ripple durations, mismatched
in_base.noRdur      = noRdur;        % no-ripple control durations, matched
in_base.noRdur_mm   = noRdur_mm;     % no-ripple control durations, mismatched
in_base.coR_replay  = coR_replay;    % co-ripple replay mask per region pair
in_base.respTimes   = respTimes;     % trial response times for fast/slow RT split
in_base.loads       = loads;         % memory load per trial

% Additional fields needed only by the shuffle function
shuff_base = in_base;
shuff_base.subjID = subjID;          % subject index per trial/pair (for within-subject shuffle)
shuff_base.nIter  = 2e2;            % number of shuffle iterations


% ---- Panel definitions -----------------------------------------------
% Each row: {label, targetRegionsA, targetRegionsB, loadParam, sameHemisphere, barColor}
panels = { ...
    'AMY ctx', 4,     [1:3], 3, [],    regionColors(4,:); ...
    'HIP ctx', 5,     [1:3], 3, [],    regionColors(5,:); ...
    'Ipsilateral', [1:5], [1:5], 3, true,  regionColors(6,:); ...
    'Contralateral',[1:5],[1:5], 3, false, contraColor(7,:)   ...
};

figure('Position', [2 632 755 299]);
ha = tight_subplot(1, 4, [0.25, 0.05], [.05 .05], [.075 .05]);

for iPanel = 1:4

    label         = panels{iPanel, 1};
    trgA          = panels{iPanel, 2};
    trgB          = panels{iPanel, 3};
    lParam        = panels{iPanel, 4};
    sameHemi      = panels{iPanel, 5};
    barColor      = panels{iPanel, 6};

    % -- aggregateRegionData: fast/slow RT split and no-ripple control --
    aggIn                = in_base;
    aggIn.targetRegionsA = trgA;
    aggIn.targetRegionsB = trgB;
    aggIn.loadParam      = lParam;
    aggIn.sameHemisphere = sameHemi;
    aggOut = aggregateRegionData(aggIn);

    fastRate  = aggOut.replay_fast_match(1)     / sum(aggOut.replay_fast_match)     * 100;
    slowRate  = aggOut.replay_slow_match(1)     / sum(aggOut.replay_slow_match)     * 100;
    noRipRate = aggOut.replay_norip_match(1)    / sum(aggOut.replay_norip_match)    * 100;

    % -- aggregateRegionDataShuff: shuffle null distribution --
    shuffIn                = shuff_base;
    shuffIn.targetRegionsA = trgA;
    shuffIn.targetRegionsB = trgB;
    shuffIn.loadParam      = lParam;
    shuffIn.sameHemisphere = sameHemi;
    shuffOut = aggregateRegionDataShuff(shuffIn);

    % -- Plot --
    X = [1, 2, 3, 3.5];
    Y = [fastRate, NaN, slowRate, mean(shuffOut.shuff)];

    axes(ha(iPanel))

    % Fast RT bar (solid)
    b = bar(X(1), Y(1), bW); hold on;
    b.FaceColor = barColor;
    b.FaceAlpha = 1;
    b.LineWidth = 1;

    % Slow RT bar (striped)
    b = bar_striped(X(2), Y(3), bW, 0.02, 5); hold on;
    b.FaceColor = barColor;
    b.FaceAlpha = 1;
    b.LineWidth = 1;

    % Shuffle mean bar (translucent)
    b = bar(X(4), Y(4), bW); hold on;
    b.FaceColor = barColor;
    b.FaceAlpha = 0.5;
    b.LineWidth = 1;

    % Shuffle distribution ? filled polygon to the left of X(4)
    [N, bn] = histcounts(shuffOut.shuff, 0:0.005:0.40, 'Normalization', 'probability');
    bn = movmean(bn, 2);
    bn(1) = [];
    N  = smoothdata(N, 'gaussian', 5);
    N  = N / max(N);
    xx = (X(1) - bW/2 - 0.2) - N(N > 0);
    xx(end+1) = xx(1);
    yy = bn(N > 0);
    yy(end+1) = yy(1);

    fl = fill(xx, yy, 'r'); hold on;
    fl.LineWidth  = 1.5;
    fl.EdgeColor  = barColor;
    fl.FaceColor  = [0.7 0.7 0.7];

    % Axes formatting
    ax = gca;
    ax.XTick             = [-0.1, 1:4];
    ax.XTickLabelRotation = 50;
    ax.YAxis(1).Color    = [0 0 0];
    ax.YAxis(1).Limits   = [0, max([Y(1), Y(3), Y(4)]) * 1.1];
    ax.YAxis(1).FontSize = 11;
    ax.LineWidth         = 1;
    xlim([-1 4.5])
    box off
    title(label, 'FontSize', 9)

    % -- Statistics --
    [~, pFastSlow] = chiSquaredTest([ ...
        [aggOut.replay_fast_match(1), aggOut.replay_fast_match(2)]; ...
        [aggOut.replay_slow_match(1), aggOut.replay_slow_match(2)]]);
    fprintf('%s  fast vs slow RT: p = %.5f\n', label, pFastSlow)

    coTot = aggOut.replay_fast_match + aggOut.replay_slow_match;
    [~, pCoNo] = chiSquaredTest([ ...
        [coTot(1),                       coTot(2)]; ...
        [aggOut.replay_norip_match(1),   aggOut.replay_norip_match(2)]]);
    fprintf('%s  co vs no-ripple: p = %.5f\n', label, pCoNo)

end

fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirFigs, 'ReplayFastSlow_shuff_fn.pdf'))

