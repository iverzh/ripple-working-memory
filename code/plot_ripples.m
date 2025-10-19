%% Ripple Characteristics Analysis
% This script analyzes ripple oscillations in human intracranial EEG data
% during a working memory (Sternberg) task.
% It computes ripple metrics, co-occurrence rates across brain regions,
% and generates visualization figures for publication.

close all 
clc
clear

%% Setup paths
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/ripple-detection/code/CCG_tools'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))

%% Define directories
dataDirectoryOrig = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';
channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';

% Create output directories if they don't exist
if ~isfolder(exportDir); mkdir(exportDir); end
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

%% Load file lists
flst = dir(fullfile(dataDirectoryOrig, '*LFP_micro*'));
subj_list_full = strrep({flst.name}', '_LFP_micro.mat', '');

unitfiles = {dir(fullfile(dataDirectoryOrig, '*unit*')).name}';
LFPfilesMacr = {dir(fullfile(dataDirectoryOrig, '*macro*')).name}';
LFPfilesMicr = {dir(fullfile(dataDirectoryOrig, '*micro*')).name}';
taskfiles = {dir(fullfile(dataDirectoryOrig, '../task','*task*')).name}';

%% Analysis parameters
recordingState = 'wake';
location = 'NC';
tag = 'vHG';
sfreq = 1e3;
modifier = '1kHz_template_z25';

% Brain regions
regionsAll = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
              'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'};

% Cross-correlogram parameters
win = 1000;  % Window size (ms)
binSize = 1;  % Bin size (ms)
hist_bins = 2 * win / binSize;
nIter = 200;
FDRwin = 500;  % +/- window (ms) for FDR significance

%% Preallocate variables
countChPerRegion = zeros(1, length(regions));
countUnitPerRegion = zeros(1, length(regions));

% Ripple LFP data by region
ripple_LFP_region_all = cell(length(regions), 1);
for iRa = 1:length(regions)
    ripple_LFP_region_all{iRa} = nan(1e6, length(-500:500));       
end

% Region counters
cRegion = ones(1, length(regions));
cchRegion = ones(1, length(regions));
cUnitRegion = ones(1, length(regions));

% Trial data
subjID = nan(6.5e3, 1);
correct_trials_all = nan(6.5e3, 1);
load_trials_all = nan(6.5e3, 1);
probe_in_out = nan(6.5e3, 1);
respLatency = nan(6.5e3, 1);

% Ripple metrics
densityAll = [];
freqAll = [];
durAll = [];
ampAll = [];
locAll = [];
subjAll = [];

% Unit properties
burstIndex = [];
troughToPeak = [];
FR = [];
uType = [];
uSubj = [];
uID = [];
uChanAll = [];

% ERSP data
ERSPall = [];
cnt = 1;

% Co-ripple data
cTrial = 1;
coRnum = nan(1e5, length(regionsAll), length(regionsAll));
coRnumShuff = nan(1e5, length(regionsAll), length(regionsAll));
coRdur = nan(1e5, length(regionsAll), length(regionsAll));
coR_PRTH = zeros(hist_bins-1, 2, 1e5, length(regionsAll), length(regionsAll));
coR_PRTH_null = zeros(hist_bins-1, 1e5, length(regionsAll), length(regionsAll));
ccgNrip = zeros(2, 1e5, length(regionsAll), length(regionsAll));
ccgCount = ones(length(regionsAll));

%% Main processing loop - iterate through subjects
for subj = 1:length(subj_list_full)
    tic
    subject = subj_list_full{subj};
    fprintf('Processing %s ... \n', subject)
    
    %% Load micro LFP data and channel labels
    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectoryOrig, LFPfilesMicr{f}), 'times', 'chan_locations');
    LFPtimeMicr = micrObj.times;
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    
    micrObj = load(fullfile(dataDirectory, sprintf('%s_1kHz_unitsRemoved.mat', subject)));
    
    %% Parse channel locations and hemispheres
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@(X) isempty(X), locations)) = [];
    
    % Convert to standardized region names
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
    locationsHem = locations;
    
    % Add hemisphere labels
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locationsHem(hemi) = cellfun(@(X) ['L' X], locationsHem(hemi), 'UniformOutput', false);
    hemi = contains(hem, 'right');
    locationsHem(hemi) = cellfun(@(X) ['R' X], locationsHem(hemi), 'UniformOutput', false);
    
    %% Load task data
    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectoryOrig, '../task', taskfiles{f}));
    
    %% Load ripple statistics
    tag = [recordingState, '_', location, '_', modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    
    % Format channel information
    LFPchanNum = cellfun(@(X) str2double(X), rippleStats.chanLabels);
    rippleStats.chanNum = rippleStats.chanLabels;
    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations)
        error('Error formatting channel labels.');
    end
    
    rippleStats.recordingType = repmat({'micro'}, [1 length(microObj.rippleStats.locs)]);
    microObj.rippleStats.window(cellfun(@(X) isempty(X), microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', ...
                                 microObj.rippleStats.window, 'UniformOutput', false);
    rippleStats.microTimes = LFPtimeMicr;
    
    %% Create ripple mask
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask, 1)
        if rippleStats.density{chRipp} <= 1
            continue
        end
        
        if strcmp(rippleStats.recordingType{chRipp}, 'macro')
            continue
        elseif strcmp(rippleStats.recordingType{chRipp}, 'micro')
            times = rippleStats.microTimes;
        end
        
        % Skip special electrodes
        if contains(rippleStats.chanLabels{chRipp}, 'SPE')
            rippMask(chRipp, :) = nan(1, length(rippMask));
            continue
        end
        
        % Collect ripple metrics
        densityAll = [densityAll rippleStats.density{chRipp}];
        freqAll = [freqAll mean(rippleStats.oscFreq{chRipp})];
        durAll = [durAll mean(rippleStats.duration{chRipp})];
        ampAll = [ampAll mean(rippleStats.rippleAmp{chRipp})];
        locAll = [locAll locations(chRipp)];
        subjAll = [subjAll subj];
        
        %% Load ERSP data
        try
            ERSPdir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/Figures/RutishauserLab/Sternberg/%s/NC/TF_wake_1kHz_template_z25/Data', subject);
            load(fullfile(ERSPdir, sprintf('%s_%s_TimeFreq.mat', subject, rippleStats.chanNum{chRipp})));
            erspTimes = times;
            ERSPall(:, :, cnt) = ersp;
            cnt = cnt + 1;
        catch
            ERSPall(:, :, cnt) = nan(size(ersp));
            cnt = cnt + 1;
        end
        
        %% Populate ripple mask
        iS = round(rippleStats.window{chRipp}(:, 1) * 1e3);
        iE = round(rippleStats.window{chRipp}(:, 2) * 1e3);
        
        for ii = 1:length(iE)
            if any([iS(ii) iE(ii)] <= 0)
                continue
            end
            rippMask(chRipp, iS(ii):iE(ii)) = 1;
        end
    end
    
    % Pad ending with NaN
    rippMask(:, end:end+3e3) = nan;
    rippAll = sum(rippMask);
    
    % Load LFP data and pad
    data = micrObj.data;
    data(:, end:end+3e3) = nan;
    
    %% Load and process unit data
    units = LoadSpikeTimes(subject, 'RutishauserLab', 'Sternberg');
    if isempty(units)
        continue
    end
    
    % Filter for valid unit types
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:, 3));
    uChan = cell2mat(units(:, 1));
    uChanAll = [uChanAll uChan'];
    
    % Convert unit locations to standardized names
    uLocations = units(U, 5);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
    
    % Collect unit metrics
    troughToPeak = [troughToPeak, cell2mat(units(:, 6))'];
    burstIndex = [burstIndex, cell2mat(units(:, 7))'];
    uSubj = [uSubj repmat(subj_list_full(subj), [1, size(units, 1)])];
    uID = [uID 1:size(units, 1)];
    uType = [uType units(:, 3)'];
    
    %% Create spike mask
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001;  % seconds
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    
    for iU = 1:nNeuron 
        X = units{unitsAll(iU), 2};
        [N, EDGES] = histcounts(X, binEdges);
        FR = [FR sum(N) / rippleStats.recordingLength * 1000];
        spikeMask(iU, :) = N;
        spikeMask_sm(iU, :) = smoothdata(N, 'gaussian', 100);
    end
    
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian', 100);
    else
        unitAll = smoothdata(spikeMask, 'gaussian', 100);
    end
    unitAll = zscore(unitAll);
    
    %% Collect ripple LFP by region
    for iR = 1:length(regions)
        chRegion = find(contains(rippleStats.chanLabels, regions(iR)));
        uRegion = find(contains(uLocations, regions(iR)));
        cchRegion(iR) = cchRegion(iR) + length(chRegion);
        cUnitRegion(iR) = cUnitRegion(iR) + length(uRegion);
        
        for ii = 1:length(chRegion)
            ch = chRegion(ii);
            rippleCenter = rippleStats.locs{ch};
            for iRpl = 1:10:length(rippleCenter)
                rip = data(ch, rippleCenter(iRpl)-500:rippleCenter(iRpl)+500);
                ripple_LFP_region_all{iR}(cRegion(iR), :) = rip;
                cRegion(iR) = cRegion(iR) + 1;
            end
        end
    end
    
    %% Create task mask for trial epochs
    taskMask = false(3, length(rippMask));
    
    for iT = 1:length(trials.start_time)
        trialLoad = trials.loads(iT);
        image1Time = round(trials.timestamps_Encoding1(iT) * 1e3);
        image2Time = round(trials.timestamps_Encoding2(iT) * 1e3);
        image3Time = round(trials.timestamps_Encoding3(iT) * 1e3);
        maintenanceTime = round(trials.timestamps_Maintenance(iT) * 1e3);
        probeTime = round(trials.timestamps_Probe(iT) * 1e3);
        respTime = round(trials.timestamps_Response(iT) * 1e3);
        
        % Define trial epochs
        taskMask(1, image1Time-1000:image1Time-1) = true;  % Baseline
        taskMask(2, image1Time:image1Time+1e3) = true;  % Image 1
        taskMask(3, probeTime+1:probeTime+1e3) = true;  % Probe 
        
        if trialLoad == 3
            taskMask(2, image1Time:image1Time+1e3) = true;  % Image 1
            taskMask(2, maintenanceTime:probeTime) = true;  % Maintenance 
        end
    end
    
    %% Calculate co-ripple statistics across regions
    rippMask(isnan(rippMask)) = 0;
    
    for iRa = 1:length(regionsAll)
        for iRb = iRa:length(regionsAll)
            regionA = find(contains(locationsHem, regionsAll(iRa)));
            regionB = find(contains(locationsHem, regionsAll(iRb)));
            
            if isempty(regionA) || isempty(regionB)
                continue
            end
            
            for iiA = 1:length(regionA)
                for iiB = 1:length(regionB)
                    chA = regionA(iiA);
                    chB = regionB(iiB);
                    
                    if chA == chB
                        continue
                    end
                    
                    if isempty(mask2bounds(rippMask(chB, :))) || isempty(mask2bounds(rippMask(chA, :)))
                        continue
                    end
                    
                    % Calculate co-ripple events
                    coR = mask2bounds(rippMask(chA, :) & rippMask(chB, :));
                    coR = coR(:, 2) - coR(:, 1);
                    coRnum(ccgCount(iRa, iRb), iRa, iRb) = sum(coR > 25);
                    coRdur(ccgCount(iRa, iRb), iRa, iRb) = length(rippleStats.locs{chA});
                    
                    ccgCount(iRa, iRb) = ccgCount(iRa, iRb) + 1;
                end
            end
        end
    end
    
    toc
end

% Clean up unused trial data
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) = [];
probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('Done processing.\n')

%%   VISUALIZATION AND ANALYSIS

%% Figure: Ripple LFP waveforms by region (zoomed view)
regionColors = brewermap(12, 'Dark2');
figure('Position', [6 1353 950 260/2]);

for iR = 1:length(regions)
    subplot(1, 5, iR)
    
    % Create shaded background for ripple period
    pa = patch([50 -50 -50 50], [-11 -11 16 16], 'k');
    hold on;
    pa.FaceColor = [245 204 85]/255;
    pa.FaceAlpha = 0.2;
    pa.LineStyle = '--';
    pa.LineWidth = 1;
    
    vl = vline(0);
    vl.LineWidth = 1;
    vl.LineStyle = '-';
    
    % Plot mean +/- SEM
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = std(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan(ripple_LFP_region_all{iR}(:, 1))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap');
    
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR, :);
    bl1.LineWidth = 1;
    bf.FaceColor = regionColors(iR, :);
    bf.FaceAlpha = 0.3;
    
    xlim([-150 150])
    ylim([-12 20])
    ax = gca;
    ax.XTick = [-150, 150];
    ax.YTick = [-10, 0 20];
    ax.FontSize = 12;
    ax.LineWidth = 1;
    
    fprintf('%s: %i ripples\n', regions{iR}, sum(~isnan(ripple_LFP_region_all{iR}(:, 1))))
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_zoom_%s.pdf', tag)))

%% Figure: Ripple LFP waveforms (ultra-zoomed, no axes)
figure('Position', [6 1353 950 220/2]);

for iR = 1:length(regions)
    subplot(1, 5, iR)
    
    % Create shaded background
    pa = patch([50 -50 -50 50], [-11 -11 16 16], 'k');
    hold on;
    pa.FaceColor = [245 204 85]/255;
    pa.FaceAlpha = 0.2;
    pa.LineStyle = '--';
    pa.LineWidth = 1;
    
    vl = vline(0);
    vl.LineWidth = 1;
    vl.LineStyle = '-';
    
    % Plot mean +/- SEM
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = std(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan(ripple_LFP_region_all{iR}(:, 1))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap');
    
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR, :);
    bl1.LineWidth = 1;
    bf.FaceColor = regionColors(iR, :);
    bf.FaceAlpha = 0.3;
    
    xlim([-52 52])
    ylim([-12 20])
    ax = gca;
    ax.XTick = [-150, 150];
    axis off
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_2zoom_%s.pdf', tag)))

%% Figure: Ripple LFP waveforms (full window)
figure('Position', [-64 1149 2353 424]);

for iR = 1:length(regions)
    subplot(1, 5, iR)
    
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = std(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan(ripple_LFP_region_all{iR}(:, 1))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap');
    
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR, :);
    bf.FaceColor = regionColors(iR, :);
    bf.FaceAlpha = 0.3;
    
    xlim([-500 500])
    ylabel('amplitude [uV]')
    xlabel('time from ripple center')
    title(regions{iR})
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_%s.pdf', tag)))

%% Figure: Channel and unit distribution pie charts
figure('Position', [1211 1267 781 294]); 

% Channel locations
subplot(1, 3, 1)
h = pie(cchRegion, regions);
patchHandles = findobj(h, 'Type', 'patch');
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = regionColors(ii, :);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1, pos(2), num2str(cchRegion(ii)));
    hold on;
end
title('channel locations')

% Unit locations
subplot(1, 3, 2)
h = pie(cUnitRegion, regions);
patchHandles = findobj(h, 'Type', 'patch');
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = regionColors(ii, :);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1, pos(2), num2str(cUnitRegion(ii)));
    hold on;
end
title('unit locations')

% Unit types
subplot(1, 3, 3)
types = {'pyr', 'int'}; 
c = [86, 180, 233; 204, 121, 167] / 255;
h = pie([sum(strcmp(uType, 'pyr')) sum(strcmp(uType, 'int'))], types);
patchHandles = findobj(h, 'Type', 'patch');
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = c(ii, :);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1, pos(2), num2str(sum(strcmp(uType, types{ii}))));
    hold on;
end
title('unit types')

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'LFP_Unit_Count.pdf'))

%% Figure: Unit classification (trough-to-peak vs burst index)
figure; 
validIdx = troughToPeak >= 45 & isfinite(burstIndex);
pl = plot(troughToPeak(validIdx), log10(burstIndex(validIdx) - min(burstIndex)), '.');
hold on;
pl.Color = c(1, :);

validIdx = troughToPeak < 45 & isfinite(burstIndex);
pl = plot(troughToPeak(validIdx), log10(burstIndex(validIdx) - min(burstIndex)), '.');
pl.Color = c(2, :);

ylabel('burstIndex')
xlabel('troughToPeak [samples]')

%% Figure: Unit classification (trough-to-peak vs firing rate)
figure; 
validIdx = troughToPeak >= 45 & isfinite(burstIndex);
pl = plot(troughToPeak(validIdx), FR(validIdx), '.');
hold on;
pl.Color = c(1, :);

validIdx = troughToPeak < 45 & isfinite(burstIndex);
pl = plot(troughToPeak(validIdx), FR(validIdx), '.');
pl.Color = c(2, :);

ylabel('FR')
xlabel('troughToPeak [samples]')

%% Figure: Burst index distribution by cell type
figure; 
validIdx = troughToPeak > 45 & isfinite(burstIndex);
histogram(log10(burstIndex(validIdx) - min(burstIndex)), -2:0.1:3, 'Normalization', 'probability');
hold on;

validIdx = troughToPeak < 45 & isfinite(burstIndex);
histogram(log10(burstIndex(validIdx) - min(burstIndex)), -2:0.1:3, 'Normalization', 'probability')

%% Figure: Event-Related Spectral Perturbation (ERSP)
fs = rippleStats.fs;
figure('Position', [503 555 round(367*430/528) 318]);

erspMu = mean(ERSPall(:, erspTimes >= -1.0e3 & erspTimes < 1.0e3, :), 3, 'omitnan');
imagesc([-1.0e3 1.0e3], [1 fs/2], erspMu);
hold on

ax = gca;
ax.YDir = 'normal';
ax = imgca;
set(ax, 'YScale', 'log')
ax.LineWidth = 1;

yticks([1 2 4 10 16 50 70 100 200 fs/2])
caxis([-1.7 1.7])
ylim([1 fs/2])
ylim([50 200])

v = vline(0, 'k:');
v.LineWidth = 1.5;

xlabel('time from ripple center [ms]')
ylabel('frequency (Hz)')
c = colorbar;
c.Label.String = 'ERSP (dB)';
c.Label.Interpreter = 'tex';
c.Ticks = [-1.5 1.5];
c.LineWidth = 1;

set(gca, 'FontSize', 15)
set(gcf, 'Color', [1 1 1])
colormap(fliplr(slanCM('fusion', 501)))

savepdf(gcf, fullfile(exportDirFigs, 'LFP_ripple_ERSP_z.pdf'))

%% Figure: Ripple metrics violin plots by region
close all
figure('Position', [21 997 809 251]);

% Ripple density
subplot(1, 4, 1)
vp = violinplot(densityAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR, :);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple density [/min]';
box off

% Ripple duration
subplot(1, 4, 2)
vp = violinplot(durAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR, :);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple duration [ms]';
ylim([50 120]);
box off

% Ripple amplitude
subplot(1, 4, 3)
vp = violinplot(ampAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR, :);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple amplitude [uV]';
ylim([0 25]);
box off

% Ripple frequency
subplot(1, 4, 4)
vp = violinplot(freqAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR, :);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple frequency [Hz]';
ylim([80 100]);
box off

% Print summary statistics
fprintf('\nRipple frequency by region:\n')
cellfun(@(X) fprintf('%s: %.2f +/- %.2f Hz\n', X, mean(freqAll(strcmp(locAll, X))), ...
    std(freqAll(strcmp(locAll, X)))), locs, 'UniformOutput', false);

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'rippleMetrics.pdf'))

%% Figure: Ripple metrics by subject
figure('Position', [1 68 2056 1099]);

% Density by subject
subplot(4, 1, 1)
vp = violinplot(densityAll, subjAll);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple density [/min]';
box off

% Duration by subject
subplot(4, 1, 2)
vp = violinplot(durAll, subjAll);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple duration [ms]';
ylim([50 120]);
box off

% Amplitude by subject
subplot(4, 1, 3)
vp = violinplot(ampAll, subjAll);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple amplitude [uV]';
ylim([0 25]);
box off

% Frequency by subject
subplot(4, 1, 4)
vp = violinplot(freqAll, subjAll);
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0;
    vp(ii).ShowWhiskers = 0;
end
ax.YLabel.String = 'ripple frequency [Hz]';
ylim([80 100]);
box off

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'rippleMetrics_bySubj.pdf'))


%% Analysis: Co-ripple rates and anatomical connectivity
clc

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
pdfExport = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));

% Anatomical parcellation mapping
ctxParc = {'ACC', 'SMA', 'OFC', 'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};

% Define regions for analysis
iRa = 4; 
iRapl = find(contains(regionsAll, regions([1:5])));
iRbpl = find(contains(regionsAll, regions([1:5])));

figure; 
c = 1;
coRrate = [];
x = [];  % Tract length
yRip = [];  % Co-ripple rate
yShuff = [];  % Shuffled co-ripple rate
hemA = [];
hemB = [];
withinBundle = [];

% Loop through all region pairs
for iRa = 1:length(regionsAll)
    for iRb = iRa:length(regionsAll)
        fprintf('%s %s\n', regionsAll{iRa}, regionsAll{iRb})
        
        % Calculate co-ripple rates
        co = coRnum(1:ccgCount(iRa, iRb)-1, iRa, iRb);
        N = coRdur(1:ccgCount(iRa, iRb)-1, iRa, iRb);
        temp = co ./ N;
        
        mu = mean(temp);
        sem = std(temp) / sqrt(length(temp));
        
        % Extract hemisphere and parcel information
        parcelA = regionsAll{iRa}(2:end);
        parcelB = regionsAll{iRb}(2:end);
        
        hemA = [hemA repmat({regionsAll{iRa}(1)}, [1 length(temp)])];
        hemB = [hemB repmat({regionsAll{iRb}(1)}, [1 length(temp)])];
        samebund = iRa == iRb;
        withinBundle = [withinBundle repmat(samebund, [1 length(temp)])];
        
        % Map to Broadmann areas for anatomical connectivity
        broadmanA = broadman{strcmp(ctxParc, parcelA)};
        HCPa = sprintf('%s_%s', hemA{end}, broadmanA);
        broadmanB = broadman{strcmp(ctxParc, parcelB)};
        HCPb = sprintf('%s_%s', hemB{end}, broadmanB);
        
        % Get tract length from connectivity matrix
        d = tractLenStruct.tractLengths(strcmp(tractLenStruct.parcelIDs, HCPa), ...
                                        strcmp(tractLenStruct.parcelIDs, HCPb));
        
        x = [x repmat(d, [1 length(temp)])];
        yRip = [yRip temp'];
        
        c = c + 1;
    end
end

% Bin data by tract length quantiles
clr = brewermap(10, 'Paired');
yRipMu = [];
yShufMu = [];
yRipSEM = [];
yShufSEM = [];
xMu = [];
xBin = quantile(x, 5);
xBin = [0 xBin max(x)];

for iB = 2:length(xBin)
    ii = x >= xBin(iB-1) & x <= xBin(iB);
    yRipMu = [yRipMu mean(yRip(ii))];
    yShufMu = [yShufMu mean(yShuff(ii))];
    yRipSEM = [yRipSEM std(yRip(ii)) / sqrt(sum(ii))];
    yShufSEM = [yShufSEM std(yShuff(ii)) / sqrt(sum(ii))];
    xMu = [xMu mean([xBin(iB-1) xBin(iB)])];
end

% Plot co-ripple rate vs tract length
[bl2, bf] = boundedline(xMu, yRipMu, yRipSEM, 'rs-');
hold on;
bl2.Color = clr(8, :);
bl2.MarkerFaceColor = clr(8, :);
bf.FaceColor = clr(8, :);
bf.FaceAlpha = 0.3;

[bl2, bf] = boundedline(xMu, yShufMu, yShufSEM, 'bs-');
bl2.Color = clr(1, :);
bl2.MarkerFaceColor = clr(1, :);
bf.FaceColor = clr(1, :);
bf.FaceAlpha = 0.3;

fig = gcf;
fig.Color = 'w';

%% Figure: Co-ripple rates - within vs between bundles/hemispheres
figure('Position', [100 368 round(357*(2/3)) round(155)]);

% Within bundle co-ripple rates
subplot(1, 3, 1)
groups = nan(1, length(yRip));
groups(withinBundle == 1) = 1;
groups(withinBundle == 0) = [];

bx = boxplot(yRip(withinBundle == 1), groups, 'notch', 'off', 'Symbol', '', 'width', 0.3);
ylim([0 0.45])
ax = gca;
ax.FontSize = 7;
ax.LineWidth = 1;
set(bx, 'LineWidth', 1);

% Customize box appearance
boxes = findobj(gca, 'Tag', 'Box');
fillColor = [0.7, 0.9, 1];
edgeColor = [0, 0, 0.8];

for i = 1:length(boxes)
    xdata = get(boxes(i), 'XData');
    ydata = get(boxes(i), 'YData');
    p = patch(xdata, ydata, fillColor, 'FaceAlpha', 0.5, ...
              'EdgeColor', edgeColor, 'LineWidth', 1.0);
    uistack(boxes(i), 'top');
end
box off
xlim([0.5 1.5])

% Between bundle co-ripple rates (same vs different hemisphere)
subplot(1, 3, [2 3])
groups = nan(1, length(yRip));
groups(strcmp(hemA, hemB) & withinBundle == 0) = 1;  % Same hemisphere
groups(~strcmp(hemA, hemB) & withinBundle == 0) = 2;  % Different hemisphere
groups(withinBundle == 1) = [];

bx = boxplot(yRip(withinBundle == 0), groups, 'notch', 'off', 'Symbol', '');
ax = gca;
ax.FontSize = 7;
ax.LineWidth = 1;
set(bx, 'LineWidth', 1);

% Customize box appearance
boxes = findobj(gca, 'Tag', 'Box');
for i = 1:length(boxes)
    xdata = get(boxes(i), 'XData');
    ydata = get(boxes(i), 'YData');
    p = patch(xdata, ydata, fillColor, 'FaceAlpha', 0.5, ...
              'EdgeColor', edgeColor, 'LineWidth', 1.0);
    uistack(boxes(i), 'top');
end

ylim([0 0.1]);
box off

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'coRipRates.pdf'))