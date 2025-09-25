

close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath(    '/space/seh10/6/halgdev/projects/iverzh/ripples/code/ripple-detection/code/CCG_tools'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%
dataDirectoryOrig = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
flst = dir(fullfile(dataDirectoryOrig, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;


matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

unitfiles = dir(fullfile(dataDirectoryOrig, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectoryOrig, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectoryOrig, '*micro*'));

LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectoryOrig, '../task','*task*'));
taskfiles = {taskfiles.name}';
bpFiles = dir(fullfile(dataDirectoryOrig, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectoryOrig, '../out', '*micro*'));
micrFiles = {micrFiles.name}';

recordingState = 'wake';
location = 'NC';
tag = 'vHG';
sfreq = 1e3;
copmuteTF = false; 
regionsAll = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end,3:end-1));   
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2,3:end-1));
bundleNotes = table2cell(bundleNotes(3:end,3:end-1));   

%CCG
win=1000;
binSize = 1; %ms
hist_bins = 2*win/binSize;
nIter = 200;
FDRwin        = 500; % +- window (in ms) to compute FDR significance


%%


%prealocate variables

countChPerRegion = zeros(1, length(regions));
countUnitPerRegion = zeros(1, length(regions));
ripple_LFP_region_all = cell(length(regions), 1);
for iRa = 1:length(regions)
    ripple_LFP_region_all{iRa} = nan(1e6, length(-500:500));       

end
cRegion = ones(1,length(regions));
cchRegion = ones(1,length(regions));
cUnitRegion = ones(1,length(regions));
subjID = nan(6.5e3, 1);
correct_trials_all = nan(6.5e3, 1);
load_trials_all = nan(6.5e3, 1);
probe_in_out = nan(6.5e3, 1);
respLatency = nan(6.5e3, 1);
densityAll = [];
freqAll = [];
durAll = [];
ampAll = [];
locAll = [];
subjAll = [];
burstIndex = [];
troughToPeak = [];
FR = [];
uType = [];
uSubj = [];
uID = [];
uChanAll = [];
ERSPall = [];
% tic
cTrial = 1;
cnt = 1;
coRnum = nan(1e5, length(regionsAll),length(regionsAll));
coRnumShuff = nan(1e5, length(regionsAll),length(regionsAll));
coRdur = nan(1e5, length(regionsAll),length(regionsAll));
coR_PRTH = zeros(hist_bins-1, 2, 1e5, length(regionsAll),length(regionsAll));
coR_PRTH_null = zeros(hist_bins-1, 1e5, length(regionsAll),length(regionsAll));
ccgNrip = zeros(2, 1e5, length(regionsAll),length(regionsAll));
ccgCount = ones(length(regionsAll));
for subj =1:length(subj_list_full)
    tic
    subject = subj_list_full{subj};

    fprintf('processing %s ... \n', subject)



    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectoryOrig, LFPfilesMicr{f}), 'times', 'chan_locations');
    LFPtimeMicr = micrObj.times;
%     chan_labels = strsplit(micrObj.chan_locations, micrObj.chan_locations(38))
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    
    micrObj = load(fullfile(dataDirectory,sprintf('%s_1kHz_unitsRemoved.mat', subject)));

    
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@(X) isempty(X), locations)) = [];
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
    locationsHem = locations;
    
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locationsHem(hemi) = cellfun(@(X) ['L' X], locationsHem(hemi),  'UniformOutput', false);
    hemi = contains(hem, 'right');
    locationsHem(hemi) = cellfun(@(X) ['R' X], locationsHem(hemi),  'UniformOutput', false);
    
    

    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectoryOrig, '../task', taskfiles{f}));

    
%     modifier = 'IISremv_singleWire_v3_z25';
    modifier = '1kHz_template_z25';

%     modifier = 'IISremv_singleWire_v3';
%     modifier = 'IISremv_v3';
%     modifier = 'LG_1kHz_z25';

    tag = [recordingState,'_',location,'_',modifier];
%     filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
    filename = sprintf('%s_ripple_stats_%s.mat', subject,  tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    LFPchanNum = cellfun(@(X) str2double(X), rippleStats.chanLabels);
    rippleStats.chanNum = rippleStats.chanLabels;
    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations); error('error formatting channel labels.'); end
        
    rippleStats.recordingType = [repmat({'micro'}, [1 length(microObj.rippleStats.locs)])];
%     rippleStats.locs = [cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false)];
    microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = [cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', microObj.rippleStats.window, 'UniformOutput', false)];

    rippleStats.microTimes = LFPtimeMicr;
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    for chRipp = 1:size(rippMask,1) %sum(strcmp(rippleStats.recordingType, 'macro'))+1:size(rippMask,1) 
        if rippleStats.density{chRipp} > 1
%             if strcmp(rippleStats.recordingType{chRipp}, 'macro'); times = rippleStats.macroTimes;
            if strcmp(rippleStats.recordingType{chRipp}, 'macro'); continue;
            elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); times = rippleStats.microTimes; 
%             elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); continue; 
            end

            if contains(rippleStats.chanLabels{chRipp}, 'SPE')
                rippMask(chRipp,:) = nan(1,length(rippMask));
                continue
            end
            
            
%             chLoc = regexprep(rippleStats.chanLabels{chRipp}, '\d', '');
%             chNum = regexp(rippleStats.chanLabels{chRipp}, '\d+', 'match');
%             chNum = str2double(chNum);
%             epCh = find(contains(bundleLab, chLoc));
%             
%             bch = badChan(subj, epCh);
%             bch = str2double(split(bch,','));
%             
%             if contains(bundleNotes(subj, epCh), 'Bad') || any(ismember(bch, chNum))
%                 rippMask(chRipp,:) = nan(1,length(rippMask));
%                 continue
%             end
            

            densityAll = [densityAll rippleStats.density{chRipp}];
            freqAll = [freqAll mean(rippleStats.oscFreq{chRipp})];
            durAll = [durAll mean(rippleStats.duration{chRipp})];
            ampAll = [ampAll mean(rippleStats.rippleAmp{chRipp})];
            locAll = [locAll locations(chRipp)];
            subjAll = [subjAll subj];

            
            try
                ERSPdir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/Figures/RutishauserLab/Sternberg/%s/NC/TF_wake_1kHz_template_z25/Data', subject);
                load(fullfile(ERSPdir, sprintf('%s_%s_TimeFreq.mat',subject, rippleStats.chanNum{chRipp})));
                erspTimes = times;
                ERSPall(:,:,cnt) = ersp;
                cnt = cnt + 1;
            catch
                ERSPall(:,:,cnt) = nan(size(ersp));
                cnt = cnt + 1;
%                 warning('could not load ERSP for %s', subject)
            end
            
            
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
            
            
        
            for ii = 1:length(iE)
%                 iS(ii) = find(times == iS(ii));
%                 iE(ii) = find(times == iE(ii));{
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
        
    end

    
    
    rippMask(:, end:end+3e3) = nan; %pad the ending
    rippAll = sum(rippMask);

    
    
    data = micrObj.data;
%     figure;
%     if all(contains({'AMY', 'HIP', 'SMA'}, locations))
%         rippAMY = sum(rippMask(contains(locations, {'AMY'}), :));
%         rippHIP = sum(rippMask(contains(locations, {'HIP'}), :));
%         rippSMA = sum(rippMask(contains(locations, {'SMA'}), :));
%         
%         b = mask2bounds((rippAMY > 0) & (rippHIP > 0) & (rippSMA > 0));
%         dur = b(:,2) - b(:,1);
%         
%         
%         ch = find(rippMask(:,find(dur == max(dur), 1,'first'), 1));
%         
%         for ii = 1:length(ch)
%             dat = data(ch(ii), b(find(dur == max(dur), 1,'first'), 1)-100: b(find(dur == max(dur), 1,'first'), 2)+100);
%             dat = dat / max(abs(dat));
%             plot(dat); hold on;
%         end
%         
%         savepdf(gcf, fullfile(exportDirFigs, sprintf('%s_coRipSingleTrial.pdf', subject)))
% 
%         
%         
%     end
%     
%     close;
%     [b,a] = butter(3,[150]/(sfreq/2), 'high');
%     [bB,aB] = butter(3,[70 190]/(sfreq/2));
%     [bR,aR] = butter(3,[70 100]/(sfreq/2));
%     [bL,aL] = butter(3,[30 50]/(sfreq/2));
%     vHG_band = NaN(size(data));
%     vHG_band_pow = NaN(size(data));
% 
%     bHG_band_pow = NaN(size(data));
%     rHG_band_pow = NaN(size(data));
%     rHG_band_phi = NaN(size(data));
%     LG_band_pow = NaN(size(data));
% 
%     for ch = 1:size(data,1)
%         dat = data(ch,:);
% %         dat(nan_edge_mask) = 0;
% % 
%         for nf = 60:60:sfreq/2 %240
%             Wo = nf/(sfreq/2);
%             BW = Wo/35;
%             [bN,aN] = iirnotch(Wo, BW);
%             dat = filtfilt(bN,aN,dat);
%         end
%         data(ch,:) = dat;
%         
%         bp = filtfilt(b,a,dat);
%         pow = abs(hilbert(bp));
%         ang = angle(hilbert(bp));
% 
% %         bp(nan_edge_mask) = nan;
% %         ang(nan_edge_mask) = nan;
%         rHG_band_phi(ch,:) = ang;
%         
%         
% 
%         vHG_band(ch,:) = zscore(abs(bp));
%         vHG_band_pow(ch,:) = zscore(pow);
% 
%         bp = filtfilt(bB,aB,dat);
%         pow = abs(hilbert(bp));
%         bHG_band_pow(ch,:) = zscore(pow);
% 
%         bp = filtfilt(bR,aR,dat);
%         pow = abs(hilbert(bp));
%         rHG_band_pow(ch,:) = zscore(pow);
% 
%         bp = filtfilt(bL,aL,dat);
%         pow = abs(hilbert(bp));
%         LG_band_pow(ch,:) = zscore(pow);
% 
%     end
%     vHG_band(vHG_band < 2) = 0;
%     vHG_band(vHG_band >= 2) = 1;
%     
%     
%     LG_band_pow(:, end:end+3e3) = nan; %pad the ending
%     rHG_band_pow(:, end:end+3e3) = nan; %pad the ending
%     vHG_band_pow(:, end:end+3e3) = nan; %pad the ending
    data(:, end:end+3e3) = nan; %pad the ending
%     
%     data = rHG_band_pow;


    
    
    units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    uLocations = units(U,5);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
        
    troughToPeak = [troughToPeak, cell2mat(units(:,6))'];
    burstIndex = [burstIndex, cell2mat(units(:,7))'];
    uSubj = [uSubj repmat(subj_list_full(subj), [1, size(units,1)])];
    uID = [uID 1:size(units,1)];
    uType = [uType units(:,3)'];
    
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001; %seconds
    
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    for iU = 1:nNeuron 
        X = units{unitsAll(iU),2};
    %     ISIx = [0 diff(X')];
    %     xx = randperm(length(ISIx), length(ISIx));
    %     Xshuff =[X(1)];
    %     for ii = 2:length(ISIx); Xshuff(ii) = Xshuff(end) + ISIx(xx(ii)); end
        [N,EDGES] = histcounts(X,binEdges);
        FR = [FR sum(N) / rippleStats.recordingLength * 1000];
    %     if FR > 0
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian',100);
    %     end
    end
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian',100);
    else
        unitAll = smoothdata(spikeMask, 'gaussian',100);
    end
    unitAll = zscore(unitAll); % - mean(unitAll);

    
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
    
    taskMask = false(3, length(rippMask));

    for iT = 1:length(trials.start_time)

        trialLoad = trials.loads(iT);
        image1Time =  round(trials.timestamps_Encoding1(iT)*1e3);
        image2Time =  round(trials.timestamps_Encoding2(iT)*1e3);
        image3Time =  round(trials.timestamps_Encoding3(iT)*1e3);
        maintenanceTime =  round(trials.timestamps_Maintenance(iT)*1e3);
        probeTime =  round(trials.timestamps_Probe(iT)*1e3);
        respTime =  round(trials.timestamps_Response(iT)*1e3);

        taskMask(1, image1Time-1000:image1Time-1) = true; %baseline
        taskMask(2, image1Time:image1Time+1e3)    = true; %img1
        taskMask(3, probeTime+1:probeTime+1e3)    = true; %probe 

        if trialLoad == 3
            taskMask(2, image1Time:image1Time+1e3)    = true; %img1
            taskMask(2, maintenanceTime:probeTime)    = true; %maintenance 
 
        end



    end
    rippMask(isnan(rippMask)) = 0;
%     rippMask = logical(rippMask);
    for iRa = 1:length(regionsAll)
        for iRb = iRa:length(regionsAll)
            regionA = find(contains(locationsHem, regionsAll(iRa)));
            regionB = find(contains(locationsHem, regionsAll(iRb)));
            
            if isempty(regionA) || isempty(regionB); continue; end
            for iiA = 1:length(regionA)
                for iiB = 1:length(regionB)
% 
                    chA = regionA(iiA); chB = regionB(iiB);
                    if chA == chB; continue; end
                    
                    if isempty(mask2bounds(rippMask(chB,:))) || isempty(mask2bounds(rippMask(chA,:)))
                        continue
                    end

                    coR = mask2bounds(rippMask(chA,:) & rippMask(chB,:));
                    coR = coR(:,2) - coR(:,1);
                    coRnum(ccgCount(iRa, iRb),iRa, iRb) = sum(coR > 25);
                    coRdur(ccgCount(iRa, iRb),iRa, iRb) = length(rippleStats.locs{chA});
%                     
%                     [event_PRTH, null_PRTH] = computePRTH(rippleStats, recordingState, chA, win, binSize, ...
%                                               1, '', 'NC', 1, subject, chB,taskMask(2,:),0,0);
%                                           
%                                           
%                     coR_PRTH_null(:,  ccgCount(iRa, iRb), iRa, iRb) = null_PRTH';
%                     coR_PRTH(:,1,ccgCount(iRa, iRb), iRa, iRb) = event_PRTH';
%                     locs = rippleStats.locs{chA};
%                     ccgNrip(1, ccgCount(iRa, iRb), iRa, iRb) = length(locs(taskMask(2,locs)));
%                     
%                     [event_PRTH, null_PRTH] = computePRTH(rippleStats, recordingState, chA, win, binSize, ...
%                                               1, '', 'NC', 1, subject, chB,taskMask(3,:),0,0);
%                                           
%                                           
%                     coR_PRTH_null(:,  ccgCount(iRa, iRb), iRa, iRb) = null_PRTH';
%                     coR_PRTH(:,2,ccgCount(iRa, iRb), iRa, iRb) = event_PRTH';
%                     ccgNrip(2, ccgCount(iRa, iRb), iRa, iRb) = length(locs(taskMask(3,locs)));
%                     
                    ccgCount(iRa, iRb) = ccgCount(iRa, iRb) + 1;
%                      
%                                           
%                     
% %                     for iter = 1:1
% %                         shuffA = shuffleMask(rippMask(chA,:));
% %                         shuffB = shuffleMask(rippMask(chB,:));
% %                         coR = mask2bounds(shuffA & shuffB);
% %                         coR = coR(:,2) - coR(:,1);
% % %                         coRnumShuff(iRa, iRb) = coRnumShuff(iRa, iRb) + sum(coR > 25);
% %                         coRnumShuff(ccgCount(iRa, iRb),iRa, iRb) = sum(coR > 25);
% 
% %                     end
                end
            end
%             
%             
%                 
%                 
%             
        end
    end
 

%     save(fullfile(exportDir, 'ripDetections', sprintf('%s_rippleStats.mat', subject)), 'rippleStats' ,'-v7.3')
    toc
end
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) =[];

probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('done processing .. \n')
toc

% [row, col] = find(isnan(recog_trials_all));
% recog_trials_all_resp = nan(size(recog_trials_all,1), 1.5e3+1);
% for iT = 1:size(recog_trials_all_resp,1)
%     I = min(col(row == iT));
%     if I  == 1; continue; end
%     trial = recog_trials_all(iT, I-1500:I);
%     recog_trials_all_resp(iT, 1:length(trial)) = trial;
% 
% 
% end

%%
close all
regionColors =  brewermap(12, 'Dark2');
figure('Position',[6 1353 950 260/2]);


for iR = 1:length(regions)
    subplot(1,5,iR)
    
    ylim([-12 20])
    pa = patch([50 -50 -50 50], [-11 -11 16 16], 'k'); hold on;
    pa.FaceColor = [245 204 85]/255;
    pa.FaceAlpha = 0.2;
    pa.LineStyle = '--';
    pa.LineWidth = 1;
    vl = vline(0); vl.LineWidth = 1; vl.LineStyle = '-';
    
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = mean(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan((ripple_LFP_region_all{iR}(:,1)))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bl1.LineWidth =1;

    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
%     text(-75, 19, sprintf('%i ripples', sum(~isnan((ripple_LFP_region_all{iR}(:,1)))))); hold on;
    xlim([-150 150])
    ylim([-12 20])
    ax = gca;
    ax.XTick = [-150, 150];
    ax.YTick = [-10, 0 20];
    ax.FontSize = 12;
    ax.LineWidth = 1;
%     ylabel('amplitude [uV]')
%     xlabel('time from ripple center')
    sum(~isnan((ripple_LFP_region_all{iR}(:,1))))
%     title(regions{iR})
 
end


fig = gcf;
fig.Color = 'w';


savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_zoom_%s.pdf', tag)))
figure('Position',[6 1353 950 220/2]);

for iR = 1:length(regions)
    subplot(1,5,iR)
    
    ylim([-12 20])
    pa = patch([50 -50 -50 50], [-11 -11 16 16], 'k'); hold on;
    pa.FaceColor = [245 204 85]/255;
    pa.FaceAlpha = 0.2;
    pa.LineStyle = '--';
    pa.LineWidth = 1;
    vl = vline(0); vl.LineWidth = 1; vl.LineStyle = '-';
    
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = mean(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan((ripple_LFP_region_all{iR}(:,1)))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bl1.LineWidth =1;

    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
%     text(-75, 19, sprintf('%i ripples', sum(~isnan((ripple_LFP_region_all{iR}(:,1)))))); hold on;
    xlim([-52 52])
    ylim([-12 20])
    ax = gca;
    ax.XTick = [-150, 150];
%     ylabel('amplitude [uV]')
%     xlabel('time from ripple center')
    sum(~isnan((ripple_LFP_region_all{iR}(:,1))))
%     title(regions{iR})

    axis off

 
end
axis off

fig = gcf;
fig.Color = 'w';


savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_2zoom_%s.pdf', tag)))


figure('Position', [-64 1149 2353 424]);


for iR = 1:length(regions)
    subplot(1,5,iR)
    
    mu = mean(ripple_LFP_region_all{iR}, 'omitnan');
    sem = mean(ripple_LFP_region_all{iR}, 'omitnan') / sqrt(sum(~isnan((ripple_LFP_region_all{iR}(:,1)))));
    [bl1, bf] = boundedline(-500:500, mu, sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    xlim([-500 500])
    
    ylabel('amplitude [uV]')
    xlabel('time from ripple center')
    title(regions{iR})
 
end


fig = gcf;
fig.Color = 'w';


savepdf(gcf, fullfile(exportDirFigs, sprintf('LFP_%s.pdf', tag)))


figure('Position', [1211 1267 781 294]); 
subplot(1,3,1)
h = pie(cchRegion, regions);
% Extract text handles and set labels
patchHandles = findobj(h, 'Type', 'patch'); % Find text objects
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = regionColors(ii,:);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1,pos(2), num2str(cchRegion(ii))); hold on;
end
title('channel locations')

subplot(1,3,2)
h = pie(cUnitRegion, regions);
% Extract text handles and set labels
patchHandles = findobj(h, 'Type', 'patch'); % Find text objects
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = regionColors(ii,:);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1,pos(2), num2str(cUnitRegion(ii))); hold on;
end
title('unit locations')

subplot(1,3,3)
types = {'pyr', 'int'}; 
c = [86, 180, 233 ;
    204, 121, 167]/255;
h = pie([sum(strcmp(uType, 'pyr')) sum(strcmp(uType, 'int'))], types);
patchHandles = findobj(h, 'Type', 'patch'); % Find text objects
for ii = 1:length(patchHandles)
    patchHandles(ii).FaceColor = c(ii,:);
    patchHandles(ii).LineWidth = 1.5;
    pos = computePatchCentroid(patchHandles(ii).Vertices);
    text(pos(1)-0.1,pos(2), num2str(sum(strcmp(uType, types{ii})))); hold on;
end
title('unit types')

fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirFigs, 'LFP_Unit_Count.pdf'))

figure; 
pl = plot(troughToPeak(troughToPeak >= 45 & isfinite(burstIndex)), log10(burstIndex(troughToPeak >= 45 & isfinite(burstIndex))-min(burstIndex)), '.'); hold on;
pl.Color = c(1,:);
pl = plot(troughToPeak(troughToPeak < 45 & isfinite(burstIndex)), log10(burstIndex(troughToPeak < 45 & isfinite(burstIndex))-min(burstIndex)), '.'); hold on;
pl.Color = c(2,:);
ylabel('burstIndex'); xlabel('troughToPeak [samples]')

figure; 
pl = plot(troughToPeak(troughToPeak >= 45 & isfinite(burstIndex)), FR(troughToPeak >= 45 & isfinite(burstIndex)), '.'); hold on;
pl.Color = c(1,:);
pl = plot(troughToPeak(troughToPeak < 45 & isfinite(burstIndex)), FR(troughToPeak < 45 & isfinite(burstIndex)), '.'); hold on;
pl.Color = c(2,:);
ylabel('FR'); xlabel('troughToPeak [samples]')


figure; 
histogram(log10(burstIndex(troughToPeak > 45 & isfinite(burstIndex))-min(burstIndex)), -2:0.1:3, 'Normalization','probability'); hold on;
histogram(log10(burstIndex(troughToPeak < 45 & isfinite(burstIndex))-min(burstIndex)), -2:0.1:3, 'Normalization','probability')
%%
fs = rippleStats.fs;
% hand = figure('Position', [503.0000 555 367.0000 318.0000]);
hand = figure('Position', [503 555 round(367*430/528) 318]);
% for ii = 1:length(regions)
%     subplot(1,5,ii)
%     erspRegion = strcmp(locAll, regions{ii});
    erspMu = mean(ERSPall(:,erspTimes >= -1.0e3 & erspTimes < 1.0e3,:),3,'omitnan');
    im = imagesc([-1.0e3   1.0e3], [1 fs/2],erspMu);hold on
    ax = gca;
    ax.YDir = 'normal';
    ax=imgca;
    set(ax,'YScale', 'log')
    ax.LineWidth = 1;
    % xticks([xd(1) (2/3)*xd(1) (1/3)*xd(1) 0 (1/3)*xd(2) (2/3)*xd(2) xd(2)])
    yticks([1 2 4 10 16 50 70 100  200 fs/2])
%     yticks([1 2 4 10 16 60 120 200 fs/2])
%     caxis([-max(abs([prctile(erspMu(:),0.01) prctile(erspMu(:),99.99)])) max(abs([prctile(erspMu(:),0.01) prctile(erspMu(:),99.99)]))])
    caxis([-1.7 1.7])
    % caxis([0 max(abs([prctile(erspMu(:),0.01) prctile(erspMu(:),99.99)]))])
    % caxis([-3 3])
    % xlim([-fs*0.15 fs*0.950])
    ylim([1 fs/2])
    ylim([50 200])
    v  = vline(0, 'k:');
    v.LineWidth = 1.5;
    xlabel('time from ripple center [ms]')
    ylabel('frequency (Hz)')
    c=colorbar;
    c.Label.String = 'ERSP (dB)';
    c.Label.Interpreter = 'tex';
    c.Ticks = [-1.5 1.5];
    c.LineWidth = 1;
    set(gca, 'FontSize', 15)
    set(gcf, 'Color', [1 1 1])
    % hline([70 90 120])
%     colormap(fliplr(slanCM('prinsenvlag',61)))
    colormap(fliplr(slanCM('fusion',501)))
% end

savepdf(gcf, fullfile(exportDirFigs, 'LFP_ripple_ERSP_z.pdf'))

%%
close all
figure('Position', [21 997 809 251]);
subplot(1,4,1)
% h = boxplot(densityAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(densityAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR,:);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String = 'ripple density [/min]';
% ylim([quantile(densityAll, 0.01) quantile(densityAll, 0.99)]);
whiskers = findobj(gca, 'Tag', 'Outliers');
set(whiskers, 'LineWidth', 4);
box off
% cellfun(@(X) mean(densityAll(strcmp(locAll, X))), locs)
% cellfun(@(X) std(densityAll(strcmp(locAll, X))), locs)
subplot(1,4,2)
% h = boxplot(durAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(durAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR,:);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String = 'ripple duration [ms]';
ylim([50 120]);
box off
% cellfun(@(X) mean(durAll(strcmp(locAll, X))), locs)
% cellfun(@(X) std(durAll(strcmp(locAll, X))), locs)
subplot(1,4,3)
% h = boxplot(ampAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(ampAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR,:);
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String =  'ripple amplitude [uV]';
ylim([0 25]);
box off
% cellfun(@(X) mean(ampAll(strcmp(locAll, X))), locs)
% cellfun(@(X) std(ampAll(strcmp(locAll, X))), locs)
subplot(1,4,4)
% h = boxplot(freqAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(freqAll, locAll, 'GroupOrder', regions);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(regions)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinColor{1} = regionColors(iR,:);
    vp(ii).ViolinAlpha{1} = 0.8;
        vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String = 'ripple frequency [Hz]';
ylim([80 100]);
box off
cellfun(@(X) mean(freqAll(strcmp(locAll, X))), locs)
cellfun(@(X) std(freqAll(strcmp(locAll, X))), locs)


fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'rippleMetrics.pdf'))

figure('Position', [1 68 2056 1099]);
subplot(4,1,1)
% h = boxplot(densityAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(densityAll, subjAll);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;
for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end

ax.YLabel.String = 'ripple density [/min]';
% ylim([quantile(densityAll, 0.01) quantile(densityAll, 0.99)]);
whiskers = findobj(gca, 'Tag', 'Outliers');
set(whiskers, 'LineWidth', 4);
box off

subplot(4,1,2)
% h = boxplot(durAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(durAll, subjAll);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(vp)
    iR = find(strcmp(regions, locs{ii}));
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String = 'ripple duration [ms]';
ylim([50 120]);
box off

subplot(4,1,3)
% h = boxplot(ampAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(ampAll, subjAll);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String =  'ripple amplitude [uV]';
ylim([0 25]);
box off
% cellfun(@(X) mean(ampAll(strcmp(locAll, X))), locs)
% cellfun(@(X) std(ampAll(strcmp(locAll, X))), locs)
subplot(4,1,4)
% h = boxplot(freqAll, 'Orientation', 'vertical', 'Symbol', '');
vp = violinplot(freqAll, subjAll);
ax = gca;
locs = ax.XTickLabel;
ax.LineWidth = 1;
ax.FontSize = 11;

for ii = 1:length(vp)
    vp(ii).ShowData = 0;
    vp(ii).BoxWidth = 0.05;
    vp(ii).ViolinAlpha{1} = 0.8;
    vp(ii).ViolinPlot.LineWidth = 1;
    vp(ii).ViolinPlot.EdgeColor = 'k';
    vp(ii).MedianPlot.SizeData = 20;
    vp(ii).ShowBox = 0; %.125;
    vp(ii).ShowWhiskers = 0; %.125;
end
ax.YLabel.String = 'ripple frequency [Hz]';
ylim([80 100]);
box off
cellfun(@(X) mean(freqAll(strcmp(locAll, X))), locs)
cellfun(@(X) std(freqAll(strcmp(locAll, X))), locs)


fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'rippleMetrics_bySubj.pdf'))

%%
% close all
statsBinWidth = 25;  %histogram bin width (in ms) for computing all stats (significance, sidedness)
plotBinWidth  = 5; %histogram bin width (in ms) for plotting
fs            = 1000; %sample rate
% compute histogram parameters
% nIter = size(NC.subjPRTH.nullPRTH{1,2},1);
win = ceil(size(coR_PRTH,1)/2);
edges = -win+(statsBinWidth/2):statsBinWidth:win;
edgesPlot = -win+(plotBinWidth/2):plotBinWidth:win;
times = -win:win;
nBins = 2*win/statsBinWidth - 1;
nBinsPlot = 2*win/plotBinWidth - 1;
center = ceil(nBins/2);
middle = movmean(-win+(plotBinWidth/2):plotBinWidth:win,2);
middle(1) = [];
iRa = 5; 
iRapl = find(contains(regionsAll, regions{iRa}));
% iRapl = find(contains(regionsAll, regions([1:3])));
% iRbpl = find(contains(regionsAll, regions{5}));

iRbpl = find(contains(regionsAll, regions([1:3])));


figure; 
c = 1;
countsPlot = [];
countsPlotNull =[];
for Aloop = iRapl %1:length(regionsAll)
    for Bloop = iRbpl%iRa+1:length(regionsAll)
        if all(ismember(iRapl,iRbpl))
            iRa = Aloop;
            iRb = Bloop;
        else
            iRa = min([Aloop Bloop]);
            iRb = max([Aloop Bloop]);
        end
        if iRa == iRb; continue; end
        
        fprintf('%s %s\n', regionsAll{iRa}, regionsAll{iRb})
        countsTemp = squeeze(coR_PRTH(:,2,1:ccgCount(iRa,iRb)-1, iRa, iRb)); %./ccgNrip(1, 1:ccgCount(iRa,iRb)-1,iRa, iRb)/plotBinWidth*1e3; 
        countsTemp(:,ccgCount(iRa,iRb):end) = [];
%         countsTemp = smoothdata(countsTemp, 1, 'gaussian',25);
        
        temp = zeros(size(countsTemp,2),nBinsPlot);
        for e = 1:length(edgesPlot)-1
           ii = times > edgesPlot(e) & times < edgesPlot(e+1);
           temp(:,e) = sum(countsTemp(ii,:));
        end
        
%         countsTemp = coR_PRTH(:,1:ccgCount(iRa,iRb)-1, iRa, iRb)./ccgNrip(1:ccgCount(iRa,iRb)-1,iRa, iRb)'/plotBinWidth*1e3; countsTemp(:,ccgCount(iRa,iRb):end) = [];
%         countsTemp = smoothdata(countsTemp, 1, 'gaussian',25);
%         tempNull = zeros(size(countsTemp,2),nBinsPlot);
% 
%         for e = 1:length(edgesPlot)-1
%            ii = times > edgesPlot(e) & times < edgesPlot(e+1);
%            tempNull(:,e) = sum(countsTemp(ii,:));
%         end
%         countsPlotNull = [countsPlotNull; tempNull];
        countsPlot = [countsPlot; temp];

        
        c = c+1;
    end
end
fprintf('\n')

% countsPlot = countsPlot;
SEM = std(countsPlot,1,'omitnan') / sqrt(size(countsPlot,1));
boundedline(middle, sum(countsPlot,1,'omitnan'), SEM); hold on;
% bar(middle, mean(countsPlot,1)); hold on;
% ylim([0.9 1.2])
% SEM = std(countsPlotNull,1) / sqrt(size(countsPlotNull,1));
% plot(middle, mean(countsPlotNull,1), 'r-');
% hl = hline(1);
% hl.LineStyle = '-';
vline(0)

xlim([-500 500])

%%
clc

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
pdfExport = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));
fs = 1e3; %1e3, 400;
ctxParc = {'ACC', 'SMA', 'OFC' ,'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};

iRa = 4; 
iRapl = find(contains(regionsAll, regions{iRa}));
iRapl = find(contains(regionsAll, regions([1:5])));
% iRbpl = find(contains(regionsAll, regions{5}));

iRbpl = find(contains(regionsAll, regions([1:5])));


figure; 
c = 1;
coRrate = [];
x = [];
yRip = [];
yShuff = [];
hemA = [];
hemB = [];
withinBundle = [];
for iRa = 1:length(regionsAll)
    for iRb = iRa:length(regionsAll)
       
        fprintf('%s %s\n', regionsAll{iRa}, regionsAll{iRb})
        co = coRnum(1:ccgCount(iRa,iRb)-1, iRa, iRb);
%         shuf = coRnumShuff(1:ccgCount(iRa,iRb)-1, iRa, iRb);
        N = coRdur(1:ccgCount(iRa,iRb)-1, iRa, iRb);
        temp = co./N;
%         tempShuf = shuf./N;
%         tempShuf(isnan(tempShuf)) = 0;

        mu = mean(temp);
        sem = std(temp) / sqrt(length(temp));
        parcelA = regionsAll{iRa}(2:end);
        parcelB = regionsAll{iRb}(2:end);

        hemA = [hemA repmat({regionsAll{iRa}(1)}, [1 length(temp)])];
        hemB = [hemB repmat({regionsAll{iRb}(1)}, [1 length(temp)])];
        samebund = iRa == iRb;
        withinBundle = [withinBundle repmat(samebund, [1 length(temp)])];

        broadmanA = broadman{strcmp(ctxParc, parcelA)};
        HCPa = sprintf('%s_%s', hemA{end}, broadmanA);
        broadmanB = broadman{strcmp(ctxParc, parcelB)};
        HCPb = sprintf('%s_%s', hemB{end}, broadmanB);


        d = tractLenStruct.tractLengths(strcmp(tractLenStruct.parcelIDs, HCPa), ...
                                        strcmp(tractLenStruct.parcelIDs, HCPb));
%         errorbar(d, mu, sem,'o'); hold on;
%         plot(repmat(d, [1 length(temp)]), temp,'.'); hold on;
        x = [x repmat(d, [1 length(temp)])];
        yRip = [yRip temp'];
%         yShuff = [yShuff tempShuf'];
        
%         if any(tempShuf > 0.25); error(); end
        
%         x = [x d];
%         y = [y mu];
        
        c = c+1;
    end
end

clr = brewermap(10,'Paired');

yRipMu = [];
yShufMu = [];
yRipSEM = [];
yShufSEM = [];
xMu = [];
xBin = quantile(x,5);
xBin = [0 xBin];
xBin = [xBin max(x)];
for iB = 2:length(xBin)
    ii = x >= xBin(iB-1) & x <= xBin(iB);
    yRipMu = [yRipMu mean(yRip(ii))];
    yShufMu = [yShufMu mean(yShuff(ii))];
    yRipSEM = [yRipSEM std(yRip(ii))/sqrt(sum(ii))];
    yShufSEM = [yShufSEM std(yShuff(ii))/sqrt(sum(ii))];
    xMu = [xMu mean([xBin(iB-1) xBin(iB)])];
end

% plot(x, yRip,'.'); hold on;

[bl2, bf] = boundedline(xMu, yRipMu, yRipSEM, 'rs-'); hold on;
bl2.Color = clr(8,:);
bl2.MarkerFaceColor = clr(8,:);
bf.FaceColor = clr(8,:);
bf.FaceAlpha = 0.3;

[bl2, bf] = boundedline(xMu, yShufMu, yShufSEM, 'bs-'); hold on;
bl2.Color = clr(1,:);
bl2.MarkerFaceColor = clr(1,:);
bf.FaceColor = clr(1,:);
bf.FaceAlpha = 0.3;
fig = gcf;
fig.Color = 'w';
% ylim([0 0.1])

figure('Position', [100 368 round(357*(2/3)) round(155)]);

subplot(1,3,1)
groups = nan(1, length(yRip));
groups(withinBundle == 1) = 1;
groups(withinBundle == 0) = [];

bx = boxplot(yRip(withinBundle == 1), groups, 'notch', 'off', 'Symbol', '', 'width', 0.3);
ylim([0 0.45])
ax = gca;
ax.FontSize = 7;
ax.LineWidth = 1;
set(bx, 'LineWidth', 1);

% Find the box objects in the current axes (each box has the tag 'Box')
boxes = findobj(gca, 'Tag', 'Box');

% Define your desired colors
fillColor = [0.7, 0.9, 1];    % light blue fill
edgeColor = [0, 0, 0.8];       % dark blue edge

% Loop over each box object to overlay a patch
for i = 1:length(boxes)
    % Get the box coordinates
    xdata = get(boxes(i), 'XData');
    ydata = get(boxes(i), 'YData');
    
    % Create a patch to fill the box with your desired color.
    % FaceAlpha adjusts transparency; set to 1 for opaque.
    p = patch(xdata, ydata, fillColor, 'FaceAlpha', 0.5, ...
              'EdgeColor', edgeColor, 'LineWidth', 1.0);
    
    % (Optional) Bring the original box lines to the top so they aren?t hidden.
    uistack(boxes(i), 'top');
end
box off
xlim([0.5 1.5])
subplot(1,3,[2 3])
groups = nan(1, length(yRip));
groups(strcmp(hemA,hemB) & withinBundle == 0) = 1;
groups(~strcmp(hemA,hemB) & withinBundle == 0) = 2;
groups(withinBundle == 1) = [];

bx = boxplot(yRip(withinBundle == 0), groups,'notch', 'off','Symbol', '');
ax = gca;
ax.FontSize = 7;
ax.LineWidth = 1;
set(bx, 'LineWidth', 1);

% Find the box objects in the current axes (each box has the tag 'Box')
boxes = findobj(gca, 'Tag', 'Box');

% Define your desired colors
fillColor = [0.7, 0.9, 1];    % light blue fill
edgeColor = [0, 0, 0.8];       % dark blue edge

% Loop over each box object to overlay a patch
for i = 1:length(boxes)
    % Get the box coordinates
    xdata = get(boxes(i), 'XData');
    ydata = get(boxes(i), 'YData');
    
    % Create a patch to fill the box with your desired color.
    % FaceAlpha adjusts transparency; set to 1 for opaque.
    p = patch(xdata, ydata, fillColor, 'FaceAlpha', 0.5, ...
              'EdgeColor', edgeColor, 'LineWidth', 1.0);
    
    % (Optional) Bring the original box lines to the top so they aren?t hidden.
    uistack(boxes(i), 'top');
end
% ax.YLabel.String = 'ripple duration [ms]';
ylim([0 0.1]);
box off

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'coRipRates.pdf'))

