

close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;


matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));

LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '../task/*task*'));
taskfiles = {taskfiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../../out', '*micro*'));
micrFiles = {micrFiles.name}';

recordingState = 'wake';
location = 'NC';
% tag = 'vHG';
sfreq = 1e3;
copmuteTF = false; 
% regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', 'LSPE', ...
%            'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP', 'RSPE'};
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end,3:end-1));   
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2,3:end-1));
bundleNotes = table2cell(bundleNotes(3:end,3:end-1));   


%%

computeCoRipple = false;
rippleFR = nan(length(subj_list_full), 100);
CoRippleRatesRegionAll = nan(length(regions), length(regions), length(subj_list_full));
win = [-9 25] * 1e3;
winResp = [-10.5 1] * 1e3;
winLength = 12905; %ms

%prealocate variables
recog_trials_all = nan(6.5e3, length(win(1):win(2)));
recog_trials_units_all = nan(6.5e3, length(win(1):win(2)));

recog_trials_resp_all = nan(6.5e3, length(winResp(1):winResp(2)));
recog_trials_resp_avg_all = nan(6.5e3,1);
recog_trials_units_resp_all = nan(6.5e3, length(winResp(1):winResp(2)));
whole_trial_all = nan(6.5e3, winLength);
whole_trial_units_all = nan(6.5e3, winLength);
whole_trials_region_all = nan(6.5e3, winLength, length(regions));
resp_region_all = nan(6.5e3, length(winResp(1):winResp(2)), length(regions));
resp_units_region_all = nan(6.5e3, length(winResp(1):winResp(2)), length(regions));
whole_trials_region_all_shuff = nan(6.5e3, winLength, length(regions));
whole_trials_ERP_region_all  = nan(6.5e3, winLength, length(regions));
whole_trials_BP_region_all  = nan(6.5e3, winLength, length(regions), 3);
whole_trials_units_sm_region_all  = nan(6.5e3, winLength, length(regions));
whole_trials_units_region_all  = nan(6.5e3, winLength, length(regions));
whole_trials_cross_region_all = cell(length(regions));
whole_trials_cross_region_all_shuff = cell(length(regions));
probe_trials_units_all = cell(3, length(regions));
probe_trials_all = cell(3, length(regions));
nChRegions = zeros(6.5e3, length(regions));
countChPerRegion = zeros(1, length(regions));
countUnitPerRegion = zeros(1, length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        whole_trials_cross_region_all{iRa, iRb} = nan(6.5e3, winLength);       
        whole_trials_cross_region_all_shuff{iRa, iRb} = nan(6.5e3, winLength);       
        for ii = 1:size(whole_trials_unit_cross_region_all,3)
            whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(6.5e3, winLength);
        end
    end
end

subjID = nan(6.5e3, 1);
correct_trials_all = nan(6.5e3, 1);
load_trials_all = nan(6.5e3, 1);
probe_in_out = nan(6.5e3, 1);
respLatency = nan(6.5e3, 1);
spikePerBundleRip = nan(6.5e3, 1);
densityAll = [];
freqAll = [];
durAll = [];
ampAll = [];
uChanAll = [];
tic
cTrial = 1;
countBundU = 1;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};

    fprintf('processing %s ... \n', subject)



    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    LFPtimeMicr = micrObj.times;
%     chan_labels = strsplit(micrObj.chan_locations, micrObj.chan_locations(38))
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    
    
    
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@(X) isempty(X), locations)) = [];
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
   
    
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locations(hemi) = cellfun(@(X) ['L' X], locations(hemi),  'UniformOutput', false);
    hemi = contains(hem, 'right');
    locations(hemi) = cellfun(@(X) ['R' X], locations(hemi),  'UniformOutput', false);
    
    

    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, '../task/',taskfiles{f}));

    
    modifier = '1kHz_template_z25';

    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    LFPchanNum = cellfun(@(X) str2double(X), rippleStats.chanLabels);

    rippleStats.chanLabels = locations;
    
    if length(rippleStats.locs) ~= length(locations); error('error formatting channel labels.'); end
        
    rippleStats.recordingType = [repmat({'micro'}, [1 length(microObj.rippleStats.locs)])];
    rippleStats.locs = [cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false)];
    microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = [cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', microObj.rippleStats.window, 'UniformOutput', false)];

    rippleStats.microTimes = LFPtimeMicr;
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    rippMaskShuff = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    for chRipp = 1:size(rippMask,1) %sum(strcmp(rippleStats.recordingType, 'macro'))+1:size(rippMask,1) 
        if rippleStats.density{chRipp} > 1
            

            if contains(rippleStats.chanLabels{chRipp}, 'SPE')
                rippMask(chRipp,:) = nan(1,length(rippMask));
                continue
            end
            

            densityAll = [densityAll rippleStats.density{chRipp}];
            freqAll = [freqAll mean(rippleStats.oscFreq{chRipp})];
            durAll = [durAll mean(rippleStats.duration{chRipp})];
            ampAll = [ampAll mean(rippleStats.rippleAmp{chRipp})];
            
            
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);        
        
            for ii = 1:length(iE)
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
            rippMaskShuff(chRipp,:) = shuffleMask(rippMask(chRipp,:));

        end
        
    end
    
    rippAll = sum(rippMask);

    rippMask(:, end:end+3e3) = nan; %pad the ending
    rippAll(end:end+3e3) = nan; %pad the ending
    rippMaskShuff(:, end:end+3e3) = nan; %pad the ending  
    
    data = micrObj.lfp_data;
    data(:, end:end+3e3) = nan; %pad the ending

    units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    
    uLocations = units(U,end-2);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
    hemi = contains(uLocations, 'left');
    uLocations(hemi) = cellfun(@(X) ['L' X], uLocations(hemi),  'UniformOutput', false);
    hemi = contains(uLocations, 'right');
    uLocations(hemi) = cellfun(@(X) ['R' X], uLocations(hemi),  'UniformOutput', false);
    uLocations = strrep(uLocations, '_left', '');
    uLocations = strrep(uLocations, '_right', '');
        
        
    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001; %seconds
    
    binEdges = 0:binWidth:(length(rippMask)/1e3);
    spikeMask = nan(nNeuron, length(binEdges)-1);
    spikeMask_sm = nan(nNeuron, length(binEdges)-1);
    for iU = 1:nNeuron 
        X = units{unitsAll(iU),2};
        [N,EDGES] = histcounts(X,binEdges);
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian',100);
    end
    
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian',100);
    else
        unitAll = smoothdata(spikeMask, 'gaussian',100);
    end
    unitAll = zscore(unitAll);
    
    baselineMask = false(1, length(rippMask));

    for iT = 1:length(trials.start_time)
            trialLoad = trials.loads(iT);
            image1Time =  round(trials.timestamps_Encoding1(iT)*1e3);
            baselineMask(image1Time-1000:image1Time-1) = true; %baseline
    end
    
    

    unit_trials_region = nan(length(trials.start_time), 12905, length(regions));
    dat_region = cell(1,length(regions));
    dat_region_shuff = cell(1,length(regions));
    erp_dat_region = cell(1,length(regions));
    u_dat_region_sm = cell(1,length(regions));
    u_dat_region = cell(1,length(regions));
    resp_per_chan = nan(size(rippMask,1), 140);
    for iT = 1:length(trials.start_time)
        respTime = round(trials.timestamps_Response(iT)*1e3);
        probeTime = round(trials.timestamps_Probe(iT)*1e3);
        respLatency(cTrial)= respTime-probeTime;

        %overal co-rippling
        trialData = nan(1, win(2)-win(1)+1);
        dat = rippAll(probeTime+win(1):respTime);
        if length(dat) > length(trialData)
            trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;
        end
        recog_trials_all(cTrial, :) = trialData;

        dat = sum(spikeMask_sm(:,respTime+winResp(1):respTime+winResp(2)),1);% / size(spikeMask_sm, 1);
        recog_trials_units_resp_all(cTrial, :) = dat;
        
        if sum(dat) == 0; continue; end

        
        trialLoad = trials.loads(iT);
        image1Time =  round(trials.timestamps_Encoding1(iT)*1e3);
        image2Time =  round(trials.timestamps_Encoding2(iT)*1e3);
        image3Time =  round(trials.timestamps_Encoding3(iT)*1e3);
        maintenanceTime =  round(trials.timestamps_Maintenance(iT)*1e3);
        probeTime =  round(trials.timestamps_Probe(iT)*1e3);
        
        dat = rippAll(respTime+winResp(1):respTime+winResp(2));
        recog_trials_resp_all(cTrial, :) = dat;
        
        dat = rippMask(:, respTime-2000:respTime);
        dat = mean(dat,2) ./ mean(rippMask(:, baselineMask),2);
        
        resp_per_chan(:,iT) = dat;

        if trialLoad == 3
            im1seg = rippAll(image1Time-1000:image1Time+2e3);
            im2seg = rippAll(image2Time-0  :image2Time+2e3);
            im3seg = rippAll(image3Time-0  :image3Time+2e3);
            maintenanceSeg = rippAll(maintenanceTime-0:maintenanceTime+1250);
            probeSeg = rippAll(probeTime-1250:probeTime+3000);
            
            
            trialFull = [im1seg, nan(1,100), ...
                         im2seg, nan(1,100), ...
                         im3seg, nan(1,100), ...
                         maintenanceSeg, nan(1,100), ...
                         probeSeg];

            nCh =  size(rippMask,1);
            im1seg = rippMask(:, image1Time-1000:image1Time+2e3);
            im2seg = rippMask(:, image2Time-0  :image2Time+2e3);
            im3seg = rippMask(:, image3Time-0  :image3Time+2e3);
            maintenanceSeg = rippMask(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = rippMask(:, probeTime-1250:probeTime+3000);
            
            
            ripTrialFullChan = [im1seg, nan(nCh,100), ...
                         im2seg, nan(nCh,100), ...
                         im3seg, nan(nCh,100), ...
                         maintenanceSeg, nan(nCh,100), ...
                         probeSeg];
                     
            nUnit =  size(spikeMask,1);
            im1seg = spikeMask_sm(:, image1Time-1000:image1Time+2e3);
            im2seg = spikeMask_sm(:, image2Time-0  :image2Time+2e3);
            im3seg = spikeMask_sm(:, image3Time-0  :image3Time+2e3);
            maintenanceSeg = spikeMask_sm(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = spikeMask_sm(:, probeTime-1250:probeTime+3000);
            
            
            unitTrialFullChan = [im1seg, nan(nUnit,100), ...
                         im2seg, nan(nUnit,100), ...
                         im3seg, nan(nUnit,100), ...
                         maintenanceSeg, nan(nUnit,100), ...
                         probeSeg];

        elseif trialLoad == 1
            im1seg = rippAll(image1Time-1000:image1Time+2e3);
            im2seg = nan(1, 2e3+0+1);
            im3seg = nan(1, 2e3+0+1);
            maintenanceSeg = rippAll(maintenanceTime-0:maintenanceTime+1250);
            probeSeg = rippAll(probeTime-1250:probeTime+3000);
            
            trialFull = [im1seg, nan(1,100), ...
                         im2seg, nan(1,100), ...
                         im3seg, nan(1,100), ...
                         maintenanceSeg, nan(1,100), ...
                         probeSeg];
                     
            nCh =  size(rippMask,1);
            im1seg = rippMask(:,image1Time-1000:image1Time+2e3);
            im2seg = nan(nCh, 2e3+0+1);
            im3seg = nan(nCh, 2e3+0+1);
            maintenanceSeg = rippMask(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = rippMask(:, probeTime-1250:probeTime+3000);
            
            
            ripTrialFullChan = [im1seg, nan(nCh,100), ...
                         im2seg, nan(nCh,100), ...
                         im3seg, nan(nCh,100), ...
                         maintenanceSeg, nan(nCh,100), ...
                         probeSeg];
                     
            nUnit =  size(spikeMask,1);
            im1seg = spikeMask_sm(:,image1Time-1000:image1Time+2e3);
            im2seg = nan(nUnit, 2e3+0+1);
            im3seg = nan(nUnit, 2e3+0+1);
            maintenanceSeg = spikeMask_sm(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = spikeMask_sm(:, probeTime-1250:probeTime+3000);
            
            
            unitTrialFullChan = [im1seg, nan(nUnit,100), ...
                         im2seg, nan(nUnit,100), ...
                         im3seg, nan(nUnit,100), ...
                         maintenanceSeg, nan(nUnit,100), ...
                         probeSeg];
        end
        
        whole_trial_all(cTrial,:) = trialFull; 
        whole_trial_units_all(cTrial,:) = sum(unitTrialFullChan, 1) / size(unitTrialFullChan,1);
        





        %within region co-rippling and unit firing
        for iR = 1:length(regions)
            chRegion = contains(rippleStats.chanLabels, regions(iR));
            uRegion = contains(uLocations, regions(iR));
            nChRegions(cTrial, iR) = sum(chRegion);
            if isempty(dat_region{iR}) && sum(chRegion) > 0
                 dat_region{iR} = sum(rippMask(chRegion,:),1);
                 dat_region_shuff{iR} = sum(rippMaskShuff(chRegion,:),1);
                 erp_dat_region{iR} = mean(data(chRegion,:),1);
                
                 bundles = unique(rippleStats.chanLabels(chRegion));
                 for iBund = 1:length(bundles)
                     
                     chBundle = contains(rippleStats.chanLabels, bundles(iBund));
                     uBundle  = find(contains(uLocations, bundles(iBund)));
                     rippBundleMask = sum(rippMask(chBundle,:)) > 0;
                     rippBnds = mask2bounds(rippBundleMask);

                     countSpikeRipp = zeros(size(spikeMask,1),1);
                
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

            if isempty(u_dat_region_sm{iR}) && sum(uRegion) > 0
                 u_dat_region_sm{iR} = sum(spikeMask_sm(uRegion,:),1)/ sum(uRegion);
                 u_dat_region{iR} = sum(spikeMask(uRegion,:),1) / sum(uRegion);
            elseif sum(uRegion) == 0
                 u_dat_region_sm{iR} = nan(1, length(rippMask));
                 u_dat_region{iR} = nan(1, length(rippMask));
            end
            
            
            
            
            if trialLoad == 3
                im1seg = dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = dat_region{iR}(image2Time-0  :image2Time+2e3);
                im3seg = dat_region{iR}(image3Time-0  :image3Time+2e3);
                maintenanceSeg = dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = dat_region{iR}(probeTime-1250:probeTime+3000);


                trialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];
                         
                im1seg = dat_region_shuff{iR}(image1Time-1000:image1Time+2e3);
                im2seg = dat_region_shuff{iR}(image2Time-0  :image2Time+2e3);
                im3seg = dat_region_shuff{iR}(image3Time-0  :image3Time+2e3);
                maintenanceSeg = dat_region_shuff{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = dat_region_shuff{iR}(probeTime-1250:probeTime+3000);


                trialFullShuff = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];
                         
                im1seg = erp_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = erp_dat_region{iR}(image2Time-0  :image2Time+2e3);
                im3seg = erp_dat_region{iR}(image3Time-0  :image3Time+2e3);
                maintenanceSeg = erp_dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = erp_dat_region{iR}(probeTime-1250:probeTime+3000);


                ERPtrialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];
    
                im1seg = u_dat_region_sm{iR}(image1Time-1000:image1Time+2e3);
                im2seg = u_dat_region_sm{iR}(image2Time-0  :image2Time+2e3);
                im3seg = u_dat_region_sm{iR}(image3Time-0  :image3Time+2e3);
                maintenanceSeg = u_dat_region_sm{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = u_dat_region_sm{iR}(probeTime-1250:probeTime+3000);


                uTrialFullsm = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];

                im1seg = u_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = u_dat_region{iR}(image2Time-0  :image2Time+2e3);
                im3seg = u_dat_region{iR}(image3Time-0  :image3Time+2e3);
                maintenanceSeg = u_dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = u_dat_region{iR}(probeTime-1250:probeTime+3000);


                uTrialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];

            elseif trialLoad == 1
                im1seg = dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = dat_region{iR}(probeTime-1250:probeTime+3000);


                trialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];
                         
                im1seg = dat_region_shuff{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = dat_region_shuff{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = dat_region_shuff{iR}(probeTime-1250:probeTime+3000);


                trialFullShuff = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];
                         
                im1seg = erp_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = erp_dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = erp_dat_region{iR}(probeTime-1250:probeTime+3000);


                ERPtrialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];

                im1seg = u_dat_region_sm{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = u_dat_region_sm{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = u_dat_region_sm{iR}(probeTime-1250:probeTime+3000);


                uTrialFullsm = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];

                im1seg = u_dat_region{iR}(image1Time-1000:image1Time+2e3);
                im2seg = nan(1, 2e3+0+1);
                im3seg = nan(1, 2e3+0+1);
                maintenanceSeg = u_dat_region{iR}(maintenanceTime-0:maintenanceTime+1250);
                probeSeg = u_dat_region{iR}(probeTime-1250:probeTime+3000);


                uTrialFull = [im1seg, nan(1,100), ...
                             im2seg, nan(1,100), ...
                             im3seg, nan(1,100), ...
                             maintenanceSeg, nan(1,100), ...
                             probeSeg];


             end

            unit_trials_region(iT,:,iR) = uTrialFull;           
            whole_trials_region_all(cTrial, :, iR) = trialFull;
            resp_region_all(cTrial, :, iR) = dat_region{iR}(respTime+winResp(1):respTime+winResp(2));
            resp_units_region_all(cTrial, :, iR) = u_dat_region_sm{iR}(respTime+winResp(1):respTime+winResp(2));
            whole_trials_region_all_shuff(cTrial, :, iR) = trialFullShuff;
            whole_trials_ERP_region_all(cTrial, :, iR) = ERPtrialFull;
            whole_trials_units_sm_region_all(cTrial, :, iR) = uTrialFullsm;
            whole_trials_units_region_all(cTrial, :, iR) = uTrialFull;
            
            
            uRegion = find(uRegion); %convert to index
            for ii = 1:length(uRegion)
                baslnUnits = mean(spikeMask(uRegion(ii),image1Time-1e3:image1Time));
                maintUnits = mean(spikeMask(uRegion(ii),maintenanceTime:maintenanceTime+2e3));
                probeUnits = mean(spikeMask(uRegion(ii),probeTime:probeTime+1250));
                probe_trials_units_all{trialLoad,iR} = [probe_trials_units_all{trialLoad,iR}; [countUnitPerRegion(iR)+ii baslnUnits maintUnits probeUnits]];
            end
            
            chRegion = find(chRegion); %convert to index
            for ii = 1:length(chRegion)
                baslnRipples = mean(rippMask(chRegion(ii),image1Time-1e3:image1Time));
                maintRipples = mean(rippMask(chRegion(ii),maintenanceTime:maintenanceTime+2e3));
                probeRipples = mean(rippMask(chRegion(ii),probeTime:probeTime+1250));
                probe_trials_all{trialLoad,iR} = [probe_trials_all{trialLoad,iR}; [countChPerRegion(iR)+ii baslnRipples maintRipples probeRipples]];
            end

        end
        
         %across region co-rippling
        for iRa = 1:length(regions)
            for iRb = iRa:length(regions)

                datArip =  whole_trials_region_all(cTrial,:,iRa);
                datBrip =  whole_trials_region_all(cTrial,:,iRb);
                
                datAripShuff =  whole_trials_region_all_shuff(cTrial,:,iRa);
                datBripShuff =  whole_trials_region_all_shuff(cTrial,:,iRb);
                
                if iRa == iRb
                    datAB = datArip;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                    whole_trials_cross_region_all_shuff{iRa, iRb}(cTrial, :) = datAripShuff;
                    whole_trials_unit_cross_region_all{iRa, iRb, 3}(cTrial, :) = unit_trials_region(iT,:,iR);
                    continue
                else

                    datAB = datArip + datBrip;
                    datAB(datArip <= 0 | datBrip <= 0) = 0;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                    
                    datAB = datAripShuff + datBripShuff;
                    datAB(datAripShuff <= 0 | datBripShuff <= 0) = 0;
                    whole_trials_cross_region_all_shuff{iRa, iRb}(cTrial, :) = datAB;

                end
                
                
                uRegionA = uChan(contains(uLocations, regions(iRa)));
                uRegionB = uChan(contains(uLocations, regions(iRb)));
                
                ovA = nan(1, length(uRegionA));
                ovAcell = arrayfun(@(X) find(LFPchanNum == X), uRegionA', 'UniformOutput', false);
                ovA(cellfun(@(X) ~isempty(X), ovAcell)) = cell2mat(ovAcell);
                ovB = nan(1, length(uRegionB));
                ovBcell = arrayfun(@(X) find(LFPchanNum == X), uRegionB', 'UniformOutput', false);
                ovB(cellfun(@(X) ~isempty(X), ovBcell)) = cell2mat(ovBcell);                
                if sum(~isnan(ovA)) < 1 || sum(~isnan(ovB)) < 1; continue; end

% 
%                 datArip =  ripTrialFullChan(ovA(~isnan(ovA)),:);
%                 datBrip =  ripTrialFullChan(ovB(~isnan(ovB)),:);
                
                

                
                                                                 
               
            end
        end

        trialData = nan(1, win(2)-win(1)+1);
        dat = unitAll(probeTime+win(1):respTime);
        if length(dat) > length(trialData)
            trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;
        end
        recog_trials_units_all(cTrial,:)= trialData;
        
        

        correct_trials_all(cTrial) = trials.response_accuracy(iT);
        load_trials_all(cTrial) = trials.loads(iT);

        probe_in_out(cTrial) = trials.probe_in_out(iT);
        subjID(cTrial) = subj;

        cTrial = cTrial + 1;
        
        
    end
    
    
    
    for iR = 1:length(regions)
            chRegion = contains(rippleStats.chanLabels, regions(iR));
            countChPerRegion(iR) = countChPerRegion(iR) + sum(chRegion);
            
            uRegion = contains(uLocations, regions(iR));
            countUnitPerRegion(iR) = countUnitPerRegion(iR) + sum(uRegion);
    end
    
    fprintf('number of total trials: %i, number of subject trials: %i \n\n', cTrial - 1, length(trials.start_time))
 
end
correct_trials_all(cTrial:end) = [];
load_trials_all(cTrial:end) =[];

probe_in_out(cTrial:end) = [];
subjID(cTrial:end) = [];

fprintf('done processing .. \n')
toc


%%
close all

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

            baselineRipp = ripTrial(1:1e3) > 0;
            baselineUnit = uTrial(1:1e3);

            baseSpike = baseSpike + sum(baselineUnit(baselineRipp), 'omitnan');
            baseRipDur = baseRipDur + sum(baselineRipp) / 1e3;

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

[P,H,STATS] = signrank((taskSpikePerBund ./ taskRipDurPerBund)  - (baseSpikePerBund ./ baseRipDurPerBund), 'tail', 'right')

%%
regionColors =  brewermap(12, 'Dark2');
plot_condition = load_trials_all == 3;
taskMarkers = [500 2600 4700 6800 9400]+500;
coR = nan(length(plot_condition), length(regions)+1, 6);
FR = nan(length(plot_condition), length(regions)+1, 6);
pLoadCoR = nan(length(regions)+1,6);
pLoadFR  = nan(length(regions)+1,6);
pBslnCoR = nan(length(regions)+1,6,2);
pBslnFR = nan(length(regions)+1,6,2);
for iR = 1:length(regions)+1
    for iT = 1:length(plot_condition)
        if iR == 6
            if plot_condition(iT)

                coR(iT, iR, 1) = std(whole_trial_all(iT, 1:taskMarkers(1)), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trial_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                coR(iT, iR, 3) = mean(whole_trial_all(iT, taskMarkers(2):taskMarkers(3)), 'omitnan');
                coR(iT, iR, 4) = mean(whole_trial_all(iT, taskMarkers(3):taskMarkers(4)), 'omitnan');
                coR(iT, iR, 5) = mean(whole_trial_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trial_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');

                FR(iT, iR, 1) = mean(whole_trial_units_all(iT, 1:taskMarkers(1)), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trial_units_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                FR(iT, iR, 3) = mean(whole_trial_units_all(iT, taskMarkers(2):taskMarkers(3)), 'omitnan');
                FR(iT, iR, 4) = mean(whole_trial_units_all(iT, taskMarkers(3):taskMarkers(4)), 'omitnan');
                FR(iT, iR, 5) = mean(whole_trial_units_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trial_units_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');

            else
                coR(iT, iR, 1) = std(whole_trial_all(iT, 1:taskMarkers(1)), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trial_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                coR(iT, iR, 3) = nan;
                coR(iT, iR, 4) = nan;
                coR(iT, iR, 5) = mean(whole_trial_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trial_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');

                FR(iT, iR, 1) = mean(whole_trial_units_all(iT, 1:taskMarkers(1)), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trial_units_all(iT, taskMarkers(1):taskMarkers(2)), 'omitnan');
                FR(iT, iR, 3) = nan;
                FR(iT, iR, 4) = nan;
                FR(iT, iR, 5) = mean(whole_trial_units_all(iT, taskMarkers(4):taskMarkers(5)), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trial_units_all(iT, taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');
                
            end
        else
            if plot_condition(iT)
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
            else
                coR(iT, iR, 1) = mean(whole_trials_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                coR(iT, iR, 2) = mean(whole_trials_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                coR(iT, iR, 3) = nan;
                coR(iT, iR, 4) = nan;
                coR(iT, iR, 5) = mean(whole_trials_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                coR(iT, iR, 6) = mean(whole_trials_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');

                FR(iT, iR, 1) = mean(whole_trials_units_sm_region_all(iT, 1:taskMarkers(1), iR), 'omitnan');
                FR(iT, iR, 2) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(1):taskMarkers(2), iR), 'omitnan');
                FR(iT, iR, 3) = nan;
                FR(iT, iR, 4) = nan;
                FR(iT, iR, 5) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(4):taskMarkers(5), iR), 'omitnan');
                FR(iT, iR, 6) = mean(whole_trials_units_sm_region_all(iT, taskMarkers(5):taskMarkers(5)+1e3, iR), 'omitnan');
            end
        end
    end
    
    for iTm = [1,2,5,6]
        
        % load 3 vs load 1
        dat = coR(:, iR, iTm);
        tab = table(subjID(~isnan(dat)), plot_condition(~isnan(dat)), dat(~isnan(dat)), ...
          'VariableNames', {'subject',   'condition',                 'rates'});
        lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
        pLoadCoR(iR, iTm) = lme.Coefficients(2,6);


        dat = FR(:, iR, iTm);
        tab = table(subjID(~isnan(dat)), plot_condition(~isnan(dat)), dat(~isnan(dat)), ...
          'VariableNames', {'subject',   'condition',                 'rates'});
        lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
        pLoadFR(iR, iTm) = lme.Coefficients(2,6);


        % load vs baseline
        dat = coR(plot_condition, iR, iTm);
        tab = table(subjID(plot_condition), dat,  coR(plot_condition, iR, 1), ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
        pBslnCoR(iR, iTm,2) = ranksum( coR(plot_condition, iR, iTm),  coR(:, iR, 1));
        
        dat = coR(~plot_condition, iR, iTm);
        tab = table(subjID(~plot_condition), dat,  coR(~plot_condition, iR, iTm),  ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
        pBslnCoR(iR, iTm,1) = ranksum( coR(~plot_condition, iR, iTm),  coR(:, iR, 1));
        
        if iTm == 2
            dat = coR(plot_condition, iR, iTm+1);
            tab = table(subjID(plot_condition), dat,  coR(plot_condition, iR, iTm+1),  ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
            pBslnCoR(iR, iTm+1,2) = ranksum( coR(plot_condition, iR, iTm+1),  coR(:, iR, 1));
        
            dat = coR(plot_condition, iR, iTm+1);
            tab = table(subjID(plot_condition), dat,  coR(plot_condition, iR, iTm), ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
            pBslnCoR(iR, iTm+2,2) = ranksum( coR(plot_condition, iR, iTm+2),  coR(:, iR, 1));

            
        end
        
        dat = FR(plot_condition, iR, iTm);
        tab = table(subjID(plot_condition), dat,  FR(plot_condition, iR, 1), ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
        pBslnFR(iR, iTm,2) = ranksum( FR(plot_condition, iR, iTm),  FR(:, iR, 1));
        
        dat = FR(~plot_condition, iR, iTm);
        tab = table(subjID(~plot_condition), dat,  FR(~plot_condition, iR, iTm),  ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
        pBslnFR(iR, iTm,1) = ranksum( FR(~plot_condition, iR, iTm),  FR(:, iR, 1));
        
        if iTm == 2
            dat = FR(plot_condition, iR, iTm+1);
            tab = table(subjID(plot_condition), dat,  FR(plot_condition, iR, iTm+1),  ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
            pBslnFR(iR, iTm+1,2) = ranksum( FR(plot_condition, iR, iTm+1),  FR(:, iR, 1));
        
            dat = FR(plot_condition, iR, iTm+1);
            tab = table(subjID(plot_condition), dat,  FR(plot_condition, iR, iTm), ...
          'VariableNames', {'subject',  'rates', 'rates_baseline'});
            pBslnFR(iR, iTm+2,2) = ranksum( FR(plot_condition, iR, iTm+2),  FR(:, iR, 1));

            
        end
            
           
        
    end    
    
    
    
end
FR = FR *1e3;
pLoadCoR_adj = nan(size(pLoadCoR));
[h, crit_p, adj_ci_cvrg, pLoadCoR_adj(~isnan(pLoadCoR))]=fdr_bh(pLoadCoR(~isnan(pLoadCoR)),0.05,'pdep','yes');
pLoadFR_adj = nan(size(pLoadCoR));
[h, crit_p, adj_ci_cvrg, pLoadFR_adj(~isnan(pLoadFR))]=fdr_bh(pLoadFR(~isnan(pLoadFR)),0.05,'pdep','yes');

pBslnCoR_adj = nan(size(pBslnCoR));
[h, crit_p, adj_ci_cvrg, pBslnCoR_adj(~isnan(pBslnCoR))]=fdr_bh(pBslnCoR(~isnan(pBslnCoR)),0.05,'pdep','yes');
pBslnFR_adj = nan(size(pBslnFR));
[h, crit_p, adj_ci_cvrg, pBslnFR_adj(~isnan(pBslnFR))]=fdr_bh(pBslnFR(~isnan(pBslnFR)),0.05,'pdep','yes');

figure('Position',[669 805 285 round(729*(6/5))]);
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
    h = daboxplot(dat,'groups',group_num, 'mean',1,'outlier', 0, 'color',cBox, 'xshift', 0, 'boxalpha', 0.2, 'boxwidth', 1.3, 'meansize', 4, 'mediansize', 2); hold on;
    ax = gca;
    ax.XTick = [1:6];
    ax.XTickLabel = {'B','E1','E2','E3','M','P'};
    hl = hline(median(coR(:,iR,1), 'all', 'omitnan'));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    
    pVals = find(pBslnCoR_adj(iR,:,1) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(h.gpos(1, iTm), y_max,repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnCoR_adj(iR,iTm,1), 'maxNum',1)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
    pVals = find(pBslnCoR_adj(iR,:,2) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;

        text(mean(h.gpos(2, iTm)), y_max,repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnCoR_adj(iR,iTm,2), 'maxNum',1)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
    pVals = find(pLoadCoR_adj(iR,:) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.3/1.27 * y_max;
        line([h.gpos(1, iTm) h.gpos(2, iTm)], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5); hold on
        text(mean(h.gpos(:, iTm)), y_max + 0.05/1.27 * y_max,repmat('*', 1, countLeadingZerosAfterDecimal(pLoadCoR_adj(iR,iTm), 'maxNum',3)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegionsBoxPlot_%s.pdf', tag)))

figure('Position',[669 805 285 round(729*(6/5))]);
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
    h = daboxplot(dat,'groups',group_num, 'mean',1,'outlier', 0, 'color',cBox, 'xshift', 0, 'boxalpha', 0.2, 'boxwidth', 1.3, 'meansize', 4, 'mediansize', 2); hold on;
    ax = gca;
    ax.XTick = [1:6];
    ax.XTickLabel = {'B','E1','E2','E3','M','P'};
    hl = hline(median(FR(:,iR,1), 'all', 'omitnan'));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    pVals = find(pBslnFR_adj(iR,:,1) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(h.gpos(1, iTm), y_max,repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnFR_adj(iR,iTm,1), 'maxNum',1)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
    pVals = find(pBslnFR_adj(iR,:,2) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]);
        y_max = y_max + 0.15/1.27 * y_max;
        text(mean(h.gpos(2, iTm)), y_max,repmat(char(8224), 1, countLeadingZerosAfterDecimal(pBslnFR_adj(iR,iTm,2), 'maxNum',1)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
    pVals = find(pLoadFR_adj(iR,:) < 0.05);
    
    for iP = 1:length(pVals)
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,1,3).YData, h.wh(iTm,2,3).YData]) + 0.30;
        y_max = y_max + 0.3/1.27 * y_max;
        line([h.gpos(1, iTm) h.gpos(2, iTm)], [y_max y_max], 'Color', 'k', 'LineWidth', 1.5); hold on
        text(mean(h.gpos(:, iTm)), y_max + 0.05/1.27 * y_max ,repmat('*', 1, countLeadingZerosAfterDecimal(pLoadFR_adj(iR,iTm), 'maxNum',3)), 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
    
end
fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegionsBoxPlot_%s.pdf', tag)))

%%



%%
close all
taskMarkers = [500 2600 4700 6800 9400]+500;
respLatency = respLatency(~isnan(respLatency));
Lim = [-1.5e3 2e3];
sm = 100;
win = [-3e3 500];
figure('Position', [904 414 round(285/2) 522]);

subplot(3,1,1) 

plot_condition = true(1, length(correct_trials_all));
plot_condition = correct_trials_all == 1 ;%& load_trials_all == 3; % & respLatency < 1000;
plot_condition_2 = correct_trials_all == 0;% & load_trials_all == 3; % & respLatency < 1000;
plot_condition = load_trials_all == 3; % & respLatency < 1000;
plot_condition_2 = load_trials_all == 1; % & respLatency < 1000;
% plot_condition = respLatency < 1500 & load_trials_all == 3;
% plot_condition_2 = respLatency > 1500 & load_trials_all == 3;
plot_condition = probe_in_out == 1;
plot_condition_2 = probe_in_out == 0;

% plot_condition = load_trials_all == 3;
% plot_condition_2 = load_trials_all == 1;
% plot_condition = probe_in_out == 1;
% plot_condition = load_trials_all == 3;

recog_trials_resp_all(isinf(recog_trials_resp_all)) = nan;
% plot_data = smoothdata(resp_region_all(plot_condition,:,3), 2, 'gaussian',200);
plot_data = smoothdata(recog_trials_resp_all(plot_condition_2,:), 2, 'gaussian',sm);

[bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2)), 'k', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;

% plot_data = smoothdata(c(plot_condition_2,:,3), 2, 'gaussian',200);
plot_data = smoothdata(recog_trials_resp_all(plot_condition,:), 2, 'gaussian',sm);

[bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
xlim([win(1) win(2)])

subplot(3,1,2)
histogram(-respLatency(plot_condition), win(1):100:win(2), 'Normalization', 'probability'); hold on;
histogram(-respLatency(plot_condition_2), win(1):100:win(2), 'Normalization', 'probability')



subplot(3,1,3) 

% plot_condition = true(1, length(correct_trials_all));
% plot_condition = probe_in_out == 1;
% plot_condition = load_trials_all == 3;

plot_data = smoothdata(recog_trials_units_resp_all(plot_condition,:), 2, 'gaussian',sm);

[bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(recog_trials_units_resp_all(plot_condition_2,:), 2, 'gaussian',sm);

[bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2)), 'k', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
xlim([win(1) win(2)])

figure('Position', [893 434 175 196]);
h1 = histogram(-respLatency(plot_condition), win(1):100:win(2), 'Normalization', 'probability'); hold on;
h1.FaceColor = 'r';
h2 = histogram(-respLatency(plot_condition_2), win(1):100:win(2), 'Normalization', 'probability');
h2.FaceColor = 'k';
fig = gcf;
fig.Color = 'w';
ylabel('prop.')
xlabel('time from KP [ms]')
savepdf(gcf, fullfile(exportDirFigs, sprintf('respTime_load_%s.pdf', tag)))

figure('Position',[669 805 285 round(729/5)]);
% plot_condition = true(1, length(correct_trials_all));
plot_condition = load_trials_all == 3 & respLatency > median(respLatency(load_trials_all == 3));
plot_condition_2 = load_trials_all == 3 & respLatency <= median(respLatency(load_trials_all == 3));
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');
whole_trial_all 
plot_data    = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian',100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(whole_trial_all(plot_condition_2,:), 2, 'gaussian',100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2));
plot_data_mu(isnan(whole_trial_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
vl = vline([500 2600 4700 6800 9400]-500); 
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_all(:,1:1e3), 'all', 'omitnan'));
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-1e3 max(taskMarkers) + 2e3])
fig = gcf;
fig.Color = 'w';
ylabel('co-ripples')
xlabel('time from trial start [ms]')
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTask_allChannels_load_%s.pdf', tag)))

figure('Position',[669 805 285 round(729/5)]);
plot_data    = smoothdata(whole_trial_units_all(plot_condition,:), 2, 'gaussian',100) * 1e3;
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_units_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_units_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(whole_trial_units_all(~plot_condition,:), 2, 'gaussian',100) * 1e3;
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_units_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_units_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
vl = vline([500 2600 4700 6800 9400]-500); 
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_units_all(:,1:1e3), 'all', 'omitnan')*1e3);
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
fig = gcf;
fig.Color = 'w';
xlim([-1e3 max(taskMarkers) + 2e3])
ylabel('unit firing')
xlabel('time from trial start [ms]')
savepdf(gcf, fullfile(exportDirFigs, sprintf('Firing_allChannels_load_%s.pdf', tag)))



figure('Position',[669 805 285 round(729/5)]);

% plot_condition = true(1, length(correct_trials_all));
plot_condition = probe_in_out == 1;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
plot_data    = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(whole_trial_all(~plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline(-1e3:length(plot_data)-1e3-1, plot_data_mu, plot_data_sem, 'k', 'nan', 'gap'); hold on;
bl2.Color = [247, 133, 177]./255;
bf.FaceColor = bl2.Color;
bf.FaceAlpha = 0.25;
vl = vline([500 2600 4700 6800 9400]-500); 
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
hl = hline(mean(whole_trial_all(:,1:1e3), 'all', 'omitnan'));
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-1e3 max(taskMarkers) + 2e3])
fig = gcf;
fig.Color = 'w';
ylabel('co-ripples')
xlabel('time from trial start [ms]')
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTask_allChannels_probe_%s.pdf', tag)))



%%
regionColors =  brewermap(12, 'Dark2');

figure('Position',[669 805 285 729]);
plot_condition = load_trials_all == 3 & respLatency > median(respLatency(load_trials_all == 3));
plot_condition_2 = load_trials_all == 3 & respLatency <= median(respLatency(load_trials_all == 3));
% plot_condition = probe_in_out == 1;
% plot_condition = logical(correct_trials_all);
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));
blall = [];
binSz = 1000;
bins = 1:binSz:size(whole_trials_region_all,2);
times = 1:size(whole_trials_region_all,2);
for iR = 1:length(regions)
    subplot(5,1,iR)



    plot_data_1 = whole_trials_region_all(plot_condition_2,:,iR);
%     plot_data_1(plot_data_1 > 0) = 1;
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian',100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',100);
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_region_all(plot_condition,:,iR);
%     plot_data_3(plot_data_3 > 0) = 1;

%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian',100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',100);
    
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
   
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s ripples response p = %.2f\n', regions{iR}, p_perm)
    
    probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
    probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
    probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s ripples maintenance p = %.2f\n', regions{iR}, p_perm)



    for iB = 1:length(bins)-1
        
        
%         data = recog_trials_all_region(plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(plot_condition)), 'o'); hold on;
%         bl1.Color = regionColors(iR,:);
        bl1.MarkerFaceColor = bl1.Color;
%     
%         bf.FaceColor = regionColors(iR,:);
% 
%         data = recog_trials_all_region(~plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(~plot_condition)), 'ko');
%         hold on;
%         bl1.MarkerFaceColor = bl1.Color;

    

    
    end
    
%     xlim([taskMarkers(end-1) taskMarkers(end)+2.0e3])
    xlim([0 taskMarkers(end)+2.0e3])

%     title(regions{iR})

%     legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


%     blall(iR) =bl;
%     xlim(Lim)
%     ylim([0.3 1])
    xlim([-1e3 max(taskMarkers) + 2e3])
    vl = vline([500 2600 4700 6800 9400]-500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
    ylabel('co-ripples')
    xlabel('time from trial start [ms]')
    
    ax = gca;
    ax.FontSize = 10;
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegions_%s.pdf', tag)))
%%
figure('Position',[669 805 round(285/2) 729]);
plot_condition = load_trials_all == 3;
% plot_condition = correct_trials_all == 1;
% plot_condition = probe_in_out == 1;
% plot_condition = logical(correct_trials_all);
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));
blall = [];
binSz = 1000;
bins = 1:binSz:size(whole_trials_region_all,2);
times = 1:size(whole_trials_region_all,2);
for iR = 1:length(regions)
    subplot(5,1,iR)

    baseline = nan(6500,1);
    
    for iS = 1:max(subjID)
        ii = find(subjID == iS);
        baseline(ii) = mean(whole_trials_region_all(ii,1:1e3,iR), 'all', 'omitnan');
        
    end
    
    

    plot_data_1 = resp_region_all(~plot_condition,:,iR); % ./ baseline(~plot_condition);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian',100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',100);
    [bl2, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = resp_region_all(plot_condition,:,iR); %./ baseline(plot_condition);;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian',sm);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',sm);

    [bl1, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
   
    
    hl = hline(mean([plot_data_mu_1(end-1e3:end), plot_data_mu_3(end-1e3:end)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;

    xlim([-3.0e3 500])

    vl = vline(0); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
%     ylabel('co-ripples / trial')
    xlabel('time from KP [ms]')
    
    ax = gca;
    ax.FontSize = 9;
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegionsResp_%s.pdf', tag)))
%%
figure('Position',[669 805 round(285/2) 729]);
plot_condition = load_trials_all == 3 & respLatency > median(respLatency(load_trials_all == 3));
plot_condition_2 = load_trials_all == 3 & respLatency <= median(respLatency(load_trials_all == 3));
% plot_condition = correct_trials_all == 1;
% plot_condition = probe_in_out == 1;
% plot_condition = logical(correct_trials_all);
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));
blall = [];
binSz = 1000;
bins = 1:binSz:size(whole_trials_region_all,2);
times = 1:size(whole_trials_region_all,2);
for iR = 1:length(regions)
    subplot(5,1,iR)

    baseline = nan(6500,1);
    
    for iS = 1:max(subjID)
        ii = find(subjID == iS);
        baseline(ii) = mean(whole_trials_units_sm_region_all(ii,1:1e3,iR), 'all', 'omitnan');
        
    end
    
    

    plot_data_1 = resp_units_region_all(~plot_condition,:,iR) * 1e3;% ./ baseline(~plot_condition);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    [bl2, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = resp_units_region_all(plot_condition,:,iR)*1e3; %./ baseline(plot_condition);
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));

    [bl1, bf] = boundedline(winResp(1):winResp(2), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
   
    
    hl = hline(mean([plot_data_mu_1(end-1e3:end), plot_data_mu_3(end-1e3:end)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;

    xlim([-3e3 500])

    vl = vline(0); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
%     ylabel('units')
    xlabel('time from KP [ms]')
    
    ax = gca;
    ax.FontSize = 9;
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegionsResp_%s.pdf', tag)))

figure('Position',[669 805 285 729]);

for iR = 1:length(regions)
    subplot(5,1,iR)



    plot_data_1 = whole_trials_units_sm_region_all(plot_condition_2,:,iR) * 1e3;
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_units_sm_region_all(plot_condition,:,iR) * 1e3;
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    
    
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 1.0;



    for iB = 1:length(bins)-1
        
        
%         data = recog_trials_all_region(plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(plot_condition)), 'o'); hold on;
%         bl1.Color = regionColors(iR,:);
        bl1.MarkerFaceColor = bl1.Color;
%     
%         bf.FaceColor = regionColors(iR,:);
% 
%         data = recog_trials_all_region(~plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(~plot_condition)), 'ko');
%         hold on;
%         bl1.MarkerFaceColor = bl1.Color;

    

    
    end
    
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s firing p = %.2f\n', regions{iR}, p_perm)

%     title(regions{iR})

%     legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')
%      xlim([taskMarkers(end-1) taskMarkers(end)+2.0e3])
     xlim([1 taskMarkers(end)+2.0e3])


%     blall(iR) =bl;
        xlim([-1e3 max(taskMarkers) + 2e3])

    vl = vline([500 2600 4700 6800 9400]-500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
%     xlim(Lim)
%     ylim([0.3 1])
    ylabel('unit firing [Hz]')
    xlabel('time from trial start [ms]')
    ax = gca;
    ax.FontSize = 10;
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegions_%s.pdf', tag)))

figure('Position', [669 805 285 729]);

for iR = 1:length(regions)
    subplot(5,1,iR)



    plot_data_1 = whole_trials_ERP_region_all(~plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_ERP_region_all(plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    
    vl = vline([500 2600 4700 6800 9400]-500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 1.0;


    
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s erp p = %.2f\n', regions{iR}, p_perm)

    title(regions{iR})
    


    legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


    blall(iR) =bl;
%     xlim(Lim)
%     ylim([0.3 1])
    ylabel('ERP [uV]')
    xlabel('time from trial start [ms]')
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('ERPTaskWithinRegions_%s.pdf', tag)))

%%
figure('Position', [1281 840 772 733]);
% c = 1;
tabCoRip = cell(5);
lmeCoRip = cell(5);
lmePval = nan(5);
for iRa = 1:length(regions)
    for iRb = iRa+1:length(regions)
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)


        plot_data_1 = whole_trials_cross_region_all{iRa,iRb}(~plot_condition,:);
%         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
        plot_data_mu_1 = mean(plot_data_1, 'omitnan');
        plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~plot_condition));
        plot_data_mu_1 = smoothdata(plot_data_mu_1, 'gaussian',100);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian',100);
        plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
        [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
        bf.FaceAlpha = 0.7;
        bl2.LineWidth = 0.5;

        plot_data_3    = whole_trials_cross_region_all{iRa,iRb}(plot_condition,:);

        
%         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
        
        plot_data_mu_3 = mean(plot_data_3, 'omitnan');
        plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(plot_condition));
        plot_data_mu_3 = smoothdata(plot_data_mu_3, 'gaussian',100);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian',100);
        plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

        [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
        bf.FaceAlpha = 0.7;
        bl1.Color = mean([regionColors(iRa,:);regionColors(iRb,:)]) ;
        bf.FaceColor = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        bl1.MarkerFaceColor = bl1.Color;
        bf.FaceAlpha = 0.3;
        bl1.LineWidth = 0.5;
        
%         if iRa ~= iRb
%             Arate    = smoothdata(whole_trials_region_all(plot_condition,:,iRa), 2, 'gaussian',500);
%             Arate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
%             Brate = smoothdata(whole_trials_region_all(plot_condition,:,iRb), 2, 'gaussian',500);
%             Brate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
%             plot_data_3_AB    = mean(Arate, 'omitnan').*mean(Brate, 'omitnan');
%             pl = plot(-1e3:length(plot_data_3)-1e3-1, plot_data_3_AB, 'k--'); hold on;
%             pl.LineWidth = 0.5;
%         end

%         pl = hline(plot_data_3_AB, '--'); hold on;
%         pl.Color = bl1.Color;


        vl = vline([500 2600 4700 6800 9400]-500); 

        for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
        
        hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));   
        hl.Color = 'k';
        hl.LineStyle = '-';
        hl.LineWidth = 0.5;
        
        
        
        
        
%         pVal = nan(1, length(bins)-1);
%         dVal = nan(2, length(bins)-1);
%         plotTimesStats = nan(1, length(bins)-1);
%         for iB = 1:length(bins)-1
%         
%         
%             data1 = trapz(whole_trials_cross_region_all{iRa,iRb}(~plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
%             data1(isnan(data1)) = [];
%             dVal(1,iB) = mean(data1);
%             data3 = trapz(whole_trials_cross_region_all{iRa,iRb}(plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
%             data3(isnan(data3)) = [];
%             dVal(2,iB) = mean(data3);
% 
%             plotTimesStats(iB) = bins(iB) + (binSz/2);
%             
%             if ~isempty(data1) && ~isempty(data3)
% %                 [~,~,p_perm] = statcond({data1',data3'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%                 [pRank,H] = ranksum(data1, data3);
%                 [H,pT,CI] = ttest2(data1,data3);
%                 pVal(iB) = pRank;
%             end
% 
% 
%            
% 
%     
% 
%     
%         end
%         adj_p = nan(size(pVal));
%         yVal = repmat(max(plot_data_mu_3), [1 length(pVal)]);
%         [h, crit_p, adj_ci_cvrg, adj_p(~isnan(pVal))]=fdr_bh(pVal(~isnan(pVal)),0.05,'pdep','no');
%         
%         bnds = mask2bounds(adj_p < 0.05);
%         
%         for iBnd = 1:size(bnds,1)
%             pl = plot(plotTimesStats(bnds(iBnd,1):bnds(iBnd,1)), yVal(bnds(iBnd,1):bnds(iBnd,1)), 'k*'); hold on;
%             pl.LineWidth = 1.5;
%         end
        

         coRresp = mean(whole_trials_cross_region_all{iRa,iRb}(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2, 'omitnan');
         coRmaintenance = mean(whole_trials_cross_region_all{iRa,iRb}(:,taskMarkers(end-1):taskMarkers(end)), 2, 'omitnan');
         coRbaseline = mean(whole_trials_cross_region_all{iRa,iRb}(:,1:taskMarkers(1)), 2, 'omitnan');
         coRresp(cTrial:end) = [];
         coRmaintenance(cTrial:end) = [];
         coRbaseline(cTrial:end) = [];
        
         tabCoRip{iRa,iRb} = table(subjID(~isnan(coRresp)), plot_condition(~isnan(coRresp)), coRresp(~isnan(coRresp)),coRmaintenance(~isnan(coRresp)),coRbaseline(~isnan(coRresp)), ...
                 'VariableNames', {'subject',               'condition',                     'rates_resp',            'rates_maintenance',            'rates_baseline'});
%          lmeCoRip{iRa,iRb} = fitlme(tabCoRip{iRa,iRb}, 'rates_resp ~ condition + (1|subject)');   
         lmeCoRip{iRa,iRb} = fitlme(tabCoRip{iRa,iRb}, 'rates_maintenance ~ condition + (1|subject)');   
         lmePval(iRa,iRb) = lmeCoRip{iRa,iRb}.Coefficients(2,6);


        
        if iRa == 1; ylabel('co-ripples'); 
        end
        
%         probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); probeResp_1(isnan(probeResp_1)) = [];
%         probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); probeResp_3(isnan(probeResp_3)) = [];
%         [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%         [pRank,H] = ttest2(coRmaintenance(plot_condition), coRmaintenance(~plot_condition));
        pct = [mean(coRmaintenance(plot_condition), 'omitnan') - mean(coRmaintenance(~plot_condition), 'omitnan')]/ mean(coRmaintenance(~plot_condition), 'omitnan');
        fprintf('%s <--> %s ripples maintenance %.2f\n', regions{iRa}, regions{iRb}, pct)
% 
%         probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_1(isnan(probeResp_1)) = [];
%         probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_3(isnan(probeResp_3)) = [];
%         [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%         [pRank,H] = ranksum(probeResp_3, probeResp_1);

%         fprintf('%s <--> %s ripples maintenance p = %.2f\n', regions{iRa}, regions{iRb},  p)

        plotMarkers = [500 2600 4700 6800 9400]-500;
%         blall(iR) =bl;
%         xlim([0 13000])
        xlim([plotMarkers(end-1)-200 plotMarkers(end)+2.0e3])
%         xlim([-1e3 max(taskMarkers) + 2e3])

        if iRb == length(regions); xlabel('time [ms]'); end
        
%         c = c+1;
    end
end
fig= gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskAcrossRegions_%s.pdf', tag)))

adj_p = nan(size(lmePval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(lmePval))]=fdr_bh(lmePval(~isnan(lmePval)),0.05,'pdep','yes');

figure('Position', [1738 1327 246 182]); 
imagesc(adj_p', [min(adj_p(:))-1e-2 1]); hold on;
colorbar;
hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
cmap = flipud(slanCM('magma'));
cmap = [1 1 1; cmap];
colormap(cmap)
[yy xx] = find(adj_p' < 0.05);
plot(xx, yy, 'r*'); hold on;



ax = gca;
ax.YTick = 1:size(adj_p,1);
ax.XTick = 1:size(adj_p,1);
ax.XTickLabel = regions;
ax.YTickLabel = regions;
fig = gcf;
fig.Color = 'w';
% savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskAcrossRegions_stats_probe_%s.pdf', tag)))
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskAcrossRegions_stats_maintenance_%s.pdf', tag)))
%%

figure('Position', [1209 1255 254 294]);
br = bar([5 12; 0 1; 2 7]'); hold on;
br(1).FaceColor = [230, 159, 0]/255;
br(2).FaceColor = [0, 158, 115]/255;
br(3).FaceColor = [213, 94, 0]/255;
ylim([0 15]);

ax = gca;
ax.XTickLabel = {'M','P'};

legend('co-R', 'co-LG', 'co-vHG', 'location', 'northwest')
box off
fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('AcrossRegionsBar_stats_%s.pdf', tag)))

%%
condNames = {'All', 'coRip', 'noRip'};
respData = cell(2, length(condNames), 10);
for cond = 1
    % Find all open figure objects
    openFigures = findobj('Type', 'figure');
    % Get the number of open figures
    numOpenFigures = numel(openFigures);
    
    tabCoFire = cell(5);

    figure('Position', [1281 291 1280 1282]);
    figure('Position', [255 142 617 1020]);
    coRegionColors = nan(10,3);
    for iRa = 1:length(regions)
        for iRb = iRa+1:length(regions)
        figure(numOpenFigures+1)
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)            
%             baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:1e3), 'all', 'omitnan')/ ...
%                              mean(whole_trials_cross_region_all{iRa,iRb}(:,1:1e3)>0, 'all', 'omitnan');
                         
            baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:1e3), 'all', 'omitnan');


%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(~plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_1 = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(~plot_condition,:))/baseline, 2, 'gaussian',500);
    %         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
            plot_data_mu_1 = mean(plot_data_1, 'omitnan');
            plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
            plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
            plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
            [bl2, bf] = boundedline(-1e3:length(plot_data_1)-1e3-1, plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
            bf.FaceAlpha = 0.7;
            
%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_3    = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(plot_condition,:))/baseline, 2, 'gaussian',500);
    %         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
            plot_data_mu_3 = mean(plot_data_3, 'omitnan');
            plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
            plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;
            plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

            [bl1, bf] = boundedline(-1e3:length(plot_data_3)-1e3-1, plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
            bf.FaceAlpha = 0.7;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
            bl1.Color = coRegionColors(c,:)  ;
            bf.FaceColor = coRegionColors(c,:) ;
            bl1.MarkerFaceColor = bl1.Color;
            bf.FaceAlpha = 0.3;


            vl = vline([500 2600 4700 6800 9400]-500); 

            for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end

%             hl = hline(baseline);
            hl = hline(1);
            hl.Color = 'k';
            hl.LineStyle = '-';
            hl.LineWidth = 1.0;


            title([regions{iRa}, ' <--> ', regions{iRb}])

        
            subjID_1 = subjID(~plot_condition); subjID_3 = subjID(plot_condition);
            
            
            coRresp = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2, 'omitnan');
            coRmaintenance = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,taskMarkers(end-1):taskMarkers(end)), 2, 'omitnan');
            coRbaseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:taskMarkers(1)), 2, 'omitnan');
            coRresp(cTrial:end) = [];
            coRmaintenance(cTrial:end) = [];
            coRbaseline(cTrial:end) = [];
            
            tabCoFire{iRa,iRb} = table(subjID(~isnan(coRresp)), plot_condition(~isnan(coRresp)), coRresp(~isnan(coRresp)),coRmaintenance(~isnan(coRresp)),coRbaseline(~isnan(coRresp)), ...
                    'VariableNames', {'subject',               'condition',                     'rates_resp',            'rates_maintenance',            'rates_baseline'});
         
        
            
            
           


            blall(iR) =bl;
            xlim([0 13000])
        %     ylim([0.3 1])
            ylabel('co-firing')
            xlabel('time from trial start [ms]')
            
            tID = tabCoFire{iRa,iRb}.subject;  
            tResp = tabCoFire{iRa,iRb}.rates_resp;  
            tCond = tabCoFire{iRa,iRb}.condition;  

            tRespID = [arrayfun(@(X) mean(tResp(tID == X & tCond == 0)), 1:43)',arrayfun(@(X) mean(tResp(tID == X & tCond == 1)), 1:43)'] ;

            figure(numOpenFigures+2); 
            subplot(5,5,c)            
            pl = parallelcoords(tRespID, 'Labels', {'load1', 'load3'}, 'Color', bl1.Color); hold on;%,'quantile',.05);
            xlim([0.25 2.75])
            ylim([min(tRespID(:)) max(tRespID(:)) ])
            
            xx = tRespID(:,1);
            er = errorbar(0.5, mean(xx, 'omitnan'), std(xx, 'omitnan')/sqrt(sum(~isnan(xx))), 'o'); hold on;
            er.Color = bl2.Color;
            xx = tRespID(:,2);
            er = errorbar(2.5, mean(xx, 'omitnan'), std(xx, 'omitnan')/sqrt(sum(~isnan(xx))), 'o'); hold on;
            er.Color = bl1.Color;
            
            [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,2));
            
            if P < 0.05; plot(1.5, max( max(tRespID(:))), 'k*'); hold on; end
            
            
            ax = gca;
            ax.XTick = [0.5 2.5];

        end
    end
    
    
    
 
    fig= gcf;
    fig.Color = 'w';
    savepdf(gcf, fullfile(exportDirFigs, sprintf('coFireTaskAcrossRegions_%s_%s.pdf', condNames{cond}, tag)))
end

%%

figure;
numGroups = size(respData, 3);


groupMeanSubj = [];

for cond = 2:3
    for ii = 1:numGroups
        subplot(2,2,cond-1)
        groupMean = mean(respData{1,cond,ii}(:,2), 'omitnan');
        groupSEM = std(respData{1,cond,ii}(:,2), 'omitnan')/sqrt(sum(~isnan(respData{1,cond,ii}(:,2))));

        % Plot the SEM
        ln = line([ii-0.11 ii-0.11], [groupMean - groupSEM, groupMean + groupSEM], ...
            'Color', 'r', 'LineWidth', 1.5); hold on;
        ln.Color = coRegionColors(ii,:)*0.5 ;
        
         % Plot the mean
        pl = plot(ii-0.11, groupMean, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
        pl.Color= 'k';
        pl.MarkerFaceColor= coRegionColors(ii,:)*0.5;
        
        groupMean = mean(respData{2,cond,ii}(:,2), 'omitnan');
        groupSEM = std(respData{2,cond,ii}(:,2), 'omitnan')/sqrt(sum(~isnan(respData{2,cond,ii}(:,2))));
         % Plot the SEM
        ln = line([ii+0.11 ii+0.11], [groupMean - groupSEM, groupMean + groupSEM], ...
            'Color', 'r', 'LineWidth', 1.5); hold on;
        ln.Color = coRegionColors(ii,:) ;
        
%         % Plot the mean
%         pl = plot(ii+0.11, groupMean, 'rs', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%         pl.Color= 'k';
%         pl.MarkerFaceColor= coRegionColors(ii,:);
%         for subj = 1:length(subj_list_full)
%             subplot(2,2,cond-1+2)
%             iSub = respData{1,cond,ii}(:,1) == subj;
%             groupMeanSubj(subj, 1) = mean(respData{1,cond,ii}(iSub,2), 'omitnan');
%             groupSEM = std(respData{1,cond,ii}(iSub,2), 'omitnan')/sqrt(sum(~isnan(respData{1,cond,ii}(iSub,2))));
% 
%             % Plot the SEM
%             ln = line([ii-0.11 ii-0.11], [groupMeanSubj(subj, 1) - groupSEM, groupMeanSubj(subj, 1) + groupSEM], ...
%                 'Color', 'r', 'LineWidth', 1.5); hold on;
%             ln.Color = coRegionColors(ii,:)*0.5 ;
% 
%              % Plot the mean
%             pl = plot(ii-0.11, groupMeanSubj(subj, 1), 'rs', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%             pl.Color= 'k';
%             pl.MarkerFaceColor= coRegionColors(ii,:)*0.5;
%             
%             iSub = respData{2,cond,ii}(:,1) == subj;
%             groupMeanSubj(subj, 2) = mean(respData{2,cond,ii}(iSub,2), 'omitnan');
%             groupSEM = std(respData{2,cond,ii}(iSub,2), 'omitnan')/sqrt(sum(~isnan(respData{2,cond,ii}(iSub,2))));
%              % Plot the SEM
%             ln = line([ii+0.11 ii+0.11], [groupMeanSubj(subj, 2) - groupSEM, groupMeanSubj(subj, 2) + groupSEM], ...
%                 'Color', 'r', 'LineWidth', 1.5); hold on;
%             ln.Color = coRegionColors(ii,:) ;
% 
%             % Plot the mean
%             pl = plot(ii+0.11, groupMean(subj, 2), 'rs', 'MarkerSize', 8, 'LineWidth', 1.5); hold on;
%             pl.Color= 'k';
%             pl.MarkerFaceColor= coRegionColors(ii,:);
%             
%         end
%         
%         
%        
    end
    
    box off
    xlim([0 numGroups+1])
%     ylim([-0.2 1.4])
    ax = gca;
    ax.XAxisLocation = 'origin';
end

figure('Position', [894 563 720 349]);
subplot(3,1,[1 2])
histogram(respLatency, 0:50:max(respLatency))
whiskerLimits = [quantile(respLatency, 0.01), quantile(respLatency, 0.95)];
xlim([0 whiskerLimits(2)]);
subplot(3,1,3)
h = boxplot(respLatency, 'Orientation', 'horizontal', 'Symbol', '');
xlabel('Response Latency')
xlim([0 whiskerLimits(2)]);

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'respLatency.pdf'))

% ylim([0.9 1.1])
% %savepdf(gcf, 'MemoryResponses_microOnly.pdf')
% figure;
% m = recog_trials_all(plot_condition,(abs(win(1)) + 500):(abs(win(1)) + 2000));
% muNovel = mean(m, 2, 'omitnan');
% s = std(muNovel(:), 'omitnan') / sqrt(length(muNovel(:)));
% errorbar(1,mean(muNovel, 'omitnan'),s, 'ro'); hold on;
% 
% m = recog_trials_all(~plot_condition,(abs(win(1)) + 500):(abs(win(1)) + 2000));
% muRecog = mean(m, 2, 'omitnan');
% s = std(muRecog(:),  'omitnan') / sqrt(length(muRecog(:)));
% errorbar(1,mean(muRecog,  'omitnan'),s, 'bo'); hold on;
%%
figure('Position', [26 1070 825 190]);
subplot(1,4,1)
h = boxplot(densityAll, 'Orientation', 'vertical', 'Symbol', '');
ax = gca;
ax.XTick = 1;
ax.XTickLabel = 'ripple density';
ylim([quantile(densityAll, 0.01) quantile(densityAll, 0.99)]);
whiskers = findobj(gca, 'Tag', 'Outliers');
set(whiskers, 'LineWidth', 4);
box off

subplot(1,4,2)
h = boxplot(durAll, 'Orientation', 'vertical', 'Symbol', '');
ax = gca;
ax.XTick = 1;
ax.XTickLabel = 'ripple duration [ms]';
ylim([quantile(durAll, 0.01) quantile(durAll, 0.99)]);
box off

subplot(1,4,3)
h = boxplot(ampAll, 'Orientation', 'vertical', 'Symbol', '');
ax = gca;
ax.XTick = 1;
ax.XTickLabel = 'ripple amplitude [uV]';
ylim([quantile(ampAll, 0.01) quantile(ampAll, 0.99)]);
box off

subplot(1,4,4)
h = boxplot(freqAll, 'Orientation', 'vertical', 'Symbol', '');
ax = gca;
ax.XTick = 1;
ax.XTickLabel = 'ripple amplitude [uV]';
ylim([quantile(freqAll, 0.01) quantile(freqAll, 0.99)]);
box off



fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'rippleMetrics.pdf'))



























