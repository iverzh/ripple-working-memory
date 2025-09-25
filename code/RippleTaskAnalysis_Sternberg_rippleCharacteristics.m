

close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%
dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess';
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;


matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));

LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '*task*'));
taskfiles = {taskfiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../out', '*micro*'));
micrFiles = {micrFiles.name}';

recordingState = 'wake';
location = 'NC';

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
winResp = [-1.5 1] * 1e3;
winLength = 12905; %ms

%prealocate variables
recog_trials_all = nan(6.5e3, length(win(1):win(2)));
recog_trials_units_all = nan(6.5e3, length(win(1):win(2)));

recog_trials_resp_all = nan(6.5e3, length(winResp(1):winResp(2)));
% recog_trials_all_region = [];

whole_trial_all = nan(6.5e3, winLength);
whole_trials_region_all = nan(6.5e3, winLength, length(regions));
whole_trials_units_region_all  = nan(6.5e3, winLength, length(regions));
whole_trials_cross_region_all = cell(length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        whole_trials_cross_region_all{iRa, iRb} = nan(6.5e3, winLength);       
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
densityAll = [];

uChanAll = [];
tic
cTrial = 1;
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
    trials = load(fullfile(dataDirectory, taskfiles{f}));

    
    modifier = 'IISremv_singleWire_v3_z25';
%     modifier = 'IISremv_singleWire_v3';
%     modifier = 'IISremv_v3';
    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
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
%             elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); continue; end
            
            iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
            iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
            
            
        
            for ii = 1:length(iE)
%                 iS(ii) = find(times == iS(ii));
%                 iE(ii) = find(times == iE(ii));
                if any([iS(ii) iE(ii)] <= 0); continue; end
                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
        
    end
    
    rippMask(:, end:end+3e3) = nan; %pad the ending


    rippAll = sum(rippMask);
    
    
    units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
    if isempty(units); continue; end
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
    uChan = cell2mat(units(:,1));
    uChanAll = [uChanAll uChan'];
    uLocations = units(U,end);
    uLocations = strrep(uLocations, 'ventral_medial_prefrontal_cortex', 'OFC');
    uLocations = strrep(uLocations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    uLocations = strrep(uLocations, 'pre_supplementary_motor_area', 'SMA');
    uLocations = strrep(uLocations, 'amygdala', 'AMY');
    uLocations = strrep(uLocations, 'hippocampus', 'HIP');
        
        
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
        FR = sum(N) / rippleStats.recordingLength * 1000;
    %     if FR > 0
        spikeMask(iU,:) = N;
        spikeMask_sm(iU,:) = smoothdata(N, 'gaussian',500);
    %     end
    end
    if nNeuron > 1
        unitAll = smoothdata(sum(spikeMask, 'omitnan'), 'gaussian',500);
    else
        unitAll = smoothdata(spikeMask, 'gaussian',500);
    end
    unitAll = zscore(unitAll); % - mean(unitAll);

    unit_trials_region = nan(length(trials.start_time), 12905, length(regions));
    dat_region = cell(1,length(regions));
    u_dat_region_sm = cell(1,length(regions));
    u_dat_region = cell(1,length(regions));
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


        dat = rippAll(respTime+winResp(1):respTime+winResp(2));
        recog_trials_resp_all(cTrial, :) = dat;
        
        if sum(dat) == 0; continue; end

        
        trialLoad = trials.loads(iT);
        image1Time =  round(trials.timestamps_Encoding1(iT)*1e3);
        image2Time =  round(trials.timestamps_Encoding2(iT)*1e3);
        image3Time =  round(trials.timestamps_Encoding3(iT)*1e3);
        maintenanceTime =  round(trials.timestamps_Maintenance(iT)*1e3);
        probeTime =  round(trials.timestamps_Probe(iT)*1e3);

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
            im1seg = spikeMask(:, image1Time-1000:image1Time+2e3);
            im2seg = spikeMask(:, image2Time-0  :image2Time+2e3);
            im3seg = spikeMask(:, image3Time-0  :image3Time+2e3);
            maintenanceSeg = spikeMask(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = spikeMask(:, probeTime-1250:probeTime+3000);
            
            
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
            im1seg = spikeMask(:,image1Time-1000:image1Time+2e3);
            im2seg = nan(nUnit, 2e3+0+1);
            im3seg = nan(nUnit, 2e3+0+1);
            maintenanceSeg = spikeMask(:, maintenanceTime-0:maintenanceTime+1250);
            probeSeg = spikeMask(:, probeTime-1250:probeTime+3000);
            
            
            unitTrialFullChan = [im1seg, nan(nUnit,100), ...
                         im2seg, nan(nUnit,100), ...
                         im3seg, nan(nUnit,100), ...
                         maintenanceSeg, nan(nUnit,100), ...
                         probeSeg];
        end
        
        whole_trial_all(cTrial,:) = trialFull ;
        





        %within region co-rippling and unit firing
        for iR = 1:length(regions)
            chRegion = contains(rippleStats.chanLabels, regions(iR));
            uRegion = contains(uLocations, regions(iR));

            if isempty(dat_region{iR}) && sum(chRegion) > 0
                 dat_region{iR} = sum(rippMask(chRegion,:),1);
            elseif sum(chRegion) == 0
                 dat_region{iR} = nan(1, length(rippMask));
            end

            if isempty(u_dat_region_sm{iR}) && sum(uRegion) > 0
                 u_dat_region_sm{iR} = sum(spikeMask_sm(uRegion,:),1);
                 u_dat_region{iR} = sum(spikeMask(uRegion,:),1);
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
            whole_trials_units_region_all(cTrial, :, iR) = uTrialFullsm;


            
             
%             if length(dat) > 16e3; continue; end
%             trialData(1:length(dat)) = dat;
%             trials_region_load1(iT,:,iR) = trialData;
        end
        
         %across region co-rippling
        for iRa = 1:length(regions)
            for iRb = iRa:length(regions)

                datArip =  whole_trials_region_all(cTrial,:,iRa);
                datBrip =  whole_trials_region_all(cTrial,:,iRb);
                
                if iRa == iRb
                    datAB = datArip;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
                    continue
                else

                    datAB = datArip + datBrip;
                    datAB(datArip <= 0 | datBrip <= 0) = 0;
                    whole_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;

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
                
                datArip =  ripTrialFullChan(contains(locations, regions(iRa)),:);
                datBrip =  ripTrialFullChan(contains(locations, regions(iRb)),:);
%                 
                nanMaskA = isnan(sum(datArip,1)); bndA = mask2bounds(nanMaskA); nanMaskA = bounds2mask(bndA,length(nanMaskA), 50)'; 
                nanMaskB = isnan(sum(datBrip,1)); bndB = mask2bounds(nanMaskB); nanMaskB = bounds2mask(bndB,length(nanMaskB), 50)';                                          
                datArip(:,nanMaskA) = 0; datArip = logical(datArip);
                datBrip(:,nanMaskB) = 0; datBrip = logical(datBrip);
                
               

                datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
                datA = sum(datA,1) > 0;
                datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                datB = sum(datB,1) > 0;

                coFireWindow = round(25/2);
                datA = movsum(datA, [coFireWindow-1 0]);
                datB = movsum(datB, [coFireWindow-1 0]);
                datAB = double((datA > 0) & (datB > 0));
%                 datAB = datAB / sum(~nanMaskA & ~nanMaskB);
                datAB(nanMaskA | nanMaskB) = nan;
                
                whole_trials_unit_cross_region_all{iRa, iRb, 1}(cTrial, :) = datAB;
                                                                 
%                 datA =  unit_trials_region(iT,:,iRa) > 0;
%                 datA(~datArip) = 0;
%                 datB =  unit_trials_region(iT,:,iRb) > 0;    
%                 datB(~datBrip) = 0;

                % co fire during co-ripples
                datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
                datA(: ,sum(datArip) < 1) = 0; datA = sum(datA, 1) > 0;
                datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                datB(:, sum(datBrip) < 1) = 0; datB = sum(datB, 1) > 0;
                
                datA = movsum(datA, [coFireWindow-1 0]);
                datB = movsum(datB, [coFireWindow-1 0]);
                datAB = double((datA > 0) & (datB > 0));
%                 datAB = datAB / sum(datArip & datBrip);

                datAB(nanMaskA | nanMaskB) = nan;

                
                whole_trials_unit_cross_region_all{iRa, iRb, 2}(cTrial, :) = datAB;
                                                                 
%                 datA =  unit_trials_region(iT,:,iRa) > 0;
%                 datA(datArip) = 0;
%                 datB =  unit_trials_region(iT,:,iRb) > 0;    
%                 datB(datBrip) = 0;

                % co fire outside of any rippling
                datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
                datA(: , sum(datArip) > 0) = 0; datA = sum(datA, 1) > 0;
                datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                datB(: , sum(datBrip) > 0) = 0; datB = sum(datB, 1) > 0;

                datA = movsum(datA, [coFireWindow-1 0]);
                datB = movsum(datB, [coFireWindow-1 0]);
                datAB = double((datA > 0) & (datB > 0));
%                 datAB = datAB / sum(~datArip & ~datBrip & ~nanMaskA & ~nanMaskB);
                datAB(nanMaskA | nanMaskB) = nan;

                
                whole_trials_unit_cross_region_all{iRa, iRb, 3}(cTrial, :) = datAB;

                
                                                                 
               
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
    
    
    fprintf('number of total trials: %i, number of subject trials: %i \n\n', cTrial - 1, length(trials.start_time))
 

%     save(fullfile(exportDir, 'ripDetections', sprintf('%s_rippleStats.mat', subject)), 'rippleStats' ,'-v7.3')

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
Lim = [-9e3 2e3]
figure('Position', [904 414 700 522]);

plot_condition = true(1, length(correct_trials_all));
% plot_condition = probe_in_out == 1;
% plot_condition = load_trials_all == 3;

plot_data = smoothdata(recog_trials_resp_all(plot_condition,:), 2, 'gaussian',200);

[bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(recog_trials_resp_all(plot_condition,:), 2, 'gaussian',200);

[bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)

figure('Position',[67 489 2490 447]);

% plot_condition = true(1, length(correct_trials_all));
plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
plot_data    = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(1:length(plot_data), plot_data_mu, plot_data_sem, 'r'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(whole_trial_all(~plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline(1:length(plot_data), plot_data_mu, plot_data_sem, 'b'); hold on;
bf.FaceAlpha = 0.7;
vl = vline([500 2600 4700 6800 9400]+500); 
for iV = 1:length(vl); vl(iV).LineWidth = 2; end


figure('Position',[67 489 2490 447]);

% plot_condition = true(1, length(correct_trials_all));
plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
plot_data    = smoothdata(whole_trial_all(plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline(1:length(plot_data), plot_data_mu, plot_data_sem, 'r'); hold on;
bf.FaceAlpha = 0.7;

plot_data = smoothdata(whole_trial_all(~plot_condition,:), 2, 'gaussian',200);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(~plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline(1:length(plot_data), plot_data_mu, plot_data_sem, 'b'); hold on;
bf.FaceAlpha = 0.7;
vl = vline([500 2600 4700 6800 9400]+500); 
for iV = 1:length(vl); vl(iV).LineWidth = 2; end

figure('Position', [904 414 700 522]);
subplot(2,3,1)
plot_condition = logical(correct_trials_all);
plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',200);
[bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',200);

[bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
legend([bl1, bl2], 'correct', 'incorrect', 'location', 'best')
xlim(Lim)
ylim([0.7 1.9])
ylabel('co-ripples')
xlabel('time from trial start [ms]')


subplot(2,3,4)  
% figure;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
ylim([-0.3 0.3])
xlim(Lim)
ylabel('unit firing [z-score]')
xlabel('time from trial start [ms]')
fig = gcf;
fig.Color = 'w';

subplot(2,3,2)
plot_condition = load_trials_all == 3;
plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',200);
[bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',200);

[bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
xlim(Lim)
ylim([0.7 1.9])
ylabel('co-ripples')
xlabel('time from trial start [ms]')

subplot(2,3,5)  
% figure;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
ylim([-0.3 0.3])
xlim(Lim)
ylabel('unit firing [z-score]')
xlabel('time from trial start [ms]')
fig = gcf;
fig.Color = 'w';

subplot(2,3,3)
plot_condition = probe_in_out == 1;
plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',300);
[bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                                std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',300);

[bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
                              std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')
xlim(Lim)
ylim([0.7 1.9])
% ylim([0 1.5])
ylabel('co-ripples')
xlabel('time from trial start [ms]')

subplot(2,3,6)  
% figure;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
bf.FaceAlpha = 0.7;
[bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
    std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
bf.FaceAlpha = 0.7;
vline(0)
ylim([-0.3 0.3])
xlim(Lim)
ylabel('unit firing [z-score]')
xlabel('time from trial start [ms]')
fig = gcf;
fig.Color = 'w';
%savepdf(gcf, fullfile(exportDirFigs, 'coRipTask.pdf'))

regionColors =  brewermap(12, 'Dark2');

figure('Position', [1 291 1280 1282]);
plot_condition = load_trials_all == 3;
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
    subplot(3,2,iR)



    plot_data_1 = whole_trials_region_all(~plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian',500);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',500);
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline(1:length(plot_data_1), plot_data_mu_1, plot_data_sem, 'k'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_region_all(plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian',500);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',500);
    
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline(1:length(plot_data_3), plot_data_mu_3, plot_data_sem, 'r'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    taskMarkers = [500 2600 4700 6800 9400]+500;
    vl = vline(taskMarkers); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 1.0;
    
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

    title(regions{iR})

    legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


    blall(iR) =bl;
%     xlim(Lim)
%     ylim([0.3 1])
    ylabel('co-ripples')
    xlabel('time from trial start [ms]')
end
fig = gcf; fig.Color = 'w';
%savepdf(gcf, fullfile(exportDirFigs, 'coRipTaskWithinRegions.pdf'))

figure('Position', [1 291 1280 1282]);

for iR = 1:length(regions)
    subplot(3,2,iR)



    plot_data_1 = whole_trials_units_region_all(~plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline(1:length(plot_data_1), plot_data_mu_1, plot_data_sem, 'k'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_units_region_all(plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline(1:length(plot_data_3), plot_data_mu_3, plot_data_sem, 'r'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    
    vl = vline([500 2600 4700 6800 9400]+500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
    
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

    title(regions{iR})

    legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


    blall(iR) =bl;
%     xlim(Lim)
%     ylim([0.3 1])
    ylabel('unit firing')
    xlabel('time from trial start [ms]')
end
fig = gcf; fig.Color = 'w';
%savepdf(gcf, fullfile(exportDirFigs, 'FiringTaskRegions.pdf'))


figure('Position', [1281 291 1280 1282]);
% c = 1;
tabCoRip = cell(5);
lmeCoRip = cell(5);
lmePval = nan(5);
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)


        plot_data_1 = whole_trials_cross_region_all{iRa,iRb}(~plot_condition,:);
%         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
        plot_data_mu_1 = mean(plot_data_1, 'omitnan');
        plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~plot_condition));
        plot_data_mu_1 = smoothdata(plot_data_mu_1, 'gaussian',500);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian',500);
        plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
        [bl2, bf] = boundedline(1:length(plot_data_1), plot_data_mu_1, plot_data_sem, 'k'); hold on;
        bf.FaceAlpha = 0.7;

        plot_data_3    = whole_trials_cross_region_all{iRa,iRb}(plot_condition,:);

        Arate    = smoothdata(whole_trials_region_all(plot_condition,:,iRa), 2, 'gaussian',500);
        Arate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
        Arate = mean(Arate, 'all' , 'omitnan');
        Brate = smoothdata(whole_trials_region_all(plot_condition,:,iRb), 2, 'gaussian',500);
        Brate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
        Brate = mean(Brate, 'all' ,'omitnan');
        plot_data_3_AB    = Arate.*Brate;
%         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
        
        plot_data_mu_3 = mean(plot_data_3, 'omitnan');
        plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(plot_condition));
        plot_data_mu_3 = smoothdata(plot_data_mu_3, 'gaussian',500);
        plot_data_sem = smoothdata(plot_data_sem, 'gaussian',500);
        plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;
        plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

        [bl1, bf] = boundedline(1:length(plot_data_3), plot_data_mu_3, plot_data_sem, 'r'); hold on;
        bf.FaceAlpha = 0.7;
        bl1.Color = mean([regionColors(iRa,:);regionColors(iRb,:)]) ;
        bf.FaceColor = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        bl1.MarkerFaceColor = bl1.Color;
        bf.FaceAlpha = 0.3;
%         pl = hline(plot_data_3_AB, '--'); hold on;
%         pl.Color = bl1.Color;


        vl = vline([500 2600 4700 6800 9400]+500); 

        for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end
        
        hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));   
        hl.Color = 'k';
        hl.LineStyle = '-';
        hl.LineWidth = 1.0;
        
        
        
        
        
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


        title([regions{iRa}, ' <--> ', regions{iRb}])

        legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
        
%         probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); probeResp_1(isnan(probeResp_1)) = [];
%         probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); probeResp_3(isnan(probeResp_3)) = [];
%         [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
        [pRank,H] = ttest2(coRmaintenance(plot_condition), coRmaintenance(~plot_condition));
% 
        fprintf('%s <--> %s ripples response p = %.2f\n', regions{iRa}, regions{iRb}, H)
% 
%         probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_1(isnan(probeResp_1)) = [];
%         probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_3(isnan(probeResp_3)) = [];
%         [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%         [pRank,H] = ranksum(probeResp_3, probeResp_1);

%         fprintf('%s <--> %s ripples maintenance p = %.2f\n', regions{iRa}, regions{iRb},  p)


        blall(iR) =bl;
        xlim([0 13000])
    %     ylim([0.3 1])
        ylabel('co-ripples')
        xlabel('time from trial start [ms]')
        
%         c = c+1;
    end
end
fig= gcf;
fig.Color = 'w';
%savepdf(gcf, fullfile(exportDirFigs, 'coRipTaskAcrossRegions.pdf'))

adj_p = nan(size(lmePval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(lmePval))]=fdr_bh(lmePval(~isnan(lmePval)),0.05,'pdep','yes');

figure('Position', [1210 1274 420 299]); 
imagesc(adj_p, [min(adj_p(:))-1e-2 1]); hold on;
colorbar;
hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
cmap = flipud(slanCM('viridis'));
cmap = [1 1 1; cmap];
colormap(cmap)
[yy xx] = find(adj_p < 0.05);
plot(xx, yy, 'r*'); hold on;



ax = gca;
ax.YTick = 1:size(adj_p,1);
ax.XTick = 1:size(adj_p,1);
ax.XTickLabel = regions;
ax.YTickLabel = regions;
fig = gcf;
fig.Color = 'w';


condNames = {'All', 'coRip', 'noRip'};
respData = cell(2, length(condNames), 10);
respTab = cell(5);
for cond = 1:3
    figure('Position', [1281 291 1280 1282]);
    coRegionColors = nan(10,3);
    for iRa = 1:length(regions)
        for iRb = iRa+1:length(regions)
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)            
%             baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:1e3), 'all', 'omitnan')/ ...
%                              mean(whole_trials_cross_region_all{iRa,iRb}(:,1:1e3)>0, 'all', 'omitnan');
                         
            baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:1e3), 'all', 'omitnan');


%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(~plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_1 = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(~plot_condition,:)), 2, 'gaussian',500);
    %         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
            plot_data_mu_1 = mean(plot_data_1, 'omitnan');
            plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
            plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
            plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
            [bl2, bf] = boundedline(1:length(plot_data_1), plot_data_mu_1, plot_data_sem, 'k'); hold on;
            bf.FaceAlpha = 0.7;
            
%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_3    = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(plot_condition,:)), 2, 'gaussian',500);
    %         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
            plot_data_mu_3 = mean(plot_data_3, 'omitnan');
            plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
            plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;
            plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

            [bl1, bf] = boundedline(1:length(plot_data_3), plot_data_mu_3, plot_data_sem, 'r'); hold on;
            bf.FaceAlpha = 0.7;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
            bl1.Color = coRegionColors(c,:)  ;
            bf.FaceColor = coRegionColors(c,:) ;
            bl1.MarkerFaceColor = bl1.Color;
            bf.FaceAlpha = 0.3;


            vl = vline([500 2600 4700 6800 9400]+500); 

            for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end

            hl = hline(baseline);
            hl.Color = 'k';
            hl.LineStyle = '-';
            hl.LineWidth = 1.0;


            title([regions{iRa}, ' <--> ', regions{iRb}])

            legend([bl1, bl2], sprintf('load3-%s', condNames{cond}), sprintf('load1-%s', condNames{cond}), 'location', 'best')
            
            pVal = nan(1, length(bins)-1);
            dVal = nan(2, length(bins)-1);
            plotTimesStats = nan(1, length(bins)-1);
            for iB = 1:length(bins)-1


                data1 = trapz(whole_trials_unit_cross_region_all{iRa,iRb}(~plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
                data1(isnan(data1)) = [];
                dVal(1,iB) = mean(data1);
                data3 = trapz(whole_trials_unit_cross_region_all{iRa,iRb}(plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
                data3(isnan(data3)) = [];
                dVal(2,iB) = mean(data3);

                plotTimesStats(iB) = bins(iB) + (binSz/2);

                if ~isempty(data1) && ~isempty(data3)
%                     [~,~,p_perm] = statcond({data1',data3'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
                    [pRank,H] = ranksum(data1, data3);
                    [H,pT,CI] = ttest2(data1,data3);
                    pVal(iB) = pRank;
                end







            end
            adj_p = nan(size(pVal));
            yVal = repmat(max(plot_data_mu_3), [1 length(pVal)]);
            [h, crit_p, adj_ci_cvrg, adj_p(~isnan(pVal))]=fdr_bh(pVal(~isnan(pVal)),0.05,'pdep','no');

            bnds = mask2bounds(adj_p < 0.05);

            for iBnd = 1:size(bnds,1)
                pl = plot(plotTimesStats(bnds(iBnd,1):bnds(iBnd,1)), yVal(bnds(iBnd,1):bnds(iBnd,1)), 'k*'); hold on;
                pl.LineWidth = 1.5;
            end
        
            subjID_1 = subjID(~plot_condition); subjID_3 = subjID(plot_condition);
            probeResp_1 = mean(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); subjID_1(isnan(probeResp_1)) = []; probeResp_1(isnan(probeResp_1)) = [];
            probeResp_3 = mean(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2); subjID_3(isnan(probeResp_3)) = []; probeResp_3(isnan(probeResp_3)) = [];
            
            coRcoF = mean(whole_trials_unit_cross_region_all{iRa,iRb}(:,taskMarkers(end):taskMarkers(end)+1.5e3), 2);
            coRcoF(cTrial:end) = [];
            
            respTab{iRa,iRb} = table(subjID, plot_condition, coRcoF, 'VariableNames', {'subject', 'condition', 'coR rates'});
            
            probeResp_1 = [subjID_1, probeResp_1];
            probeResp_3 = [subjID_3, probeResp_3];
            respData{1,cond, c} = probeResp_1; % respData(2,1,cond, c) = std(probeResp_1)/sqrt(length(probeResp_1)); 
            respData{2,cond, c} = probeResp_3;  %respData(2,2,cond, c) = std(probeResp_3)/sqrt(length(probeResp_3)); 

%             [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% 
%             fprintf('%s <--> %s cofire during %s --- response p = %.2f\n', regions{iRa}, regions{iRb}, condNames{cond}, pRank)
% 
%             probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%             probeResp_1(isnan(probeResp_1)) = [];
%             probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%             probeResp_3(isnan(probeResp_3)) = [];
% %             [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%             [pRank,H] = ranksum(probeResp_3, probeResp_1);
% 
%             fprintf('%s <--> %s cofire during %s --- maintenance p = %.2f\n', regions{iRa}, regions{iRb},  condNames{cond}, pRank)

            blall(iR) =bl;
        %     xlim(Lim)
        %     ylim([0.3 1])
            ylabel('co-firing')
            xlabel('time from trial start [ms]')

        end
    end

    fig= gcf;
    fig.Color = 'w';
%     %savepdf(gcf, fullfile(exportDirFigs, sprintf('coFireTaskAcrossRegions_%s.pdf', condNames{cond})))
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
%savepdf(gcf, fullfile(exportDirFigs, 'respLatency.pdf'))

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

figure('Position', [894.0000 135 224.0000 777.0000]);
h = boxplot(densityAll, 'Orientation', 'vertical', 'Symbol', '');

xlabel('ripple density')
% xlim([0 whiskerLimits(2)]);
fig = gcf;
fig.Color = 'w';
%savepdf(gcf, fullfile(exportDirFigs, 'density.pdf'))



























