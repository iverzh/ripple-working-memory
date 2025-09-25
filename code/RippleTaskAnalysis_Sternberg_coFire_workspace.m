

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
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

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

recordingState = 'wake';
location = 'NC';

regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
% regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end,3:end-1));   
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2,3:end-1));
bundleNotes = table2cell(bundleNotes(3:end,3:end-1));   


%%
for coFireWindow = round(25/2) %round([25 50 15 8 100 ]/2)

    computeCoRipple = false;
    rippleFR = nan(length(subj_list_full), 100);
    CoRippleRatesRegionAll = nan(length(regions), length(regions), length(subj_list_full));
    win = [-9 25] * 1e3;
    winResp = [-1.5 1] * 1e3;
    winLength = 12905; %ms
    coFireSubsample = 1;
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
    whole_trials_cross_region_dur = cell(length(regions), length(regions), 3);
    load_trials_all_unitPair = cell(length(regions), length(regions));
    unitPairID = cell(length(regions), length(regions));
    testCoFire = cell(length(regions), length(regions));

    for iRa = 1:length(regions)
        for iRb = iRa:length(regions)
            whole_trials_cross_region_all{iRa, iRb} = nan(6.5e3, winLength);       
            testCoFire{iRa, iRb} = zeros(1, ceil(winLength/coFireSubsample));       
            load_trials_all_unitPair{iRa, iRb} = nan(1, 3.5e5);    
            unitPairID{iRa, iRb} = nan(4.5e5, 3);

            for ii = 1:size(whole_trials_unit_cross_region_all,3)
    %             whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e5, ceil(winLength/coFireSubsample));
                whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e5, 4);
                whole_trials_cross_region_dur{iRa, iRb,ii} = nan(4.5e5, 4);
            end
        end
    end

    subjID = nan(6.5e3, 1);
    correct_trials_all = nan(6.5e3, 1);
    load_trials_all = nan(6.5e3, 1);
    probe_in_out = nan(6.5e3, 1);
    respLatency = nan(6.5e3, 1);
    densityAll = [];

    taskMarkers = [500 2600 4700 6800 9400]+500;
    uChanAll = [];
    tic
    cTrial = 1;
    cUnitPair = ones(length(regions), length(regions));
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
        trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));


        modifier = 'IISremv_singleWire_v3_z25';
        modifier = '1kHz_template_z25';
    %     modifier = 'IISremv_singleWire_v3';
    %     modifier = 'IISremv_v3';
        tag = [recordingState,'_',location,'_',modifier];
%         filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
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
        micrObj.lfp_data(:, end:end+3e3) = nan; %pad the ending


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
        erp_dat_region = cell(1,length(regions));
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
                     u_dat_region{iR} = sum(spikeMask(uRegion,:),1);
                elseif sum(uRegion) == 0
                     u_dat_region{iR} = nan(1, length(rippMask));
                end

                 if trialLoad == 3
                     

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




    %             if length(dat) > 16e3; continue; end
    %             trialData(1:length(dat)) = dat;
    %             trials_region_load1(iT,:,iR) = trialData;
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


                    uLFPregionA = uChan(contains(uLocations, regions(iRa)));
                    uLFPregionB = uChan(contains(uLocations, regions(iRb)));
                    uRegionA = find(contains(uLocations, regions(iRa)));
                    uRegionB = find(contains(uLocations, regions(iRb)));

                    ovA = nan(1, length(uLFPregionA));
                    ovAcell = arrayfun(@(X) find(LFPchanNum == X), uLFPregionA', 'UniformOutput', false);
                    ovA(cellfun(@(X) ~isempty(X), ovAcell)) = cell2mat(ovAcell);
                    ovB = nan(1, length(uLFPregionB));
                    ovBcell = arrayfun(@(X) find(LFPchanNum == X), uLFPregionB', 'UniformOutput', false);
                    ovB(cellfun(@(X) ~isempty(X), ovBcell)) = cell2mat(ovBcell);                
                    if sum(~isnan(ovA)) < 1 || sum(~isnan(ovB)) < 1; continue; end

        % 
        %                 datArip =  ripTrialFullChan(ovA(~isnan(ovA)),:);
        %                 datBrip =  ripTrialFullChan(ovB(~isnan(ovB)),:);

                    datArip =  ripTrialFullChan(contains(locations, regions(iRa)),:);
                    datBrip =  ripTrialFullChan(contains(locations, regions(iRb)),:);
        %                 
                    nanMaskA = isnan(sum(datArip,1)); bndA = mask2bounds(nanMaskA); nanMaskA = bounds2mask(bndA,length(nanMaskA), 100)'; 
                    nanMaskB = isnan(sum(datBrip,1)); bndB = mask2bounds(nanMaskB); nanMaskB = bounds2mask(bndB,length(nanMaskB), 100)';        
                    nanMask = nanMaskA | nanMaskB;
                    datArip(:,nanMaskA) = 0; datArip = logical(datArip);
                    datBrip(:,nanMaskB) = 0; datBrip = logical(datBrip);
                    datABcorip = sum(datArip) > 0 & sum(datBrip) > 0;
                    datABnorip = sum(datArip) == 0 & sum(datBrip) == 0;


                    nUnits = length(uLFPregionA) + length(uLFPregionB);

                    for uA = 1:length(uRegionA)
                        for uB = 1:length(uRegionB)
        %                         datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
        %                         datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                            datA =  unitTrialFullChan(uRegionA(uA),:); datB =  unitTrialFullChan(uRegionB(uB),:);
                            datA = sum(datA,1) > 0;
                            datB = sum(datB,1) > 0;

%                             datA = movsum(datA, [coFireWindow coFireWindow]);
%                             datB = movsum(datB, [coFireWindow coFireWindow]);
                            datAB = double((datA > 0) & (datB > 0));
            %                 datAB = datAB / sum(~nanMaskA & ~nanMaskB);
%                             datAB(nanMask) = nan;
                            
    %                         whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),:) = sumBins(datAB, coFireSubsample); 
                            whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),1) = countCoFire(datA, datB, 1:1e3, coFireWindow*2); % / (length(1:1e3) - sum(nanMask(1:1e3)));
                            whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),2) = countCoFire(datA, datB, taskMarkers(1):taskMarkers(2), coFireWindow*2); % /  (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),3) = countCoFire(datA, datB, taskMarkers(4):taskMarkers(5), coFireWindow*2);  % / (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
                            whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),4) = countCoFire(datA, datB, taskMarkers(5):taskMarkers(5)+1e3, coFireWindow*2); % / (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),1) = length(mask2bounds(datAB(1:1e3))); % / (length(1:1e3) - sum(nanMask(1:1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),2) = length(mask2bounds(datAB(taskMarkers(1):taskMarkers(2)))); % /  (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),3) = length(mask2bounds(datAB(taskMarkers(4):taskMarkers(5)))); % / (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),4) = length(mask2bounds(datAB(taskMarkers(5):taskMarkers(5)+1e3)));% / (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),1) = sum(datAB(1:1e3),'omitnan'); %/ (length(1:1e3) - sum(nanMask(1:1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),2) = sum(datAB(taskMarkers(1):taskMarkers(2)),'omitnan'); %/  (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),3) = sum(datAB(taskMarkers(4):taskMarkers(5)),'omitnan'); %/ (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 1}(cUnitPair(iRa, iRb),4) = sum(datAB(taskMarkers(5):taskMarkers(5)+1e3), 'omitnan');% / (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
                            whole_trials_cross_region_dur{iRa, iRb, 1}(cUnitPair(iRa, iRb),1) = (length(1:1e3) - sum(nanMask(1:1e3))); 
                            whole_trials_cross_region_dur{iRa, iRb, 1}(cUnitPair(iRa, iRb),2) = (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_cross_region_dur{iRa, iRb, 1}(cUnitPair(iRa, iRb),3) = (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
                            whole_trials_cross_region_dur{iRa, iRb, 1}(cUnitPair(iRa, iRb),4) = (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
                            testCoFire{iRa, iRb} = testCoFire{iRa, iRb} +  datAB;


                            % co fire during co-ripples
        %                         datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
        %                         datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                            datA =  unitTrialFullChan(uRegionA(uA),:); datB =  unitTrialFullChan(uRegionB(uB),:);
                            datA(:, ~datABcorip) = 0; datA = sum(datA, 1) > 0;
                            datB(:, ~datABcorip) = 0; datB = sum(datB, 1) > 0;
%                             datA(:,  sum(datArip) < 1) = 0; datA = sum(datA, 1) > 0;
%                             datB(:,  sum(datBrip) < 1) = 0; datB = sum(datB, 1) > 0;

%                             datA = movsum(datA, [coFireWindow coFireWindow]);
%                             datB = movsum(datB, [coFireWindow coFireWindow]);
%                             datAB = double((datA > 0) & (datB > 0));
            %                 datAB = datAB / sum(datArip & datBrip);

%                             datAB(nanMaskA | nanMaskB) = nan;


    %                         whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),:) = sumBins(datAB, coFireSubsample);
                            whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),1) = countCoFire(datA, datB, 1:1e3, coFireWindow*2); % / (length(1:1e3) - sum(nanMask(1:1e3)));
                            whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),2) = countCoFire(datA, datB, taskMarkers(1):taskMarkers(2), coFireWindow*2); % /  (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),3) = countCoFire(datA, datB, taskMarkers(4):taskMarkers(5), coFireWindow*2);  % / (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
                            whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),4) = countCoFire(datA, datB, taskMarkers(5):taskMarkers(5)+1e3, coFireWindow*2); % / (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),1) = length(mask2bounds(datAB(1:1e3))); % / (sum(datABcorip(1:1e3)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),2) = length(mask2bounds(datAB(taskMarkers(1):taskMarkers(2))));% / (sum(datABcorip(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),3) = length(mask2bounds(datAB(taskMarkers(4):taskMarkers(5)))); %/ (sum(datABcorip(taskMarkers(4):taskMarkers(5))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),4) = length(mask2bounds(datAB(taskMarkers(5):taskMarkers(5)+1e3))); % / (sum(datABcorip(taskMarkers(5):taskMarkers(5)+1e3))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),1) = sum(datAB(1:1e3),'omitnan'); %/ (sum(datABcorip(1:1e3)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),2) = sum(datAB(taskMarkers(1):taskMarkers(2)), 'omitnan');% / (sum(datABcorip(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),3) = sum(datAB(taskMarkers(4):taskMarkers(5)), 'omitnan'); %/ (sum(datABcorip(taskMarkers(4):taskMarkers(5))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 2}(cUnitPair(iRa, iRb),4) = sum(datAB(taskMarkers(5):taskMarkers(5)+1e3),'omitnan'); %/ (sum(datABcorip(taskMarkers(5):taskMarkers(5)+1e3))); 
                            whole_trials_cross_region_dur{iRa, iRb, 2}(cUnitPair(iRa, iRb),1) = (sum(datABcorip(1:1e3), 'omitnan')); 
                            whole_trials_cross_region_dur{iRa, iRb, 2}(cUnitPair(iRa, iRb),2) = (sum(datABcorip(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_cross_region_dur{iRa, iRb, 2}(cUnitPair(iRa, iRb),3) = (sum(datABcorip(taskMarkers(4):taskMarkers(5))));  
                            whole_trials_cross_region_dur{iRa, iRb, 2}(cUnitPair(iRa, iRb),4) = (sum(datABcorip(taskMarkers(5):taskMarkers(5)+1e3))); 
    

                            % co fire outside of any rippling
        %                         datA =  unitTrialFullChan(contains(uLocations, regions(iRa)),:); datA(isnan(ovA),:) = []; 
        %                         datB =  unitTrialFullChan(contains(uLocations, regions(iRb)),:); datB(isnan(ovB),:) = []; 
                            datA = unitTrialFullChan(uRegionA(uA),:); datB = unitTrialFullChan(uRegionB(uB),:);
                            datA(: , ~datABnorip) = 0; datA = sum(datA, 1) > 0;
                            datB(: , ~datABnorip) = 0; datB = sum(datB, 1) > 0;
%                             datA(:,  sum(datArip) > 0) = 0; datA = sum(datA, 1) > 0;
%                             datB(:,  sum(datBrip) > 0) = 0; datB = sum(datB, 1) > 0;
                            
%                             datA = movsum(datA, [coFireWindow coFireWindow]);
%                             datB = movsum(datB, [coFireWindow coFireWindow]);
%                             
%                             datAB = double((datA > 0) & (datB > 0));
            %                 datAB = datAB / sum(~datArip & ~datBrip & ~nanMaskA & ~nanMaskB);
%                             datAB(nanMaskA | nanMaskB) = nan;


    %                         whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),:) = sumBins(datAB, coFireSubsample);
                            whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),1) = countCoFire(datA, datB, 1:1e3, coFireWindow*2); % / (length(1:1e3) - sum(nanMask(1:1e3)));
                            whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),2) = countCoFire(datA, datB, taskMarkers(1):taskMarkers(2), coFireWindow*2); % /  (length(taskMarkers(1):taskMarkers(2)) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),3) = countCoFire(datA, datB, taskMarkers(4):taskMarkers(5), coFireWindow*2);  % / (length(taskMarkers(4):taskMarkers(5)) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
                            whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),4) = countCoFire(datA, datB, taskMarkers(5):taskMarkers(5)+1e3, coFireWindow*2); % / (length(taskMarkers(5):taskMarkers(5)+1e3)- sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),1) = length(mask2bounds(datAB(1:1e3))); % / (sum(datABnorip(1:1e3))- sum(nanMask(1:1e3))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),2) = length(mask2bounds(datAB(taskMarkers(1):taskMarkers(2)))); % / (sum(datABnorip(taskMarkers(1):taskMarkers(2))) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),3) = length(mask2bounds(datAB(taskMarkers(4):taskMarkers(5)))); % / (sum(datABnorip(taskMarkers(4):taskMarkers(5))) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),4) = length(mask2bounds(datAB(taskMarkers(5):taskMarkers(5)+1e3))); % / (sum(datABnorip(taskMarkers(5):taskMarkers(5)+1e3)) - sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),1) = sum(datAB(1:1e3),'omitnan'); %/ (sum(datABnorip(1:1e3))- sum(nanMask(1:1e3))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),2) = sum(datAB(taskMarkers(1):taskMarkers(2)),'omitnan'); %/ (sum(datABnorip(taskMarkers(1):taskMarkers(2))) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),3) = sum(datAB(taskMarkers(4):taskMarkers(5)),'omitnan'); %/ (sum(datABnorip(taskMarkers(4):taskMarkers(5))) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
%                             whole_trials_unit_cross_region_all{iRa, iRb, 3}(cUnitPair(iRa, iRb),4) = sum(datAB(taskMarkers(5):taskMarkers(5)+1e3),'omitnan'); %/ (sum(datABnorip(taskMarkers(5):taskMarkers(5)+1e3)) - sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));                            
                            whole_trials_cross_region_dur{iRa, iRb, 3}(cUnitPair(iRa, iRb),1) = (sum(datABnorip(1:1e3))- sum(nanMask(1:1e3))); 
                            whole_trials_cross_region_dur{iRa, iRb, 3}(cUnitPair(iRa, iRb),2) = (sum(datABnorip(taskMarkers(1):taskMarkers(2))) - sum(nanMask(taskMarkers(1):taskMarkers(2))));  
                            whole_trials_cross_region_dur{iRa, iRb, 3}(cUnitPair(iRa, iRb),3) = (sum(datABnorip(taskMarkers(4):taskMarkers(5))) - sum(nanMask(taskMarkers(4):taskMarkers(5)))); 
                            whole_trials_cross_region_dur{iRa, iRb, 3}(cUnitPair(iRa, iRb),4) = (sum(datABnorip(taskMarkers(5):taskMarkers(5)+1e3)) - sum(nanMask(taskMarkers(5):taskMarkers(5)+1e3)));

                            load_trials_all_unitPair{iRa, iRb}(cUnitPair(iRa, iRb)) = trialLoad;
                            unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),1) = subj;
                            unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),2:3) = [uA uB];
                            cUnitPair(iRa, iRb) = cUnitPair(iRa, iRb) + 1;
                        end
                    end
                end
            end
            
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
    save(sprintf('Sternbger_ripple_results_z25_%icoF_splitHem_clip_v2.mat', coFireWindow*2), '-v7.3')
end
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
% close all
regionColors =  brewermap(12, 'Dark2');
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  


% plot_condition = true(1, length(correct_trials_all));
plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');

coRegionColors = nan(10,3);
for iRa = 1:length(regionsPlot)
    for iRb = iRa+1:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        
        plot_condition_coF = [];
        subjID_uPair = [];
        uPairsID  =[];
        coRresp = [];
        noRresp = [];
        resp = [];
        coRmaintenance = [];
        noRmaintenance = [];
        maintenance = [];
        baseline = [];
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                
                if isempty(testCoFire{iiA, iiB}); continue; end
                testAll = testAll + testCoFire{iiA, iiB};
                
            end
        end
        
        figure(10);
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)
%         plot(1, sum(testAll(1:10)/length(1:1e1)), 'o'); hold on;
%         plot(2, sum(testAll(101:110))/length(taskMarkers(5):taskMarkers(5)+1e1), 'o'); hold on;
        plot(testAll)
%         xlim([0 3])
        
        
        
        
        
        
        c = c+1;

    end
end
% taskMarkers = [[500 2600 4700 6800 9400]+500]/100;
%%

close all
regionColors =  brewermap(12, 'Dark2');

regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

condNames = {'All', 'coRip', 'noRip'};
respData = cell(2, length(condNames), 10);
% Find all open figure objects
openFigures = findobj('Type', 'figure');
% Get the number of open figures
numOpenFigures = numel(openFigures);

tabCoFire = cell(5);

figure('Position', [95 140 1500 1120]);
% c = 1;
coRegionColors = nan(10,3);
Pval = nan(length(regionsPlot), length(regionsPlot), 3);
for iRa = 1:length(regionsPlot)
    for iRb = iRa+1:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        
        plot_condition_coF = [];
        subjID_uPair = [];
        uPairsID  =[];
        coRimg = [];
        noRimg = [];
        img = [];
        coRresp = [];
        noRresp = [];
        resp = [];
        coRmaintenance = [];
        noRmaintenance = [];
        maintenance = [];
        baseline = [];
        
        for iiA = iRapl
            for iiB = iRbpl
                plot_condition_T = load_trials_all_unitPair{iiA, iiB}; plot_condition_T(isnan(plot_condition_T)) = [];
                plot_condition_T = plot_condition_T==3;
                
                if isempty(unitPairID{iiA, iiB}); continue; end
                subjID_uPairT = unitPairID{iiA, iiB}(:,1);

                coRimgT         = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,3) ./ whole_trials_cross_region_dur{iiA, iiB,2}(:,3) * 1e3;
                noRimgT         = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,3) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,3) * 1e3;
                imgT            = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,3) ./ whole_trials_cross_region_dur{iiA, iiB,1}(:,3) * 1e3;
                coRrespT        = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,4) ./ whole_trials_cross_region_dur{iiA, iiB,2}(:,4) * 1e3;
                noRrespT        = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,4) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;
                respT           = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,4) ./ whole_trials_cross_region_dur{iiA, iiB,1}(:,4) * 1e3;
                coRmaintenanceT = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,2) ./ whole_trials_cross_region_dur{iiA, iiB,2}(:,2) * 1e3;
                noRmaintenanceT = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,2) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,2) * 1e3;
                maintenanceT    = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,2)./ whole_trials_cross_region_dur{iiA, iiB,1}(:,2) * 1e3;
                baselineT       = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,1) ./ whole_trials_cross_region_dur{iiA, iiB,1}(:,1) * 1e3;
                
                coRimgT(cUnitPair(iiA, iiB):end) = [];
                noRimgT(cUnitPair(iiA, iiB):end) = [];
                imgT(cUnitPair(iiA, iiB):end) = [];
                noRrespT(cUnitPair(iiA, iiB):end) = [];
                respT(cUnitPair(iiA, iiB):end) = [];
                coRrespT(cUnitPair(iiA, iiB):end) = [];
                noRrespT(cUnitPair(iiA, iiB):end) = [];
                respT(cUnitPair(iiA, iiB):end) = [];
                coRmaintenanceT(cUnitPair(iiA, iiB):end) = [];
                noRmaintenanceT(cUnitPair(iiA, iiB):end) = [];
                maintenanceT(cUnitPair(iiA, iiB):end) = [];
                baselineT(cUnitPair(iiA, iiB):end) = [];
                subjID_uPairT(cUnitPair(iiA, iiB):end) = [];
                
                plot_condition_coF = [plot_condition_coF plot_condition_T];
                subjID_uPair = [subjID_uPair subjID_uPairT'];
                
                baseline = [baseline baselineT'];

                coRimg = [coRimg coRimgT'];
                noRimg = [noRimg noRimgT'];
                img = [img imgT'];
                
                coRmaintenance = [coRmaintenance coRmaintenanceT'];
                noRmaintenance = [noRmaintenance noRmaintenanceT'];
                maintenance = [maintenance maintenanceT'];
                
                coRresp = [coRresp coRrespT'];
                noRresp = [noRresp noRrespT'];
                resp = [resp respT'];
                
                uPairs = unitPairID{iiA, iiB}; uPairs(isnan(uPairs(:,1)),:) = [];
                uPairsIDT = nan(1, size(uPairs,1));
                for iN = 1:length(uPairsIDT)
                    rowStr = sprintf('%04d', uPairs(iN,:));
                    uPairsIDT(iN) = str2double(rowStr);
                end
                uPairsID =  [uPairsID uPairsIDT];
            end
        end
        
        
        
        
%         
%         
%         coRresp = coRresp(~isnan(coRresp)); coRresp(~isfinite(coRresp)) = 0;
%         noRresp = noRresp(~isnan(coRresp)); noRresp(~isfinite(noRresp)) = 0;
%         resp = resp(~isnan(coRresp)); resp(~isfinite(resp)) = 0;
%         coRmaintenance = coRmaintenance(~isnan(coRresp)); coRmaintenance(~isfinite(coRmaintenance)) = 0;
%         noRmaintenance = noRmaintenance(~isnan(coRresp)); noRmaintenance(~isfinite(noRmaintenance)) = 0;
%         maintenance = maintenance(~isnan(coRresp)); maintenance(~isfinite(maintenance)) = 0;
%         baseline = baseline(~isnan(coRresp)); baseline(~isfinite(baseline)) = 0;
%         subjID_uPair = subjID_uPair(~isnan(coRresp));
%         uPairsID = uPairsID(~isnan(coRresp));

        coRresp(isinf(coRresp)) = 0;
        noRresp(isinf(noRresp)) = 0;
        resp(isinf(resp)) = 0;
        coRimg(isinf(coRimg)) = 0;
        noRimg(isinf(noRimg)) = 0;
        img(isinf(img)) = 0;
        coRmaintenance(isinf(coRmaintenance)) = 0;
        noRmaintenance(isinf(noRmaintenance)) = 0;
        maintenance(isinf(maintenance)) = 0;
        baseline(isinf(baseline)) = 0;
        
        
        
        

        
%         tabCoFire{iRa,iRb} = table(subjID_uPair, plot_condition_coF(~isnan(coRresp)), coRresp,coRmaintenance,baseline, ...
%                 'VariableNames', {'subject',    'condition',         'rates_resp','rates_maintenance','rates_baseline'});






%         tID = tabCoFire{iRa,iRb}.subject;  
        
        tCond = plot_condition_coF(~isnan(coRresp));  

          



        tRespID = [arrayfun(@(X) mean(baseline(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRimg(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRimg(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(img(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRmaintenance(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRmaintenance(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(maintenance(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRresp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRresp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(resp(uPairsID == X), 'omitnan'), unique(uPairsID))'] ;
               
        bMu = mean(tRespID(:,1)); tRespID = tRespID/bMu;

% %                
%        tRespID = [baseline', ...
%                   coRmaintenance', ...    
%                   noRmaintenance', ...
%                   maintenance', ...
%                   coRresp', ...
%                   noRresp', ...
%                   resp'] ;
               
%         tRespID(sum(tRespID,2) ==0, :) = [];
%         tRespID(tRespID(:,1) ==0, :) = [];
%         tRespID(~isfinite(tRespID)) = 0;
%             tRespID z= 
        figure(numOpenFigures+1); 
        c = (iRb - 1) * 5 + iRa;

%         subplot(5,5,c)            
        subplot(5,5,c)     
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
%             vs(2).BoxColor = coRegionColors(c,:);

%             vs = violinplot(tRespID, {'load1', 'load3'}, 'ShowData', false, 'EdgeColor', [1 1 1], 'ViolinAlpha', 0, 'BoxWidth', 0.1); hold on;
%         er = errorbar([1:7], mean(tRespID), std(tRespID) / sqrt(size(tRespID,1)), 'o'); hold on;
        errobarPatch(1, tRespID(:,1), 'k', 0.3);
        errobarPatch(2:3:10, tRespID(:,2:3:10), coRegionColors(c,:), 0.3);
        errobarPatch(3:3:10, tRespID(:,3:3:10), coRegionColors(c,:)*0.6, 0.3);
        errobarPatch(4:3:10, tRespID(:,4:3:10), coRegionColors(c,:)*0.9, 0.3);
        errobarPatch(7, tRespID(:,7), coRegionColors(c,:)*0.8, 0.3);
        
        vline([1.5 4.5 7.5])
%         hl = hline(1); hl.LineWidth = 0.5; hl.Color = 'k'; hl.LineStyle = '-';
        
       


        [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,10), 'tail', 'left');
        Pval(iRa,iRb,1) = P;
        [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,9), 'tail', 'left');
        Pval(iRa,iRb,2) = P;
        [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,8), 'tail', 'left');
        Pval(iRa,iRb,3) = P;
        
        
        fprintf('%s <--> %s       p = %.4f\n', regions{iRa},regions{iRb}, P)

        
%         if P < 0.05; plot(1.5, mean(tRespID(:)), 'k*'); hold on; end

        yLim = quantile(tRespID(:), [0.05 0.95]);
        ax = gca;
        ax.XTick = [1 3 6 9];
        ax.XTickLabel = {'bsl.', 'img.','maint.', 'resp'};
%             ax.YLim = [-0.0 yLim(2)];
        title([regionsPlot{iRa}, ' <--> ', regionsPlot{iRb}])
        xlim([0 11])
        ylabel('co-fire [obs / bsl.]')
        
        c = c+1;

    end
end
fig = gcf;
fig.Color = 'w';
% savepdf(gcf, fullfile(exportDirFigs, sprintf('coFireTaskAcrossRegions_%i.pdf', coFireWindow*2)))


adj_p = nan(size(Pval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(Pval))]=fdr_bh(Pval(~isnan(Pval)),0.05,'pdep','yes');

conds = {'all', 'noR', 'coR'};
figure('Position', [1205 1029 1486 292]); 
for ii = 1:size(adj_p, 3)
    subplot(1,3,ii)
    imagesc(adj_p(:,:,ii)', [min(adj_p(:))-1e-2 1]); hold on;
    colorbar;
    hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
    vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
    cmap = flipud(slanCM('viridis'));
    cmap = [1 1 1; cmap];
    colormap(cmap)
    [yy xx] = find(adj_p(:,:,ii)' < 0.05);
    plot(xx, yy, 'r*'); hold on;



    ax = gca;
    ax.YTick = 1:size(adj_p,1);
    ax.XTick = 1:size(adj_p,1);
    ax.XTickLabel = regions;
    ax.YTickLabel = regions;
    
    title(conds{ii})



end
sgtitle('co firing resp. vs baseline')
fig= gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coFireTaskAcrossRegions_%i_stats.pdf', coFireWindow*2)))

%%


close all

openFigures = findobj('Type', 'figure');
    % Get the number of open figures
numOpenFigures = numel(openFigures);
% figure('Position', [95 140 1500 1120]);
% c = 1;
coRegionColors = nan(10,3);
Pval = nan(length(regionsPlot), length(regionsPlot), 3);
for iRa = 1:length(regionsPlot)
    for iRb = iRa+1:length(regionsPlot)
        figure('Position', [958 686 313 223]);
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        
        plot_condition_coF = [];
        subjID_uPair = [];
        uPairsID  =[];
        coRimg = [];
        noRimg = [];
        img = [];
        coRresp = [];
        noRresp = [];
        resp = [];
        coRmaintenance = [];
        noRmaintenance = [];
        maintenance = [];
        baseline = [];
        
        for iiA = iRapl
            for iiB = iRbpl
                plot_condition_T = load_trials_all_unitPair{iiA, iiB}; plot_condition_T(isnan(plot_condition_T)) = [];
                plot_condition_T = plot_condition_T==3;
                
                if isempty(unitPairID{iiA, iiB}); continue; end
                subjID_uPairT = unitPairID{iiA, iiB}(:,1);

                coRimgT = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,3);
                noRimgT = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,3);
                imgT = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,3);
                coRrespT = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,4);
                noRrespT = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,4);

                
                respT = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,4);
                coRmaintenanceT = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,2);
                noRmaintenanceT = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,2);
                maintenanceT = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,2);
                baselineT = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,1);
                
                coRimgT(cUnitPair(iiA, iiB):end) = [];
                noRimgT(cUnitPair(iiA, iiB):end) = [];
                imgT(cUnitPair(iiA, iiB):end) = [];
                noRrespT(cUnitPair(iiA, iiB):end) = [];
                respT(cUnitPair(iiA, iiB):end) = [];
                coRrespT(cUnitPair(iiA, iiB):end) = [];
                noRrespT(cUnitPair(iiA, iiB):end) = [];
                respT(cUnitPair(iiA, iiB):end) = [];
                coRmaintenanceT(cUnitPair(iiA, iiB):end) = [];
                noRmaintenanceT(cUnitPair(iiA, iiB):end) = [];
                maintenanceT(cUnitPair(iiA, iiB):end) = [];
                baselineT(cUnitPair(iiA, iiB):end) = [];
                subjID_uPairT(cUnitPair(iiA, iiB):end) = [];
                
                plot_condition_coF = [plot_condition_coF plot_condition_T];
                subjID_uPair = [subjID_uPair subjID_uPairT'];
                
                baseline = [baseline baselineT'];

                coRimg = [coRimg coRimgT'];
                noRimg = [noRimg noRimgT'];
                img = [img imgT'];
                
                coRmaintenance = [coRmaintenance coRmaintenanceT'];
                noRmaintenance = [noRmaintenance noRmaintenanceT'];
                maintenance = [maintenance maintenanceT'];
                
                coRresp = [coRresp coRrespT'];
                noRresp = [noRresp noRrespT'];
                resp = [resp respT'];
                
                uPairs = unitPairID{iiA, iiB}; uPairs(isnan(uPairs(:,1)),:) = [];
                uPairsIDT = nan(1, size(uPairs,1));
                for iN = 1:length(uPairsIDT)
                    rowStr = sprintf('%04d', uPairs(iN,:));
                    uPairsIDT(iN) = str2double(rowStr);
                end
                uPairsID =  [uPairsID uPairsIDT];
            end
        end
        
        
        
        
%         
%         
%         coRresp = coRresp(~isnan(coRresp)); coRresp(~isfinite(coRresp)) = 0;
%         noRresp = noRresp(~isnan(coRresp)); noRresp(~isfinite(noRresp)) = 0;
%         resp = resp(~isnan(coRresp)); resp(~isfinite(resp)) = 0;
%         coRmaintenance = coRmaintenance(~isnan(coRresp)); coRmaintenance(~isfinite(coRmaintenance)) = 0;
%         noRmaintenance = noRmaintenance(~isnan(coRresp)); noRmaintenance(~isfinite(noRmaintenance)) = 0;
%         maintenance = maintenance(~isnan(coRresp)); maintenance(~isfinite(maintenance)) = 0;
%         baseline = baseline(~isnan(coRresp)); baseline(~isfinite(baseline)) = 0;
%         subjID_uPair = subjID_uPair(~isnan(coRresp));
%         uPairsID = uPairsID(~isnan(coRresp));

        coRresp(isinf(coRresp)) = 0;
        noRresp(isinf(noRresp)) = 0;
        resp(isinf(resp)) = 0;
        coRmaintenance(isinf(coRmaintenance)) = 0;
        noRmaintenance(isinf(noRmaintenance)) = 0;
        maintenance(isinf(maintenance)) = 0;
        baseline(isinf(baseline)) = 0;
        img(isinf(img)) = 0;
        
        
        
        

        
%         tabCoFire{iRa,iRb} = table(subjID_uPair, plot_condition_coF(~isnan(coRresp)), coRresp,coRmaintenance,baseline, ...
%                 'VariableNames', {'subject',    'condition',         'rates_resp','rates_maintenance','rates_baseline'});






%         tID = tabCoFire{iRa,iRb}.subject;  
        
        tCond = plot_condition_coF(~isnan(coRresp));  

          

        plot_condition_coF = logical(plot_condition_coF);
        

%         tRespID = [arrayfun(@(X) mean(baseline(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(img(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(img(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(maintenance(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(maintenance(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(resp(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(resp(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))'] ;
               
       
       mu = [arrayfun(@(X) mean(resp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(resp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRresp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X)  mean(noRresp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRresp(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X)  mean(coRresp(uPairsID == X), 'omitnan'), unique(uPairsID))'] ;
                             
       tRespID = [arrayfun(@(X) mean(resp(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(resp(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRresp(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
                   arrayfun(@(X) mean(noRresp(uPairsID == X & ~plot_condition_coF), 'omitnan') , unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRresp(uPairsID == X & plot_condition_coF), 'omitnan') , unique(uPairsID))', ...
                   arrayfun(@(X) mean(coRresp(uPairsID == X & ~plot_condition_coF), 'omitnan') , unique(uPairsID))'] ;
%                
%       tRespID = [arrayfun(@(X) mean(maintenance(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(maintenance(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRmaintenance(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRmaintenance(uPairsID == X & ~plot_condition_coF), 'omitnan') , unique(uPairsID))', ...
%                    arrayfun(@(X) mean(coRmaintenance(uPairsID == X & plot_condition_coF), 'omitnan') , unique(uPairsID))', ...
%                    arrayfun(@(X) mean(coRmaintenance(uPairsID == X & ~plot_condition_coF), 'omitnan') , unique(uPairsID))'] ;
               
        mu(tRespID(:,1) == 0 & tRespID(:,2) == 0, :) = [];
%         tRespID(tRespID(:,6) == 0 & tRespID(:,5) == 0, :) = [];
        tRespID(tRespID(:,1) == 0 & tRespID(:,2) == 0, :) = [];
%         tRespID(tRespID(:,1) < 1 | tRespID(:,2) < 1, :) = [];
%         tRespID = tRespID ./ mu;
%         tRespID(isnan(tRespID)) = 0;
               
%        tRespID = [arrayfun(@(X) mean(baseline(uPairsID == X), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRimg(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRimg(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRmaintenance(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRmaintenance(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRresp(uPairsID == X & plot_condition_coF), 'omitnan'), unique(uPairsID))', ...
%                    arrayfun(@(X) mean(noRresp(uPairsID == X & ~plot_condition_coF), 'omitnan'), unique(uPairsID))'] ;
               
        bMu = mean(tRespID(:,1)); %tRespID = tRespID/bMu;

% %                
%        tRespID = [baseline', ...
%                   coRmaintenance', ...    
%                   noRmaintenance', ...
%                   maintenance', ...
%                   coRresp', ...
%                   noRresp', ...
%                   resp'] ;
               
%         tRespID(sum(tRespID,2) ==0, :) = [];

%         tRespID(tRespID(:,1) ==0, :) = [];
%         tRespID(~isfinite(tRespID)) = 0;
%             tRespID z= 
%         figure(numOpenFigures+1); 
        c = (iRb - 1) * 5 + iRa;

%         subplot(5,5,c)            
        subplot(5,5,c)     
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
%             vs(2).BoxColor = coRegionColors(c,:);

%             vs = violinplot(tRespID, {'load1', 'load3'}, 'ShowData', false, 'EdgeColor', [1 1 1], 'ViolinAlpha', 0, 'BoxWidth', 0.1); hold on;
%         er = errorbar([1:7], mean(tRespID), std(tRespID) / sqrt(size(tRespID,1)), 'o'); hold on;

        subplot(1,3,1)
        errobarPatch(1, tRespID(:,1), coRegionColors(c,:), 0.3);
        errobarPatch(2, tRespID(:,2), [ 0 0 0], 0.3);
        ylabel('co-fire per trial')
        yLim = quantile(tRespID(:,1:2), [0.35 0.80], 'all');
        ylim(yLim)
        ax = gca;
        ax.XTick = 1.5;
        ax.XTickLabel = sprintf('all coF');
        

        subplot(1,3,2)
        errobarPatch(3, tRespID(:,3), coRegionColors(c,:), 0.3);
        errobarPatch(4, tRespID(:,4), [ 0 0 0], 0.3);
        yLim = quantile(tRespID(:,3:4), [0.40 0.80], 'all');
        ylim(yLim)
        ax = gca;
        ax.XTick = 3.5;
        ax.XTickLabel = sprintf('noR coF');

        subplot(1,3,3)
        errobarPatch(5, tRespID(:,5), coRegionColors(c,:), 0.3);
        errobarPatch(6, tRespID(:,6), [ 0 0 0], 0.3);
        yLim = quantile(tRespID(:,5:6), [0.40 0.80], 'all');
        ylim(yLim)
        ax = gca;
        ax.XTick = 5.5;
        ax.XTickLabel = sprintf('coR coF');

        
%         vline([1.5 4.5 7.5])
%         hl = hline(1); hl.LineWidth = 0.5; hl.Color = 'k'; hl.LineStyle = '-';
        
       

% % 
%         [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,2), 'tail', 'right');
%         Pval(iRa,iRb,1) = P;
%         [H,P,CI,STATS] = ttest(tRespID(:,3), tRespID(:,4), 'tail', 'right');
%         Pval(iRa,iRb,2) = P;
%         [H,P,CI,STATS] = ttest(tRespID(:,5), tRespID(:,6), 'tail', 'right');
%         Pval(iRa,iRb,3) = P;
        
%         [P, obs, g] = permutationTest(tRespID(:,1), tRespID(:,2), 10000, 'sidedness', 'larger');
%         Pval(iRa,iRb,1) = P;
%         [P, obs, g] = permutationTest(tRespID(:,3), tRespID(:,4), 10000, 'sidedness', 'larger');
%         Pval(iRa,iRb,2) = P;
%         [P, obs, g] = permutationTest(tRespID(:,5), tRespID(:,6), 10000, 'sidedness', 'larger');
%         Pval(iRa,iRb,3) = P;
        
%         [P,H,STATS] = signrank(tRespID(:,1), tRespID(:,2), 'tail', 'right');
%         Pval(iRa,iRb,1) = P;
%         [P,H,STATS] = signrank(tRespID(:,3), tRespID(:,4), 'tail', 'right');
%         Pval(iRa,iRb,2) = P;
%         [P,H,STATS] = signrank(tRespID(:,5), tRespID(:,6), 'tail', 'right');
%         Pval(iRa,iRb,3) = P;
        
        [stats, df, P, surrog] = statcond({tRespID(:,1)', tRespID(:,2)'},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
        Pval(iRa,iRb,1) = P;
        [stats, df, P, surrog] = statcond({tRespID(:,3)', tRespID(:,4)'},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
        Pval(iRa,iRb,2) = P;
        [stats, df, P, surrog] = statcond({tRespID(:,5)', tRespID(:,6)'},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
        Pval(iRa,iRb,3) = P;

        fprintf('%s <--> %s       p = %.4f\n', regions{iRa},regions{iRb}, P)

        
%         if P < 0.05; plot(1.5, mean(tRespID(:)), 'k*'); hold on; end

        yLim = quantile(tRespID(:), [0.05 0.95]);
        ax = gca;
%         ax.XTick = [1 3 6 9];
%         ax.XTickLabel = {'bsl.', 'img.','maint.', 'resp'};
%             ax.YLim = [-0.0 yLim(2)];
%         title([regionsPlot{iRa}, ' <--> ', regionsPlot{iRb}])
        sgtitle(sprintf('%i pairs', size(tRespID,1)))
%         xlim([0 11])
        
        c = c+1;
        
        fig = gcf;
        fig.Color = 'w';
        
        savepdf(gcf, fullfile(exportDirFigs, sprintf('coFirePerTrial_%i_%s-%s.pdf', coFireWindow*2, regionsPlot{iRa}, regionsPlot{iRb})))


    end
    

end


adj_p = nan(size(Pval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(Pval))]=fdr_bh(Pval(~isnan(Pval)),0.05,'pdep','yes');

conds = {'all', 'no-rip', 'co-rip'};
figure('Position', [1205 1029 1486 292]); 
for ii = 1:size(adj_p, 3)
    subplot(1,3,ii)
    imagesc(adj_p(:,:,ii)', [min(adj_p(:))-1e-2 1]); hold on;
    colorbar;
    hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
    vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
    cmap = flipud(slanCM('viridis'));
    cmap = [1 1 1; cmap];
    colormap(cmap)
    [yy xx] = find(adj_p(:,:,ii)' < 0.05);
    plot(xx, yy, 'r*'); hold on;



    ax = gca;
    ax.YTick = 1:size(adj_p,1);
    ax.XTick = 1:size(adj_p,1);
    ax.XTickLabel = regions;
    ax.YTickLabel = regions;
    
    title(conds{ii})



end
% sgtitle('co firing load 3 vs load 1')
fig= gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirFigs, sprintf('coFirePerTrialAcrossRegions_%i_stats.pdf', coFireWindow*2)))

%%
% close all
taskMarkers = [[500 2600 4700 6800 9400]+500]/100;

condNames = {'All', 'coRip', 'noRip'};
respData = cell(2, length(condNames), 10);
for cond = 1:3
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
            plot_condition_coF = load_trials_all_unitPair{iRa, iRb}; plot_condition_coF(isnan(plot_condition_coF)) = [];
            plot_condition_coF = plot_condition_coF==3;
            figure(numOpenFigures+1)
            c = (iRb - 1) * 5 + iRa;
            subplot(5,5,c)            
%             baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:1e3), 'all', 'omitnan')/ ...
%                              mean(whole_trials_cross_region_all{iRa,iRb}(:,1:1e3)>0, 'all', 'omitnan');
                         
            baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:10), 'all', 'omitnan');


%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(~plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_1 = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(~plot_condition_coF,:))/baseline, 2, 'gaussian',5);
    %         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
            plot_data_mu_1 = mean(plot_data_1, 'omitnan');
            plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
%             plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
%             plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
            [bl2, bf] = boundedline(1:size(plot_data_1,2)-1, plot_data_mu_1(1:end-1), plot_data_sem(1:end-1), 'k'); hold on;
            bf.FaceAlpha = 0.7;
            
%             ripCond = smoothdata(mean(whole_trials_cross_region_all{iRa,iRb}(plot_condition, :)>0, 'omitnan'), 'gaussian' , 500);
            plot_data_3    = smoothdata((whole_trials_unit_cross_region_all{iRa,iRb,cond}(plot_condition_coF,:))/baseline, 2, 'gaussian',5);
    %         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
            plot_data_mu_3 = mean(plot_data_3, 'omitnan');
            plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
%             plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;
%             plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

            [bl1, bf] = boundedline(1:size(plot_data_3,2)-1, plot_data_mu_3(1:end-1), plot_data_sem(1:end-1), 'r'); hold on;
            bf.FaceAlpha = 0.7;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
            bl1.Color = coRegionColors(c,:)  ;
            bf.FaceColor = coRegionColors(c,:) ;
            bl1.MarkerFaceColor = bl1.Color;
            bf.FaceAlpha = 0.3;


            vl = vline(taskMarkers); 

            for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end

%             hl = hline(baseline);
            hl = hline(1);
            hl.Color = 'k';
            hl.LineStyle = '-';
            hl.LineWidth = 1.0;


            title([regions{iRa}, ' <--> ', regions{iRb}])

%             legend([bl1, bl2], sprintf('load3-%s', condNames{cond}), sprintf('load1-%s', condNames{cond}), 'location', 'best')
            
%             pVal = nan(1, length(bins)-1);
%             dVal = nan(2, length(bins)-1);
%             plotTimesStats = nan(1, length(bins)-1);
%             for iB = 1:length(bins)-1
% 
% 
%                 data1 = trapz(whole_trials_unit_cross_region_all{iRa,iRb}(~plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
%                 data1(isnan(data1)) = [];
%                 dVal(1,iB) = mean(data1);
%                 data3 = trapz(whole_trials_unit_cross_region_all{iRa,iRb}(plot_condition,times >= bins(iB) & times < bins(iB+1)), 2); 
%                 data3(isnan(data3)) = [];
%                 dVal(2,iB) = mean(data3);
% 
%                 plotTimesStats(iB) = bins(iB) + (binSz/2);
% 
%                 if ~isempty(data1) && ~isempty(data3)
% %                     [~,~,p_perm] = statcond({data1',data3'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%                     [pRank,H] = ranksum(data1, data3);
%                     [H,pT,CI] = ttest2(data1,data3);
%                     pVal(iB) = pRank;
%                 end
% 
% 
% 
% 
% 
% 
% 
%             end
%             adj_p = nan(size(pVal));
%             yVal = repmat(max(plot_data_mu_3), [1 length(pVal)]);
%             [h, crit_p, adj_ci_cvrg, adj_p(~isnan(pVal))]=fdr_bh(pVal(~isnan(pVal)),0.05,'pdep','no');
% 
%             bnds = mask2bounds(adj_p < 0.05);
% 
%             for iBnd = 1:size(bnds,1)
%                 pl = plot(plotTimesStats(bnds(iBnd,1):bnds(iBnd,1)), yVal(bnds(iBnd,1):bnds(iBnd,1)), 'k*'); hold on;
%                 pl.LineWidth = 1.5;
%             end
        
            subjID_uPair = unitPairID{iRa, iRb}(:,1);
            
            coRresp = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,taskMarkers(end):taskMarkers(end)+10), 2, 'omitnan');
            coRmaintenance = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,taskMarkers(end-1):taskMarkers(end)), 2, 'omitnan');
            baseline = mean(whole_trials_unit_cross_region_all{iRa,iRb,cond}(:,1:taskMarkers(1)), 2, 'omitnan');
            coRresp(cUnitPair(iRa, iRb):end) = [];
            coRmaintenance(cUnitPair(iRa, iRb):end) = [];
            baseline(cUnitPair(iRa, iRb):end) = [];
            subjID_uPair(cUnitPair(iRa, iRb):end) = [];
            
            tabCoFire{iRa,iRb} = table(subjID_uPair(~isnan(coRresp)), plot_condition_coF(~isnan(coRresp))', coRresp(~isnan(coRresp)),coRmaintenance(~isnan(coRresp)),baseline(~isnan(coRresp)), ...
                    'VariableNames', {'subject',               'condition',                     'rates_resp',            'rates_maintenance',            'rates_baseline'});
         
        
            
            
           


            blall(iR) =bl;
            xlim([0 130])
        %     ylim([0.3 1])
            ylabel('co-firing')
            xlabel('time from trial start [ms]')
            
            tID = tabCoFire{iRa,iRb}.subject;  
            tResp = tabCoFire{iRa,iRb}.rates_resp;  
            tCond = plot_condition_coF;  
            
            uPairs = unitPairID{iRa, iRb}; uPairs(isnan(uPairs(:,1)),:) = [];
            uPairsID = nan(1, size(uPairs,1));
            for iN = 1:length(uPairsID)
                rowStr = sprintf('%d', uPairs(iN,:));
                uPairsID(iN) = str2double(rowStr);
            end
            

            tRespID = [arrayfun(@(X) mean(tResp(uPairsID == X & tCond == 0)), unique(uPairsID))',arrayfun(@(X) mean(tResp(uPairsID == X & tCond == 1)), unique(uPairsID))'] ;
            tRespID(sum(tRespID,2) ==0, :) = [];
%             tRespID z= 
            figure(numOpenFigures+2); 
            subplot(5,5,c)            
            
            vs = violinplot(tRespID, {'load1', 'load3'}, 'ShowData', false, 'EdgeColor', [1 1 1], 'ViolinAlpha', 0, 'BoxWidth', 0.1); hold on;
            vs(2).BoxColor = coRegionColors(c,:);
            
            [H,P,CI,STATS] = ttest(tRespID(:,1), tRespID(:,2), 'tail', 'left');
            
            if P < 0.05; plot(1.5, mean(tRespID(:)), 'k*'); hold on; end
            
            yLim = quantile(tRespID(:), [0.05 0.95]);
            ax = gca;
            ax.XTick = [0.5 2.5];
%             ax.YLim = [-0.05 yLim(2)];
            xlim([0 3])

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



























