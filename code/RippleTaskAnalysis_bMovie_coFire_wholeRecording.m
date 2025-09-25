

close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/eeglab2022.0/functions/'))

% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%
subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P48CS_R1', 'P48CS_R2','P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};
% epilFocus = {{' '}, {' '}, ... %P41
%              {' '}, {' '}, ... %P42
%              {'LHIP', 'LAMY'}, {'LHIP', 'LAMY'}, ... %P43
%              {'RHIP', 'RAMY'}, ... %P44
%              {'RHIP', 'RAMY'}, {'RHIP', 'RAMY'}, ... %P47
%              {'LHIP', 'LAMY'}, {'LHIP', 'LAMY'}, ... %P48
%              {'LHIP', 'LAMY'}, {'LHIP', 'LAMY'}, ... %P49
%              {' '}, {' '}, ... %P51
%              {'LHIP', 'LAMY', 'RHIP','RAMY'}, {'LHIP', 'LAMY',  'RHIP','RAMY'}, ... %P53
%              {'RHIP', 'RAMY'}, {'RHIP', 'RAMY'}, ... %P54
%              {'RHIP', 'RAMY'}, {'RHIP', 'RAMY'}, ... %P55
%              {'LHIP', 'LAMY', 'RHIP','RAMY'}, {'LHIP', 'LAMY',  'RHIP','RAMY'}, ... %P56
%              {' '}, {' '}, ... %P57
%              {' '}, ... %P58
%              {'LHIP', 'LAMY'}, ... %P60
%              {'RHIP', 'RAMY'}, {'RHIP', 'RAMY'} ... %P62
% 
%              };

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
dataDirectory =     '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/preprocess';
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
taskfiles = dir(fullfile(dataDirectory, '*task*'));
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
    testRip = zeros(2.5e3, length(regions), length(regions), 2);
    
    trialsProbe = zeros(length(regions),length(regions),1e3);

    for iRa = 1:length(regions)
        for iRb = iRa:length(regions)
            whole_trials_cross_region_all{iRa, iRb} = nan(6.5e3, winLength);       
            testCoFire{iRa, iRb} = zeros(1, ceil(winLength/coFireSubsample));       
            load_trials_all_unitPair{iRa, iRb} = nan(1, 3.5e5);    
            unitPairID{iRa, iRb} = nan(4.5e5, 3);

            for ii = 1:9
    %             whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e5, ceil(winLength/coFireSubsample));
                whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e4, 2);
                whole_trials_cross_region_dur{iRa, iRb,ii} = nan(4.5e4, 2);
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

        modifier = 'IISremv_singleWire_v3_z25';
        tag = [recordingState,'_',location,'_',modifier];
%         filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
        filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
        microObj = load(fullfile(matExportFolder, filename));
        rippleStats = microObj.rippleStats;
        LFPchanNum = 1:length(rippleStats.locs);
        
        
        f = contains(LFPfilesMicr, subject);
        micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
        LFPtimeMicr = micrObj.times;
        chan_labels = rippleStats.chanLabels; locations =chan_labels;
        
        f = contains(LFPfilesMacr, subject);
        macrObj = load(fullfile(dataDirectory, LFPfilesMacr{f}));
        LFPtimeMacr = macrObj.times;


        f = contains(taskfiles, subject);
        trials = load(fullfile(dataDirectory, taskfiles{f}));

        if length(rippleStats.locs) ~= length(locations); error('error formatting channel labels.'); end

        rippleStats.recordingType = [repmat({'micro'}, [1 length(rippleStats.locs)])];
        rippleStats.locs = [cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false)];
        rippleStats.window(cellfun(@(X) isempty(X), rippleStats.window)) = {[1 1]};
        microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
        rippleStats.window = [cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', microObj.rippleStats.window, 'UniformOutput', false)];

%         rippleStats.chanLabels = [rippleStats.chanLabels; microObj.rippleStats.chanLabels];
        rippleStats.density = [ microObj.rippleStats.density];
        rippleStats.microTimes = LFPtimeMicr;
    %     emptyNum = [emptyNum cellfun(@(X) (length(X)/(rippleStats.recordingLength/1e3/60)) < 5 , rippleStats.locs)];
        rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
        for chRipp = 1:size(rippMask,1) %sum(strcmp(rippleStats.recordingType, 'macro'))+1:size(rippMask,1) 
            if length(rippleStats.window{chRipp}) > 2
    %             if strcmp(rippleStats.recordingType{chRipp}, 'macro'); times = rippleStats.macroTimes;
                if strcmp(rippleStats.recordingType{chRipp}, 'macro'); continue;
                elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); times = rippleStats.microTimes; 
    %             elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); continue; 
                end

                if contains(rippleStats.chanLabels{chRipp}, 'SPE')
                    rippMask(chRipp,:) = nan(1,length(rippMask));
                    continue
                end

                if strcmp(rippleStats.recordingType{chRipp}, 'micro')
                    chLoc = regexprep(rippleStats.chanLabels{chRipp}, '\d', '');
                    chNum = regexp(rippleStats.chanLabels{chRipp}, '\d+', 'match');
                    chNum = str2double(chNum);
                    epCh = find(contains(bundleLab, chLoc));

                    bch = badChan(subj, epCh);
                    bch = str2double(split(bch,','));

                    if contains(bundleNotes(subj, epCh), 'Bad') || any(ismember(bch, chNum))
                        rippMask(chRipp,:) = nan(1,length(rippMask));
                        continue
                    end
                end




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

        rippMask(:, end:end+4e3) = nan; %pad the ending
        micrObj.lfp_data(:, end:end+4e3) = nan; %pad the ending


        rippAll = sum(rippMask);


        units = LoadSpikeTimes(subject,'RutishauserLab', 'bmovie');
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
            FR = sum(N) / rippleStats.recordingLength * 1000;
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
        
        taskMask = false(2, length(rippMask));
        taskMaskConf = false(2, length(rippMask));
        taskMaskUnsure = false(2, length(rippMask));
        for iT = 1:length(trials.start_time)-1

            trialConf = trials.response_confidence(iT+1);
            probeTime =  round(trials.start_time(iT+1)*1e3);
            respTime  =  round(trials.response_time(iT+1)*1e3);
            
            taskMask(1, probeTime-1000:probeTime-1) = true; %baseline 
            taskMask(2, probeTime+1:probeTime+2.5e3)  = true; %probe 
            
            if trialConf == 3
                taskMaskConf(1, probeTime-1000:probeTime-1) = true; %baseline
                taskMaskConf(2, probeTime+1:probeTime+2.5e3)    = true; %probe 
            else
                taskMaskUnsure(1, probeTime-1000:probeTime-1) = true; %baseline
                taskMaskUnsure(2, probeTime+1:probeTime+2.5e3)    = true; %probe 
            end
            
            
            
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


                uLFPregionA = uChan(contains(uLocations, regions(iRa)));
                uLFPregionB = uChan(contains(uLocations, regions(iRb)));
                uRegionA = find(contains(uLocations, regions(iRa)));
                uRegionB = find(contains(uLocations, regions(iRb)));
                if isempty(uRegionA) || isempty(uRegionA); continue; end

               

    % 
    %                 datArip =  ripTrialFullChan(ovA(~isnan(ovA)),:);
    %                 datBrip =  ripTrialFullChan(ovB(~isnan(ovB)),:);

                
    %                 
                datArip = sum(rippMask(contains(locations, regions(iRa)),:), 'omitnan') > 0;
                datBrip = sum(rippMask(contains(locations, regions(iRb)),:), 'omitnan') > 0;
                if sum(datArip) == 0 || sum(datBrip) == 0; continue; end
                
                datABcorip     = false(2, length(rippMask)); datABnorip       = false(2, length(rippMask));
                datABcoripConf = false(2, length(rippMask)); datABcoripUnsure = false(2, length(rippMask));
                datABnoripConf = false(2, length(rippMask)); datABnoripUnsure = false(2, length(rippMask));
                for ii = 1:2
                    datABcorip(ii,:) = datArip > 0 & datBrip > 0 & taskMask(ii,:);
                    datABnorip(ii,:) = datArip == 0 & datBrip == 0  & taskMask(ii,:);
                    datABcoripConf(ii,:) = datABcorip(ii,:) & taskMaskConf(ii,:);
                    datABcoripUnsure(ii,:) = datABcorip(ii,:) & taskMaskUnsure(ii,:);
                    datABnoripConf(ii,:) = datABnorip(ii,:) & taskMaskConf(ii,:);
                    datABnoripUnsure(ii,:) = datABnorip(ii,:) & taskMaskUnsure(ii,:);


                
                end                

%                 temp = datABcoripConf(2,taskMaskConf(2,:));
%                 temp = reshape(temp, [], 2.5e3);
%                 testRip(:, iRa, iRb, 1) = testRip(:, iRa, iRb, 1) + mean(temp)';
% 
%                 temp = datABcoripUnsure(2,taskMaskUnsure(2,:));
%                 temp = reshape(temp, [], 2.5e3);
%                 testRip(:, iRa, iRb, 2) = testRip(:, iRa, iRb, 2) + mean(temp)';

                nUnits = length(uLFPregionA) + length(uLFPregionB);

                for uA = 1:length(uRegionA)
                    for uB = 1:length(uRegionB)
                        
%                         datArip =  rippMask(LFPchanNum == uChan(uRegionA(uA)),:);
%                         datBrip =  rippMask(LFPchanNum == uChan(uRegionB(uB)),:);

                        for ii = 1:2
%                             datABcorip = datArip > 0 & datBrip > 0 & taskMask(ii,:);
%                             datABnorip = datArip == 0 & datBrip == 0  & taskMask(ii,:);
                            if isempty(datArip) || isempty(datBrip); continue; end

                            
                            %co-ripples
                            datA =  find(spikeMask(uRegionA(uA),datABcorip(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABcorip(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,1}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,1}(cUnitPair(iRa, iRb), ii) = sum(datABcorip(ii,:));
                            
                            datA =  find(spikeMask(uRegionA(uA),datABcoripConf(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABcoripConf(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,4}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,4}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence == 3);

                            datA =  find(spikeMask(uRegionA(uA),datABcoripUnsure(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABcoripUnsure(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,5}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,5}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence < 3);

                            
        
                            %no-ripples
                            datA =  find(spikeMask(uRegionA(uA),datABnorip(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABnorip(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,2}(cUnitPair(iRa, iRb), ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,2}((cUnitPair(iRa, iRb)), ii) = sum(datABnorip(ii,:));
                            
                            datA =  find(spikeMask(uRegionA(uA),datABnoripConf(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABnoripConf(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,6}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,6}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence == 3);

                            datA =  find(spikeMask(uRegionA(uA),datABnoripUnsure(ii,:))); datB =  find(spikeMask(uRegionB(uB),datABnoripUnsure(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,7}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,7}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence < 3);
                            
                           

                            %all co-firing
                            datA =  find(spikeMask(uRegionA(uA),taskMask(ii,:))); datB =  find(spikeMask(uRegionB(uB),taskMask(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
%                             whole_trials_unit_cross_region_all{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(datAB);
                            whole_trials_cross_region_dur{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(taskMask(ii,:));
                            
                            datABconf = taskMaskConf(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABconf)); datB =  find(spikeMask(uRegionB(uB),datABconf));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,8}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,8}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence == 3);

                            datABunsure = taskMaskUnsure(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABunsure)); datB =  find(spikeMask(uRegionB(uB),datABunsure));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,9}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,9}(cUnitPair(iRa, iRb),ii) = sum(trials.response_confidence < 3);

                        end
                        
                        unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),1) = subj;
                        unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),2:3) = [uA uB];
                        cUnitPair(iRa, iRb) = cUnitPair(iRa, iRb) + 1;
                    end
                end
            end
        end


    end






%     fprintf('number of total trials: %i, number of subject trials: %i \n\n', cTrial - 1, length(trials.start_time))


    %     save(fullfile(exportDir, 'ripDetections', sprintf('%s_rippleStats.mat', subject)), 'rippleStats' ,'-v7.3')

    
    correct_trials_all(cTrial:end) = [];
    load_trials_all(cTrial:end) =[];

    probe_in_out(cTrial:end) = [];
    subjID(cTrial:end) = [];

    fprintf('done processing .. \n')
    toc
%     save(sprintf('Sternbger_ripple_results_z25_%icoF_splitHem_template_v2.mat', coFireWindow*2), '-v7.3')
end    
%%

figure;
for iRa = 1:length(regionsPlot)
    for iRb = iRa+1:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        all = zeros(1,2e3);
        
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                
                
                all =  all + smoothdata(testRip(:,iiA,iiB,1)', 'gaussian', 500);

                
                
            end
        end
        
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        subplot(5,5,c)
        plot(all)
    end
end

%%


close all
regionColors =  brewermap(12, 'Dark2');
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

thresh = 0;
coRegionColors = nan(10,3);
coRall = [];
noRall = [];
baseAll = [];
figure('Position',[6 975 510 598])
for iRa = 1:length(regionsPlot)
    for iRb = iRa:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        coR = [];
        noR = [];
        all = [];
        base = [];
        
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                
                if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                
                tempBase= [];
                tempCoR= [];
                tempNoR= [];
                for ii = 1:2
                    if ii == 1
                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,ii) * 1e3;                
    %                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,ii) * 1e3;                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,ii) * 1e3;
                    else
                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,ii) * 1e3;                
%                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,1}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,1}(:,ii) * 1e3;                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,2}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,2}(:,ii) * 1e3;
                    end
                end 
                
%                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];
                
%                 all = [all tempAll];
                coR = [coR tempCoR];
                noR = [noR tempNoR];
                base = [base tempBase];
                
                
            end
        end
        
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        subplot(5,5,c)
        
        keep = (coR(2,:) > 0) & ...
               (noR(2,:) > 0) & ...
               (base(2,:) > 0);      
        if sum(keep) == 0; continue; end
        
        coRall = [coRall coR(:, keep)];
        noRall = [noRall noR(:, keep)];
        baseAll = [baseAll base(:, keep)];
        xlim([0.5 4.5])
        hl = hline(median(base(1,keep),'omitnan'), 'k-');
        hl.LineWidth = 0.5;
  
        errobarPatch(1, base(1,keep)', coRegionColors(c,:)*0.6, 0.2);

        
        errobarPatch(1:2, base(1:end,keep)', coRegionColors(c,:)*0.6, 0.1);
        errobarPatch([1, 2.2], coR(1:end,keep)', coRegionColors(c,:), 0.1);
        errobarPatch([1, 1.8], noR(1:end,keep)', 'k', 0.1);
%         errobarPatch(1, coR', coRegionColors(c,:), 0.3);
%         errobarPatch(2, noR', coRegionColors(c,:)*0.3, 0.3);
%         errobarPatch(3, all', coRegionColors(c,:)*0.3, 0.3);
%         errobarPatch(4, base', coRegionColors(c,:)*0.3, 0.3);
%         group_num = [ones(1, length(coR)), 2*ones(1, length(noR)), 3*ones(1, length(all))];
%         dat = [coR, noR, all]';
%         h = daboxplot(dat,'groups',group_num, ...
%         'mean',1,'outlier', 0, 'color',regionColors,'xtlabels',{'CoR', 'NoR', 'all'}); hold on;
%         [P,H] = signrank(coR, all);

%         title(sprintf('%i %e', sum(~isnan(coR(1,:))), P))
        xlim([0.5 2.5])
        ax = gca;
        ax.XTick = [1:2];
        ax.XTickLabel = {'B','P'};
        title(num2str(sum(keep)))


        
        
        
        c = c+1;

    end
end
fig = gcf;
fig.Color = 'w';


plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
thresh = 0; t = 2;
coRegionColors = nan(10,3);

figure('Position',[6 975 881 598])
% figure('Position',[6 975 881 598])
% figure('Position',[6 975 881 598])
Pval = nan(length(regionsPlot), length(regionsPlot), 3);
for iRa = 1:length(regionsPlot)
    for iRb = iRa:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        coR = [];
        noR = [];
        all = [];
        base = [];
        
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                
                if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                
                tempBase= [];
                tempCoR= [];
                tempNoR= [];
                for ii = 1:2
                    
                    tempBase(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,7+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,7+ii}(:,t);
%                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                    tempCoR(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,3+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,3+ii}(:,t);                
                    tempNoR(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,5+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,5+ii}(:,t);


                end 
                
%                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];
                
%                 all = [all tempAll];
                coR = [coR tempCoR];
                noR = [noR tempNoR];
                base = [base tempBase];
                
                
            end
        end
        
        keep = coR(1,:) > thresh | coR(2,:) > thresh;
        
        if sum(keep) < 2; continue; end
        
%         base = base(:,keep);
%         base = zscore(base(:));
%         base = reshape(base, 2, []);
        
%         coR = coR(:,keep);
%         coR = zscore(coR(:));
%         coR = reshape(coR, 2, []);
        
%         noR = noR(:,keep);
%         noR = zscore(noR(:));
%         noR = reshape(noR, 2, []);
%         [stats, df, P, surrog] = statcond({base(1,keep), base(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
%         Pval(iRa,iRb,1) = P;
%         [stats, df, P, surrog] = statcond({noR(1,keep), noR(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
%         Pval(iRa,iRb,2) = P;
%         [stats, df, P, surrog] = statcond({coR(1,keep), coR(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
%         Pval(iRa,iRb,3) = P;
        
        [H,P,CI,STATS] = ttest(base(1,:), base(2,:),'tail', 'right');
        Pval(iRa,iRb,1) = P;
        [H,P,CI,STATS] = ttest(noR(1,:), noR(2,:),'tail', 'right');
        Pval(iRa,iRb,2) = P;
        [H,P,CI,STATS] = ttest(coR(1,:), coR(2,:), 'tail', 'right');
        Pval(iRa,iRb,3) = P;
        
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        
        figure(2)
        subplot(5,5,c)
        keep = (coR(1,:) > 0/70 | coR(2,:) > 0/70) & ...
               (noR(1,:) > 0/70 | noR(2,:) > 0/70) & ...
               (base(1,:) > 0/70 | base(2,:) > 0/70);
        baseModul = [base(1,keep)-base(2,keep)]./mean([base(1,keep) base(2,keep)]) * 100;
        coRmodul = [coR(1,keep)-coR(2,keep)]./mean([coR(1,keep) coR(2,keep)]) * 100;
        noRmodul = [noR(1,keep)-noR(2,keep)]./mean([noR(1,keep) noR(2,keep)]) * 100;
        modul = [baseModul; coRmodul; noRmodul];
        groups = [ones(1, sum(keep)) 2*ones(1, sum(keep)) 3*ones(1, sum(keep))];
        
%         vp = violinplot([baseModul, coRmodul, noRmodul], groups); 
%         for iV = 1:length(vp)
%             vp(iV).ViolinPlot.Visible = 'off';
%             vp(iV).ScatterPlot.Visible = 'off';
%             vp(iV).WhiskerPlot.Visible = 'off';
%             vp(iV).BoxWidth = 0.1;
%             
%         end
        YY = mean(modul, 2);
         % Compute MAD
        mad_Y = median(abs(modul' - median(modul,2)'));

        % Compute robust standard error
        YYdev = mad_Y / sqrt(sum(keep));
        
        YYdev = std(modul,[], 2) / sqrt(sum(keep));
        xx = 1:3;
        w = 0.2;
        for ii = 1:length(xx)

            p = patch([xx(ii)-w xx(ii)+w xx(ii)+w xx(ii)-w], [YY(ii)-YYdev(ii) YY(ii)-YYdev(ii) YY(ii)+YYdev(ii) YY(ii)+YYdev(ii)], c); hold on;
            p.FaceAlpha = 0.4;
            p.FaceColor = coRegionColors(c,:);
            p.EdgeAlpha = 0;

            pl = plot([xx(ii)-w xx(ii)+w], [YY(ii) YY(ii)], '-'); hold on;
            pl.LineWidth = 2;
            pl.Color = coRegionColors(c,:);

        end
%         errobarPatch([1:2], coR(1:end,keep)', coRegionColors(c,:), 0.1);

        xlim([0.5 3.5])
%         ax = gca;
%         ax.XTick = [1:2];
%         ax.XTickLabel = {'load3', 'load1'};
        fig = gcf;
        fig.Color = 'w';
        title(num2str(sum(keep)))
        hline(0)

        

%         errobarPatch([1:2], coR(1:end,:)', coRegionColors(c,:), 0.1);
% 
%         xlim([-0.5 2.5])
%         ax = gca;
%         ax.XTick = [1:2];
%         ax.XTickLabel = {'conf', 'unsure'};
%         fig = gcf;
%         fig.Color = 'w';
%         title(sprintf('%e', Pval(iRa,iRb,3)))
%         
%         figure(3)
%         subplot(5,5,c)
% 
%         
%         errobarPatch([1:2], noR(1:end,:)', coRegionColors(c,:), 0.1);
% 
%         xlim([-0.5 2.5])
%         ax = gca;
%         ax.XTick = [1:2];
%         ax.XTickLabel = {'conf', 'unsure'};
%         fig = gcf;
%         fig.Color = 'w';
%         title(sprintf('%e', Pval(iRa,iRb,2)))
% 
%         figure(4)
%         subplot(5,5,c)
% 
%         
%         errobarPatch([1:2], base(1:end,:)', coRegionColors(c,:), 0.1);
% 
%         xlim([-0.5 2.5])
%         ax = gca;
%         ax.XTick = [1:2];
%         ax.XTickLabel = {'conf', 'unsure'};
%         fig = gcf;
%         fig.Color = 'w';
%         title(sprintf('%e', Pval(iRa,iRb,1)))

        
        

        
        
        
        c = c+1;

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

baseModul =[]
noRmodul = [];
noRmodul = [];
coR = [];
noR = [];
all = [];
base = [];
for iRa = 1:length(regionsPlot)
    for iRb = iRa:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                
                if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                
                tempBase= [];
                tempCoR= [];
                tempNoR= [];
                for ii = 1:2
                    
                    tempBase(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,7+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,7+ii}(:,t);
%                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                    tempCoR(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,3+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,3+ii}(:,t);                
                    tempNoR(ii,:) = [whole_trials_unit_cross_region_all{iiA, iiB,5+ii}(:,t)] ./ whole_trials_cross_region_dur{iiA, iiB,5+ii}(:,t);

                end 
                
%                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];
                 keep = (tempCoR(1,:) > 0/70 | tempCoR(2,:) > 0/70) ...
                      & (tempNoR(1,:) > 0/70 | tempNoR(2,:) > 0/70) & ...
                        (tempBase(1,:) > 0/70 | tempBase(2,:) > 0/70);

%                 all = [all tempAll];
                coR = [coR tempCoR(:,keep)];
                noR = [noR tempNoR(:,keep)];
                base = [base tempBase(:,keep)];
                
                
            end
        end
        
        baseModul = [baseModul [base(1,:)-base(2,:)]./mean([base(1,:) base(2,:)]) * 100];
        coRmodul = [coRmodul [coR(1,:)-coR(2,:)]./mean([coR(1,:) coR(2,:)]) * 100];
        noRmodul = [noRmodul [noR(1,:)-noR(2,:)]./mean([noR(1,:) noR(2,:)]) * 100];
%         modul = [baseModul; coRmodul; noRmodul];
        

        
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        
       
        
        

        
        
        
        c = c+1;

    end
end
figure('Position',[1224 1354 560 219])

% keep = base(1,:) > thresh & base(2,:) > thresh;

base = baseModul;
% base = zscore(base(:));
% base = reshape(base, 2, []);

coR = coRmodul;
% coR = zscore(coR(:));
% coR = reshape(coR, 2, []);

noR = noRmodul;
% noR = zscore(noR(:));
% noR = reshape(noR, 2, []);
Pval = [];   
% [stats, df, P, surrog] = statcond({base, zeros(1, length(base))},'paired','on','tail', 'both', 'method', 'perm', 'naccu', 1000, 'verbose','off');
% Pval(1) = P;
% [stats, df, P, surrog] = statcond({noR, zeros(1, length(noR))},'paired','on','tail', 'both', 'method', 'perm', 'naccu', 1000, 'verbose','off');
% Pval(2) = P;
% [stats, df, P, surrog] = statcond({coR, zeros(1, length(coR))},'paired','on','tail', 'both', 'method', 'perm', 'naccu', 100, 'verbose','off');
% Pval(3) = P;

[H,P,CI,STATS] = ttest(base, zeros(1, length(base)),'tail', 'both');
Pval(1) = P;
[H,P,CI,STATS] = ttest(noR, zeros(1, length(noR)),'tail', 'both');
Pval(2) = P;
[H,P,CI,STATS] = ttest(coR, zeros(1, length(coR)), 'tail', 'both');
Pval(3) = P;

figure('Position',[1219 1160 218 413])
errobarPatch([1], coR', coRegionColors(c-10,:), 0.2, 'DevType', 'SEM')

xlim([0.5 2.5])
ax = gca;
ax.XTick = [1:2];
ax.XTickLabel = {'conf', 'unsure'};
fig = gcf;
fig.Color = 'w';
title(sprintf('%e', Pval(3)))
% subplot(1,3,2)
errobarPatch([2], noR', coRegionColors(c-10,:), 0.2, 'DevType', 'SEM');

xlim([0.5 2.5])
ax = gca;
ax.XTick = [1:2];
% ax.XTickLabel = {'conf', 'unsure'};
fig = gcf;
fig.Color = 'w';
title(sprintf('%e', Pval(2)))
errobarPatch([3], base', coRegionColors(c-10,:), 0.2, 'DevType', 'SEM');

xlim([0.5 3.5])
ax = gca;
ax.XTick = [1:3];
ax.XTickLabel = {'conf', 'unsure'};
fig = gcf;
fig.Color = 'w';
title(sprintf('%e', Pval(1)))
