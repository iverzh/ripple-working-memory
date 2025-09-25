

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
    UnitPairType = cell(length(regions), length(regions));
    UnitPairFR = cell(length(regions), length(regions));

    trialsProbe = zeros(length(regions),length(regions),1e3);

    for iRa = 1:length(regions)
        for iRb = iRa:length(regions)
            whole_trials_cross_region_all{iRa, iRb} = nan(6.5e3, winLength);       
            testCoFire{iRa, iRb} = zeros(1, ceil(winLength/coFireSubsample));       
            load_trials_all_unitPair{iRa, iRb} = nan(1, 3.5e5);    
            unitPairID{iRa, iRb} = nan(4.5e5, 3);
            UnitPairType{iRa, iRb} = cell(2, 4.5e5);
            UnitPairFR{iRa, iRb} = nan(2, 4.5e5);


            for ii = 1:9
    %             whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e5, ceil(winLength/coFireSubsample));
                whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(4.5e4, 4);
                whole_trials_cross_region_dur{iRa, iRb,ii} = nan(4.5e4, 4);
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
        FR = [];
        for iU = 1:nNeuron 
            X = units{unitsAll(iU),2};
        %     ISIx = [0 diff(X')];
        %     xx = randperm(length(ISIx), length(ISIx));
        %     Xshuff =[X(1)];
        %     for ii = 2:length(ISIx); Xshuff(ii) = Xshuff(end) + ISIx(xx(ii)); end
            [N,EDGES] = histcounts(X,binEdges);
            FR(iU) = sum(N) / rippleStats.recordingLength * 1000;
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
        
        taskMask = false(4, length(rippMask));
        taskMaskLoad3 = false(4, length(rippMask));
        taskMaskLoad1 = false(4, length(rippMask));
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
            taskMask(3, maintenanceTime:probeTime)    = true; %maintenance 
            taskMask(4, probeTime+1:probeTime+1e3)    = true; %probe 
            
            if trialLoad == 3
                taskMaskLoad3(1, image1Time-1000:image1Time-1) = true; %baseline
                taskMaskLoad3(2, image1Time:image1Time+1e3)    = true; %img1
                taskMaskLoad3(3, maintenanceTime:probeTime)    = true; %maintenance 
                taskMaskLoad3(4, probeTime+1:probeTime+1e3)    = true; %probe 
            else
                taskMaskLoad1(1, image1Time-1000:image1Time-1) = true; %baseline
                taskMaskLoad1(2, image1Time:image1Time+1e3)    = true; %img1
                taskMaskLoad1(3, maintenanceTime:probeTime)    = true; %maintenance 
                taskMaskLoad1(4, probeTime+1:probeTime+1e3)    = true; %probe 
            end
            
            
            
        end
       

         %across region co-rippling
        for iRa = 1:length(regions)
            for iRb = iRa:length(regions)
                
%                 if iRa == iRb; error(); end
                
                
                uLFPregionA = uChan(contains(uLocations, regions(iRa)));
                uLFPregionB = uChan(contains(uLocations, regions(iRb)));
                uRegionA = find(contains(uLocations, regions(iRa)));
                uRegionB = find(contains(uLocations, regions(iRb)));
                

    %                 
                datArip = sum(rippMask(contains(locations, regions(iRa)),:), 'omitnan') > 0;
                datBrip = sum(rippMask(contains(locations, regions(iRb)),:), 'omitnan') > 0;
                
                if sum(datArip) == 0 || sum(datBrip) == 0; continue; end


                nUnits = length(uLFPregionA) + length(uLFPregionB);

                for uA = 1:length(uRegionA)
                    for uB = 1:length(uRegionB)
                        typeA = units{uRegionA(uA),3};
                        typeB = units{uRegionB(uB),3};
%                         datArip =  rippMask(LFPchanNum == uChan(uRegionA(uA)),:);
%                         datBrip =  rippMask(LFPchanNum == uChan(uRegionB(uB)),:);
                        if uLFPregionA(uA) == uLFPregionB(uB); continue; end
                        
%                         datArip = rippMask(uLFPregionA(uA) == LFPchanNum,:); datArip(isnan(datArip)) = 0;
%                         datBrip = rippMask(uLFPregionB(uB) == LFPchanNum,:); datBrip(isnan(datBrip)) = 0;
                        
                        
%                         if sum(datArip, 'all') == 0 || sum(datBrip, 'all') == 0; continue; end

                        
                        for ii = 1:4
                            datABcorip = datArip > 0 & datBrip > 0 & taskMask(ii,:);
                            datABnorip = datArip == 0 & datBrip == 0  & taskMask(ii,:);

                            
                            %co-ripples
                            datA =  find(spikeMask(uRegionA(uA),datABcorip)); datB =  find(spikeMask(uRegionB(uB),datABcorip));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,1}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,1}(cUnitPair(iRa, iRb), ii) = sum(datABcorip);
                            
                            datABcoripLoad3 = datABcorip & taskMaskLoad3(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABcoripLoad3)); datB =  find(spikeMask(uRegionB(uB),datABcoripLoad3));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,4}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,4}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 3);

                            datABcoripLoad1 = datABcorip & taskMaskLoad1(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABcoripLoad1)); datB =  find(spikeMask(uRegionB(uB),datABcoripLoad1));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,5}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,5}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 1);

                            
        
                            %no-ripples
                            datA =  find(spikeMask(uRegionA(uA),datABnorip)); datB =  find(spikeMask(uRegionB(uB),datABnorip));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,2}(cUnitPair(iRa, iRb), ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,2}((cUnitPair(iRa, iRb)), ii) = sum(datABnorip);
                            
                            datABnoripLoad3 = datABnorip & taskMaskLoad3(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABnoripLoad3)); datB =  find(spikeMask(uRegionB(uB),datABnoripLoad3));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,6}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,6}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 3);

                            datABnoripLoad1 = datABnorip & taskMaskLoad1(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABnoripLoad1)); datB =  find(spikeMask(uRegionB(uB),datABnoripLoad1));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,7}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,7}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 1);
                            
                           

                            %all co-firing
                            datA =  find(spikeMask(uRegionA(uA),taskMask(ii,:))); datB =  find(spikeMask(uRegionB(uB),taskMask(ii,:)));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
%                             whole_trials_unit_cross_region_all{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(datAB);
                            whole_trials_cross_region_dur{iRa, iRb,3}(cUnitPair(iRa, iRb), ii) = sum(taskMask(ii,:));
                            
                            datABLoad3 = taskMaskLoad3(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABLoad3)); datB =  find(spikeMask(uRegionB(uB),datABLoad3));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,8}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,8}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 3);

                            datABLoad1 = taskMaskLoad1(ii,:);
                            datA =  find(spikeMask(uRegionA(uA),datABLoad1)); datB =  find(spikeMask(uRegionB(uB),datABLoad1));
                            s = datA - datB';
                            whole_trials_unit_cross_region_all{iRa, iRb,9}(cUnitPair(iRa, iRb),ii) = sum(s(:) >= -coFireWindow*2 & s(:) <= coFireWindow*2);
                            whole_trials_cross_region_dur{iRa, iRb,9}(cUnitPair(iRa, iRb),ii) = sum(trials.loads == 1);

                        end
                        
                        unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),1) = subj;
                        unitPairID{iRa, iRb}(cUnitPair(iRa, iRb),2:3) = [uA uB];
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
close all
regionColors =  brewermap(12, 'Dark2');
regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  

thresh = 0;
coRegionColors = nan(10,3);
coRall = [];
noRall = [];
baseAll = [];
Npair = nan(5,5);
pBaseline = nan(length(regionsPlot), length(regionsPlot),3, 3);
modBaseline = nan(length(regionsPlot), length(regionsPlot),3, 3);
keepAll = cell(length(regionsPlot), length(regionsPlot));
figure('Position',[6 896 771 677])
for iRa = 1:length(regionsPlot)
    for iRb = iRa:length(regionsPlot)
        
        iRapl = find(contains(regions, regionsPlot{iRa}));
        iRbpl = find(contains(regions, regionsPlot{iRb}));
        
        coR = [];
        noR = [];
        base = [];
        fireA = [];
        fireB = [];
        testAll = zeros(1,12905);
        for iiA = iRapl
            for iiB = iRbpl
                if iiA == iiB; continue; end
                if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                
                types = UnitPairType{iiA, iiB};
                
%                 iiType = strcmp(types(1,:), 'int') &  strcmp(types(2,:), 'int');
                iiType = strcmp(types(1,:), 'pyr') |  strcmp(types(1,:), 'int');
                
                tempBase= [];
                tempCoRcount= [];
                tempNoR= [];
                for ii = 1:4 
                    if ii == 1
                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(iiType,ii) * 1e3;                
    %                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoRcount(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(iiType,ii) * 1e3;                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(iiType,ii) * 1e3;
                    else
                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(iiType,ii) * 1e3;                
%                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoRcount(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,1}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,1}(iiType,ii) * 1e3;                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,2}(iiType,ii) ./ whole_trials_cross_region_dur{iiA, iiB,2}(iiType,ii) * 1e3;
                    end
                end 
                
%                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
%                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];
                
%                 all = [all tempAll];
                coR = [coR tempCoRcount];
                noR = [noR tempNoR];
                base = [base tempBase];
                fireA = [fireA UnitPairFR{iiA,iiB}(1,iiType)];
                fireB = [fireB UnitPairFR{iiA,iiB}(2,iiType)];
                
            end
        end
        
        c = (iRb - 1) * 5 + iRa;
        coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        subplot(5,5,c)
        
        keep = coR(2,:) > thresh | coR(3,:) > thresh | coR(4,:) > thresh;
%         keep = base(2,:) > thresh | base(3,:) > thresh | base(4,:) > thresh;

        keepAll{iRa,iRb} = keep;
%         keep = fireA > 1 & fireB > 1;
%         keep = base(2,:) > thresh | base(3,:) > thresh | base(4,:) > thresh;
        Npair(iRa, iRb) = sum(keep);
%         keep = ~isnan(coR(2,:)) | ~isnan(coR(3,:)) | ~isnan(coR(4,:) > thresh);
%         keep = base(2,:) > thresh | base(3,:) > thresh | base(4,:) > thresh;
%         keep = coR(4,:) > thresh;

        for iT = 1:3
%             pBaseline(iRa, iRb, iT, 1) = signrank(base(1, keep), coR(iT+1, keep), 'tail', 'both');
%             pBaseline(iRa, iRb, iT, 3) = signrank(base(1, keep), noR(iT+1, keep), 'tail', 'both');
%             pBaseline(iRa, iRb, iT, 2) = signrank(base(1, keep), base(iT+1, keep), 'tail', 'both');
            [pBaseline(iRa, iRb, iT, 1), perm_stats, obs_stat] = pairedPermutationTest(base(1, keep)', coR(iT+1, keep)', 10000);
            [pBaseline(iRa, iRb, iT, 3), perm_stats, obs_stat] = pairedPermutationTest(base(1, keep)', noR(iT+1, keep)', 10000);
            [pBaseline(iRa, iRb, iT, 2), perm_stats, obs_stat] = pairedPermutationTest(base(1, keep)', base(iT+1, keep)', 10000);
            
            modBaseline(iRa, iRb, iT, 1) = (median(coR(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
            modBaseline(iRa, iRb, iT, 3) = (median(noR(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
            modBaseline(iRa, iRb, iT, 2) = (median(base(iT+1, keep)) - median(base(1, keep))) / median(base(1, keep));
        end

        coRall = [coRall coR(:, keep)];
        noRall = [noRall noR(:, keep)];
        baseAll = [baseAll base(:, keep)];
        xlim([0.5 4.5])
%         hl = hline(median(base(1,keep),'omitnan'), 'k-');
        hl = hline(median(base(1,keep),'omitnan'), 'k-');
        hl.LineWidth = 0.5;
  
        errobarPatch(1, base(1,keep)', 'k', 0.2);

        
        errobarPatch([2:4]-0.2, base(2:end,keep)', coRegionColors(c,:)*0.6, 0.1);
        errobarPatch([2:4]+0.2-0.2, coR(2:end,keep)', coRegionColors(c,:), 0.1);
        errobarPatch([2:4]-0.2-0.2, noR(2:end,keep)', coRegionColors(c,:)*0.3, 0.1);
%         errobarPatch(1, coR', coRegionColors(c,:), 0.3);
%         errobarPatch(2, noR', coRegionColors(c,:)*0.3, 0.3);
%         errobarPatch(3, all', coRegionColors(c,:)*0.3, 0.3);
%         errobarPatch(4, base', coRegionColors(c,:)*0.3, 0.3);
%         group_num = [ones(1, length(coR)), 2*ones(1, length(noR)), 3*ones(1, length(all))];
%         dat = [coR, noR, all]';
%         h = daboxplot(dat,'groups',group_num, ...
%         'mean',1,'outlier', 0, 'color',regionColors,'xtlabels',{'CoR', 'NoR', 'all'}); hold on;
%         [P,H] = signrank(coR, all);

%         title(sprintf('%i', sum(keep)))
        xlim([0.5 4.7])
%         ylim([0 0.7])
        ax = gca;
        ax.XTick = [1:4];
        ax.XTickLabel = {'B','E1','M','P'};
        
        if iRa == 1; ylabel('co-fire rate [Hz]'); end
        
%         title(num2str(Npair(iRa,iRb)))
        
        pct = [median(coR(3,keep), 'omitnan') - median(base(1,keep), 'omitnan')]/ median(base(1,keep), 'omitnan');
        fprintf('%s <--> %s ripples maintenance %.2f\n', regionsPlot{iRa}, regionsPlot{iRb}, pct)
        pBaseline(iRa, iRb, 2, 1)
        
        axRng = range(ax.YLim);
        ax.YLim(1) = ax.YLim(1) - 0.2*axRng;
        
        c = c+1;

    end
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegions_%s.pdf', tag)))

[p,tbl,stats] = kruskalwallis([baseAll(4,:); coRall(3,:); noRall(4,:)]') 

%%
adj_p = nan(size(pBaseline));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pBaseline))]=fdr_bh(pBaseline(~isnan(pBaseline)),0.05,'pdep','yes');

% for ii = 1:3
%     figure('Position', [1219 1378 1004 195]);
%     for jj = 1:3
%         subplot(1,3,jj)
%         imagesc(adj_p(:,:,jj,ii)', [min(adj_p(:))-1e-2 1]); hold on;
%         colorbar;
%         hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
%         vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
%         cmap = flipud(slanCM('magma'));
%         cmap = [1 1 1; cmap];
%         colormap(cmap)
%         [yy xx] = find(adj_p(:,:,jj,ii)' < 0.05);
%         plot(xx, yy, 'r*'); hold on;
% 
%         
%     end
% end
% fig = gcf;
% fig.Color = 'w';

figure('Position', [1219 1378 600 600]);
ha = tight_subplot(3,3,[0.05 ,0.05],[.05 .05],[.07 .07]);
for ii = 1:3
    for jj = 1:3
        adj_pFL = fliplr(adj_p(:,:,jj,ii));
        modBaselineFL = fliplr(modBaseline(:,:,jj,ii));
        c = (ii - 1) * 3 + jj;
        axes(ha(c))
%         imagesc(modBaseline(:,:,jj,ii)', [-0.25 0.5]); hold on;
%         colorbar; caxis([-0.25 0.5]);
        
        pArray = [1 0.8 0.5 0.2 0.05 1e-3 1e-4 1e-5];
        pSize  = [ 2   3   6   10   15   20   25];
        
        % Interpolate between the colors
        cmap = slanCM('RdBu');
        cmapRGB = [cmap(1,:); [0.9 0.9 0.9]; cmap(end,:)];
        xq = linspace(0, 1, 256*3); % Query points for interpolation
        colorPositions = [0, 1/3, 1]; 
        cmap = interp1(colorPositions, cmapRGB, xq, 'linear');
        [xx, yy] = find(~isnan(adj_pFL));
        for p = 1:length(xx)
            pl = plot(xx(p), yy(p), 'rs'); hold on;
            
            pVal = adj_pFL(xx(p), yy(p)); 
            iiP = arrayfun(@(X) pArray(X) > pVal & pArray(X+1) <= pVal, 1:length(pArray)-1);
            if pVal == 1; iiP(1) = true;
            elseif pVal < 1e-5; iiP(end) = true; end
            pl.MarkerSize = pSize(iiP);
            
            pl.MarkerFaceColor = interp1(linspace(-0.25, 0.5, 256*3), cmap, modBaselineFL(xx(p), yy(p)), 'nearest','extrap');
            pl.MarkerEdgeColor = interp1(linspace(-0.25, 0.5, 256*3), cmap, modBaselineFL(xx(p), yy(p)), 'nearest','extrap');
            
            if xx(p) + yy(p) == 6
                pl = plot(xx(p), yy(p), 'rs'); hold on;
                pl.MarkerEdgeColor = 'g';
                pl.MarkerFaceColor = 'none';
                pl.MarkerSize = max(pSize)+3;
                pl.LineWidth = 1.5;
            end
            
        end
        
        
        % Define color positions (normalized between 0 and 1)

        % Define the RGB values of the three colors

        
        cmap = [1 1 1; cmap];
        
%         [xx, yy] = find(adj_pFL < 0.05 & abs(modBaselineFL) < 0.3);
%         pl = plot(xx, yy, 'k*'); hold on;
%         pl.MarkerSize = 4;
        [xx, yy] = find(adj_pFL < 0.05);
        if ~isempty(xx)
            pl = plot(xx, yy, '*'); hold on;
            pl.MarkerSize = 4;
            pl.Color = [0.8 0 0];
        end
        ylim([0.5 5.5])
        xlim([0.5 5.5])
%         hl = hline(-(1/2):size(adj_p,1)+(3/2), 'k-'); hold on;
%         vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
       box off
       ax = gca;
       ax.YTickLabel = fliplr(regionsPlot);
       ax.XTickLabel = regionsPlot;
       
       fig = gcf;
       fig.Color = 'w';

        
    end
end
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_%s.pdf', tag)))

figure('Position',[1214 1349 172 212]);
imagesc(modBaseline(:,:,jj,ii)', [-0.25 0.5]); hold on;
h = colorbar; 
caxis([-0.25 0.5]*100);
h.Ticks = [-0.25:0.25:0.5]*100;
colormap(cmap)
fig = gcf;
fig.Color = 'w';
h.Location = 'southoutside'; % Set horizontal colorbar

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_modlegend_%s.pdf', tag)))


figure('Position', [1219 1378-((ii-1)*200) 600 202]);
 for p = 1:length(pSize)
    pl = plot(1+p*(1.3/3), 1, 'rs'); hold on;
    pl.MarkerSize = pSize(p);
    pl.MarkerFaceColor = cmap(end,:);
    pl.MarkerEdgeColor = cmap(end,:);
            
 end
  axis off
 xlim([0 p])
 fig = gcf;
 fig.Color = 'w';
 savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsStats_plegend_%s.pdf', tag)))

 
 
cBox = [];
cBox(2,:) = [1 0 0];
cBox(3,:) = [1 0 0]*0.3;
cBox(1,:) = [1 0 0]*0.6;
figure('Position',[370 1318 362 209])
 group_num = [ones(1, length(coRall)), 2*ones(1, length(noRall)), 3*ones(1, length(baseAll))];
dat = [ coRall(2:end,:),  baseAll(2:end,:),noRall(2:end,:)]';
h = daboxplot(dat,'groups',group_num, 'mean',1,'outlier', 0, 'color',cBox, 'xshift', 1, 'boxalpha', 0.5); hold on;
hB = daboxplot( log10(baseAll(1,:)+1)', 'mean',1,'outlier', 0, 'color','k', 'boxalpha', 0.5, 'boxwidth', 0.6); hold on;

ax = gca;
ax.XTick = [1:4];
ax.XTickLabel = {'B','E1','M','P'};
ax.LineWidth = 1.0;
fig = gcf;
fig.Color = 'w';
xlim([0.5 4.5])
ylabel('co-fire rate [Hz]');
hl = hline(mean(baseAll(1,:)));
hl.LineStyle = '-';



%%

% close all
figure('Position',[682 1313 425 202])


dat = nan(size(baseAll,1), size(baseAll,2), 3);
pVal = nan(3,3);
effectsize = nan(3,3);
for t = 2:4
    [pVal(1,t-1), ~, effectsize(1,t-1)] = pairedPermutationTest(baseAll(1,:)', baseAll(t,:)', 2000);
    [pVal(2,t-1), ~, effectsize(2,t-1)] = pairedPermutationTest(baseAll(1,:)', coRall(t,:)', 2000);
    [pVal(3,t-1), ~, effectsize(3,t-1)] = pairedPermutationTest(baseAll(1,:)', noRall(t,:)',2000);
%     [pVal(1,t-1)] = pairedPermutationTest(baseAll(1,:)', baseAll(t,:)', 1000);
%     [pVal(2,t-1)] = pairedPermutationTest(baseAll(1,:)', coRall(t,:)', 1000);
%     [pVal(3,t-1)] = pairedPermutationTest(baseAll(1,:)', noRall(t,:)', 1000);
    bs = baseAll(1,:);
    xx = coRall(t,:);
    dat(t-1, :, 1) = (xx-bs)/ median(bs)*100;
    effectsize(1,t-1) = median(dat(t-1, :, 1));
    xx = baseAll(t,:);
    dat(t-1, :, 2) = (xx-bs)/ median(bs)*100;
    effectsize(2,t-1) = median(dat(t-1, :, 2));
    xx = noRall(t,:);
    dat(t-1, :, 3) = (xx-bs)/ median(bs)*100; 
    effectsize(3,t-1) = median(dat(t-1, :, 3));

end
% effectsize = effectsize/median(bs) * 100;
p_adj = nan(size(pVal));
[~, crit_p, adj_ci_cvrg, p_adj(~isnan(pVal))]=fdr_bh(pVal(~isnan(pVal)),0.05,'pdep','yes');


xlim([1.5 4.5])

hl = hline(0); hold on;
hl.LineStyle = '-';
hl.Color = 'k';
cBox = [      clr(8,:);
    mean([clr(1,:); clr(8,:)]);
              clr(1,:)];

% errobarPatch(hB.gpos(1), [baseAll(1,:)-median(baseAll(1,:))]'*100, 'k', 0.15); hold on;

for ii = 1:3
%     [bl, bf] = errobarPatch(h.gpos(ii, 1), dat(1,:,ii)', cBox(ii,:) , 0.10, 'DevType', 'IQR'); hold on;
    bf = bar(h.gpos(ii, 1), median(dat(1,:,ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii,:);
    
%     [bl, bf] = errobarPatch(h.gpos(ii, 2), dat(2,:,ii)', cBox(ii,:) , 0.10, 'DevType', 'IQR'); hold on;
    bf = bar(h.gpos(ii, 2), median(dat(2,:,ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii,:);

%     [bl, bf] =  errobarPatch(h.gpos(ii, 3), dat(3,:,ii)', cBox(ii,:) , 0.10, 'DevType', 'IQR'); hold on;
    bf = bar(h.gpos(ii, 3), median(dat(3,:,ii)), 0.20); hold on;
    bf.EdgeAlpha = 1;
    bf.EdgeColor = 'k';
    bf.LineWidth = 1;
    bf.FaceColor = cBox(ii,:);

    
    
    quantile(dat(1,:,ii), [0.25 0.75])
    medians = mean(dat(1,:,ii));
    mad_values = std(dat(1,:,ii), 1)/length(dat(1,:,ii));
%     b = bar(h.gpos(ii, 1), medians, 0.1); hold on;

    % Get x coordinates for error bars
%     x = b.XEndPoints;

    % Add error bars
%     errorbar(x, medians, mad_values, 'k', 'linestyle', 'none', 'LineWidth', 1.2);


    pVals = find(p_adj(ii,:) < 0.05);
    for iP = 1:length(pVals)
     
        iTm = pVals(iP);
        y_max = max([h.wh(iTm,ii,:).YData]);
        y_max = y_max + 0.10/1.27 * y_max;

%         text(mean(h.gpos(ii, iTm)), median(dat(iTm,:,ii)) + 2, '*', 'FontSize', 10, 'HorizontalAlignment', 'center'); hold on
    end
end


ax = gca;
ax.XTick = [1:4];
ax.XTickLabel = {'B','E1','M','P'};
ax.LineWidth = 1.0;
fig = gcf;
fig.Color = 'w';
xlim([1.5 4.5])
% vl = vline([2.5 3.5]);
% for v = 1:length(vl); vl(v).LineStyle = '--'; vl(v).LineWidth = 1; end
% ylim([-0.1 0.6]*100)

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskAcrossRegionsAll_%s.pdf', tag)))
%%
xx = median(baseAll(1,:));
yy = median(noRall(4,:));
yy = median(baseAll(4,:));

(yy -xx)/xx


%%
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));
fs = 1e3; %1e3, 400;
ctxParc = {'ACC', 'SMA', 'OFC' ,'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};
clc
close all
figure('Position', [1224 1316 328 257]); 
c = 1;
coRrate = [];
x = [];
xMu = [];
yCo = [];
yNo = [];
yCoMu = [];
yNoMu = [];
yCoSEM = [];
yNoSEM = [];

iRaAll = [];
iRbAll = [];
hemA =[];
hemB = [];
parcA = [];
parcB = [];
withinBundle = [];
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        parcelA = regions{iRa}(2:end);
        parcelB = regions{iRb}(2:end);

%         if ~strcmp(hemA,hemB); continue; end
%         if any(contains({parcelA, parcelB},{'AMY', 'HIP'})); continue; end


        
        types = UnitPairType{iRa, iRb};
        iiType = strcmp(types(1,:), 'pyr') | strcmp(types(1,:), 'int');
        fireA = UnitPairFR{iRa,iRb}(1,iiType);
        fireB = UnitPairFR{iRa,iRb}(2,iiType);
        tempCoRcount = [];
        tempCoRdur = [];
        tempNoRcount = [];
        tempNoRdur = [];
        tempBaseCount = [];
        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb,1}(iiType,:),2);
        tempCoRcount = [tempCoRcount temp'];
        temp = sum(whole_trials_cross_region_dur{iRa, iRb,1}(iiType,:),2);              
        tempCoRdur = [tempCoRdur temp'];

        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb,2}(iiType,:),2);
        tempNoRcount = [tempNoRcount temp'];
        temp = sum(whole_trials_cross_region_dur{iRa, iRb,2}(iiType,:),2);
        tempNoRdur = [tempNoRdur temp'];

        temp = sum(whole_trials_unit_cross_region_all{iRa, iRb,3}(iiType,:),2);
        tempBaseCount = [tempBaseCount temp'];

        
        tCo = ((tempCoRcount ./ tempCoRdur) * 1e3) ./ geomean([fireA; fireB]) ; 
        tNo = ((tempNoRcount ./ tempNoRdur) * 1e3) ./ geomean([fireA; fireB]);
        
        
                                    
        z = [(tCo(tempBaseCount >= 0)) (tNo(tempBaseCount >= 0))];
        tCo = z(1:length(z)/2);
        tNo = z((length(z)/2)+1:end);
        x = [x repmat(d, [1 length(tCo)])];

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
%         yCoMu = [yCoMu mean(z(1:length(tempCoRcount)))];
%         yNoMu = [yNoMu mean(z(length(tempCoRcount)+1:end))];
        

        
%         z = zscore([tCo(isfinite(tCo) & isfinite(tNo)), tNo(isfinite(tCo) & isfinite(tNo))]);
%         yCoMu = [yCoMu mean(z(1:length(z)/2))];
%         yNoMu = [yNoMu mean(z(length(z)/2+1:end))];
        yCoMu = [yCoMu mean(tCo)];
        yNoMu = [yNoMu mean(tNo)];

        tx = sprintf('%s %s', regions{iRa}, regions{iRb});
%         pl = plot(d, mean(tempCoR),'rs'); hold on;
%         pl.MarkerSize = 9;
%         pl.LineWidth =1.5;
%         pl.MarkerFaceColor = regionColors(strcmp(parcelA, regionsPlot),:);
% %         pl = plot(d+4, mean(tempCoR),'ks'); hold on;
% %         pl.MarkerSize = 9;
% %         pl.MarkerFaceColor = regionColors(strcmp(parcelB, regionsPlot),:);
% %         pl.LineWidth =1.5;
%         
%         tx = sprintf('%s %s', regions{iRa}, regions{iRb});
%         pl = plot(d, mean(tempNoR),'ks'); hold on;
%         pl.MarkerSize = 9;
%         pl.LineWidth =1.5;
%         pl.MarkerFaceColor = regionColors(strcmp(parcelB, regionsPlot),:);
% %         pl = plot(d+4, mean(tempNoR),'ks'); hold on;
% %         pl.MarkerSize = 9;
% %         pl.MarkerFaceColor = regionColors(strcmp(parcelB, regionsPlot),:) * 0.6;
% %         pl.LineWidth =1.5;


%         text(x(end), mean([yCo(end) yNo(end)]), tx); hold on;
        
    end
end
% z = zscore([yCo, yNo]);
% yCo = z(1:length(yCo));
% yNo = z(length(yCo)+1:end);

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
    pctDelta = yCo(ii); %[yCo(ii)-yNo(ii)];
    yCoMu = [yCoMu mean(pctDelta, 'omitnan')];
    yNoMu = [yNoMu mean(yNo(ii))];
    yCoSEM = [yCoSEM std(pctDelta, 'omitnan')/sqrt(sum(ii))];
    yNoSEM = [yNoSEM std(yNo(ii))/sqrt(sum(ii))];
    xMu = [xMu mean([xBin(iB-1) xBin(iB)])];
end



clr = brewermap(10,'Paired');


% plot(x, yCo,'r.'); hold on;
% plot(x(isfinite(yNo) & isfinite(yCo)), yNo(isfinite(yNo) & isfinite(yCo)),'k.'); hold on;
% plot(x(isfinite(yNo) & isfinite(yCo)), yCo(isfinite(yNo) & isfinite(yCo)),'r.'); hold on;
[bl2, bf] = boundedline(xMu, yCoMu, yCoSEM, 'rs-'); hold on;
bl2.Color = clr(8,:);
bl2.MarkerFaceColor = clr(8,:);
bf.FaceColor = clr(8,:);
bf.FaceAlpha = 0.3;

[bl2, bf] = boundedline(xMu, yNoMu, yNoSEM, 'bs-'); hold on;
bl2.Color = clr(1,:);
bl2.MarkerFaceColor = clr(1,:);
bf.FaceColor = clr(1,:);
bf.FaceAlpha = 0.3;

ax = gca;
ax.FontSize = 11;
ylh = ylabel('co-fire / $\sqrt{FR_a FR_b}$');
% ylh = ylabel('co-fire [z score]');
ylh.Interpreter = 'latex';

xlabel('distance [mm]')

fig = gcf;
fig.Color = 'w';


% plot(xMu(isfinite(yNoMu) & isfinite(yCoMu)), yNoMu(isfinite(yNoMu) & isfinite(yCoMu)),'ko'); hold on;
% plot(xMu(isfinite(yNoMu) & isfinite(yCoMu)), yCoMu(isfinite(yNoMu) & isfinite(yCoMu)),'ro'); hold on;
xlim([0 250])
[RHO,PVAL] = corr(xMu', yCoMu');
[RHO,PVAL] = corr(xMu', yNoMu');
[RHO,PVAL] = corr(x', yNo','type', 'Kendall')
% [RHO,PVAL] = corr(x(isfinite(yNo) & isfinite(yCo))', [yCo(isfinite(yNo) & isfinite(yCo)) - yNo(isfinite(yNo) & isfinite(yCo))]')
[RHO,PVAL] = corr(x', yCo', 'type', 'Kendall')
[RHO,PVAL] = corr(x', [yCo-yNo]', 'type', 'Kendall')

[H,P] = signrank(yCo,yNo, 'tail', 'both');
q = quantile((yCo-yNo)./yNo, [0.25 0.5 0.75]);

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireDistanceGeoMean_%s.pdf', tag)))
% savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireDistance_%s.pdf', tag)))


% Create a new figure
figure('Position', [100 368 round(535*(2/3)) round(232*(2/3))]);
num_variables = 6;
num_conditions = 2;
variable_names = cell(num_variables, 1);

% Create custom x positions for the boxes
% This is the key modification - we're creating custom positions with smaller gaps within variables
pos = zeros(num_variables * num_conditions, 1);
gap_within = 0.3;  % Small gap between conditions of the same variable
gap_between = 0.8; % Larger gap between different variables

for var = 1:num_variables
    base_pos = (var-1) * (gap_between);
    pos((var-1)*num_conditions + 1) = base_pos;
    pos((var-1)*num_conditions + 2) = base_pos + gap_within;
end

% Prepare data for the boxplot
all_data = [];
all_groups = [];

datCond1 = yCo(withinBundle == 1);
all_data = [all_data, datCond1];
idx = (1-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 1);
all_data = [all_data, datCond2];
idx = (1-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'AMY'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond1];
idx = (2-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'AMY'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond2];
idx = (2-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H,P] = signrank(datCond1,datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H) 

datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'HIP'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond1];
idx = (3-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'HIP'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}));
all_data = [all_data, datCond2];
idx = (3-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H,P] = signrank(datCond1,datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H) 

datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'HIP', 'AMY'}) & ismember(parcB, {'HIP', 'AMY'}));
all_data = [all_data, datCond1];
idx = (4-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'HIP', 'AMY'}) & ismember(parcB, {'HIP', 'AMY'}));
all_data = [all_data, datCond2];
idx = (4-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H,P] = signrank(datCond1,datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H) 

datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}) & strcmp(hemA, hemB));
all_data = [all_data, datCond1];
idx = (5-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}) & strcmp(hemA, hemB));
all_data = [all_data, datCond2];
idx = (5-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H,P] = signrank(datCond1,datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H) 

datCond1 = yCo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}) & ~strcmp(hemA, hemB));
all_data = [all_data, datCond1];
idx = (6-1)*num_conditions + 1;
all_groups = [all_groups; repmat(idx, length(datCond1), 1)];

datCond2 = yNo(withinBundle == 0 & ismember(parcA, {'OFC', 'ACC', 'SMA'}) & ismember(parcB, {'OFC', 'ACC', 'SMA'}) & ~strcmp(hemA, hemB));
all_data = [all_data, datCond2];
idx = (6-1)*num_conditions + 2;
all_groups = [all_groups; repmat(idx, length(datCond2), 1)];

[H,P] = signrank(datCond1,datCond2, 'tail', 'both');
Q = quantile((datCond1-datCond2)/median(datCond2), [0.25 0.5 0.75]);
fprintf('%i [%.2f %.2f %.2f] %e\n', length(datCond2), Q(1), Q(2), Q(3), H) 
% Create the boxplot with custom positions
% bx = boxplot(all_data, all_groups, 'positions', pos, 'width', 0.2, 'Colors', clr([8,1],:), 'Symbol', '', 'Whisker', 0);
bx = boxplot(all_data, all_groups, 'positions', pos, 'width', 0.2, 'Colors', [0 0 0], 'Symbol', '', 'Whisker', 0);

% Create custom x-tick positions and labels
tick_pos = zeros(num_variables, 1);
tick_labels = cell(num_variables, 1);

for var = 1:num_variables
    tick_pos(var) = pos((var-1)*num_conditions + 1) + gap_within/2;
    tick_labels{var} = variable_names{var};
end

% Set tick positions and labels
set(gca, 'XTick', tick_pos);
set(gca, 'XTickLabel', tick_labels);
set(gca, 'XLim', [min(pos)-0.5, max(pos)+0.5]);
set(gca, 'LineWidth', 0.75);
set(gca, 'FontSize', 7);
set(bx, 'LineWidth', 0.75);
% Add a legend for the conditions
h = findobj(gca, 'Tag', 'Box');
% legend([h(end), h(end-1)], condition_names, 'Location', 'Best');

% Customize box colors
for var = 1:num_variables
    h = findobj(gca, 'Tag', 'Box');
    for cond = 1:num_conditions
        idx = (var-1)*num_conditions + (num_conditions-cond+1);
        if cond == 1
            patch(get(h(idx), 'XData'), get(h(idx), 'YData'), clr(8,:), 'FaceAlpha', 0.5, 'LineWidth', 0.75);
        else
            patch(get(h(idx), 'XData'), get(h(idx), 'YData'), clr(1,:), 'FaceAlpha', 0.5, 'LineWidth', 0.75);
        end
    end
end



% Add title and labels
ylabel('co-firing rate [Hz]', 'FontSize', 7);

% Add statistical comparison lines
hold on;
% for var = 1:num_variables
%     [h, p] = ttest2(data{1}{var}, data{2}{var});
%     if p < 0.05
%         x1 = pos((var-1)*num_conditions + 1);
%         x2 = pos((var-1)*num_conditions + 2);
%         
%         % Find the maximum y value for this variable pair
%         max_val = max([data{1}{var}; data{2}{var}]);
%         
%         % Draw the significance line and asterisk
%         plot([x1, x2], [max_val*1.1, max_val*1.1], 'k-', 'LineWidth', 1);
%         text((x1+x2)/2, max_val*1.15, '*', 'FontSize', 15, 'HorizontalAlignment', 'center');
%     end
% end
hold off;

box off

fig = gcf;
fig.Color = 'w';


ylim([-0.1 1])

qWithin = quantile(yCo(withinBundle == 1), [0.25 0.5 0.75]);
qAcross = quantile(yCo(withinBundle == 0), [0.25 0.5 0.75]);

[H,P] = ranksum(yCo(withinBundle == 1),yCo(withinBundle == 0), 'tail', 'both');
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireSummary_%s.pdf', tag)))

%%
close all
plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
thresh = 0; 
coRegionColors = nan(10,3);

EdgeColors = [clr(8,:);
            mean([clr(1,:); clr(8,:)]);
              clr(1,:)];
          
% figure('Position',[1 891 round(525*1085/883) round(670*1334/1037)])
% figure('Position',[1 891 round(670*1085/883) round(670*1334/1037)])
figure('Position',[1 891 round(670*1085/883) round(670*1085/883)])
Pval = nan(length(regionsPlot), length(regionsPlot), 3, 3);
y_max = nan(length(regionsPlot), length(regionsPlot),3,3);
y_min = nan(length(regionsPlot), length(regionsPlot),3,3);
modulAll = nan(length(regionsPlot), length(regionsPlot),3,3);
modulAllPairs = cell(1,2);
bl_all = cell(length(regionsPlot), length(regionsPlot), 3, 3);
bf_all = cell(length(regionsPlot), length(regionsPlot), 3, 3);
spacing = 0.4;
for t = 4:-1:3
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
                    if iiA == iiB; continue; end
                    if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                    types = UnitPairType{iiA, iiB};
    %                 iiType = strcmp(types(1,:), 'pyr') &  strcmp(types(2,:), 'pyr');
                    iiType = strcmp(types(1,:), 'pyr') |  strcmp(types(1,:), 'int');
                    fprintf('%s %s \n', regions{iiA}, regions{iiB})
                    tempBase= [];
                    tempCoRcount= [];
                    tempNoR= [];
                    for ii = 1:2

                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,7+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,7+ii}(iiType,t);
    %                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoRcount(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,3+ii}(iiType,t);                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,5+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,5+ii}(iiType,t);

                    end 

    %                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];

    %                 all = [all tempAll];
                    coR = [coR tempCoRcount];
                    noR = [noR tempNoR];
                    base = [base tempBase];


                end
            end

            c = (iRb - 1) * 5 + iRa;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);

            figure(1)
            subplot(5,5,c)
            

%             keep = coR(1,:) > 0/70 | coR(2,:) > 0/70;
            keep = keepAll{iRa,iRb};
            
            
            baseModul = [base(1,keep)-base(2,keep)]./mean([base(1,keep) base(2,keep)]) * 100;
            coRmodul = [coR(1,keep)-coR(2,keep)]./mean([coR(1,keep) coR(2,keep)]) * 100;
            noRmodul = [noR(1,keep)-noR(2,keep)]./mean([noR(1,keep) noR(2,keep)]) * 100;
            modul = [coRmodul; baseModul; noRmodul];
            
            modulAllPairs{t-2} = [modulAllPairs{t-2}, modul];
    %         vp = violinplot([baseModul, coRmodul, noRmodul], groups); 
    %         for iV = 1:length(vp)
    %             vp(iV).ViolinPlot.Visible = 'off';
    %             vp(iV).ScatterPlot.Visible = 'off';
    %             vp(iV).WhiskerPlot.Visible = 'off';
    %             vp(iV).BoxWidth = 0.1;
    %             
    %         end
            YY = mean(modul, 2, 'omitnan');
             % Compute MAD

            YYdev = std(modul,[], 2, 'omitnan') / sqrt(sum(keep));
            xx = [(t-2)*3+1:(t-2)*3+3]+((t-2)*spacing);
            w = 0.30;
            cBox = [];
            cBox(1,:) = coRegionColors(c,:);
            cBox(2,:) = coRegionColors(c,:)*0.6;
            cBox(3,:) = coRegionColors(c,:)*0.3;
            for ii = 1:length(xx)

                bf_all{iRa,iRb,t-2,ii} = patch([xx(ii)-w xx(ii)+w xx(ii)+w xx(ii)-w], [YY(ii)-YYdev(ii) YY(ii)-YYdev(ii) YY(ii)+YYdev(ii) YY(ii)+YYdev(ii)], c); hold on;
                bf_all{iRa,iRb,t-2,ii}.FaceAlpha = 0.9;
                bf_all{iRa,iRb,t-2,ii}.FaceColor = cBox(ii,:);
                bf_all{iRa,iRb,t-2,ii}.EdgeAlpha = 1;
                bf_all{iRa,iRb,t-2,ii}.EdgeColor = 'k';  %EdgeColors(ii,:);

                bl_all{iRa,iRb,t-2,ii} = plot([xx(ii)-w xx(ii)+w], [YY(ii) YY(ii)], '-'); hold on;
                bl_all{iRa,iRb,t-2,ii}.LineWidth = 2;
                bl_all{iRa,iRb,t-2,ii}.Color = 'k'; %EdgeColors(ii,:); %cBox(ii,:);

            end
            y_max(iRa,iRb,1, t-1) = YY(1)+YYdev(1);
            y_max(iRa,iRb,2, t-1) = YY(2)+YYdev(2);
            y_max(iRa,iRb,3, t-1) = YY(3)+YYdev(3);
            y_min(iRa,iRb,1, t-1) = YY(1)-YYdev(1);
            y_min(iRa,iRb,2, t-1) = YY(2)-YYdev(2);
            y_min(iRa,iRb,3, t-1) = YY(3)-YYdev(3);
    %         errobarPatch([1:2], coR(1:end,keep)', coRegionColors(c,:), 0.1);

            xlim([3.3 9+(2*spacing)+w])
    %         ylim([-20 40])
    %         ax = gca;
    %         ax.XTick = [1:2];
    %         ax.XTickLabel = {'load3', 'load1'};
            fig = gcf;
            fig.Color = 'w';
    %         title(num2str(sum(keep)))
            hl = hline(0);
            hl.LineStyle = '-';
            hl.Color = 'k';
            hl.LineWidth = 1;
            uistack(hl, 'bottom')
            ax = gca;
            ax.XTick = [(3-2)*3+2+((3-2)*spacing) (4-2)*3+2+((4-2)*spacing)];
   box off
            
    %         figure(3)
    %         subplot(5,5,c)
    % 
    %         
    %         errobarPatch([1:2], noR(1:end,keep)', coRegionColors(c,:), 0.1);
    % 
    %         xlim([-1.5 3.5])
    %         ax = gca;
    %         ax.XTick = [1:2];
    %         ax.XTickLabel = {'load3', 'load1'};
    %         fig = gcf;
    %         fig.Color = 'w';
    %         
    %         figure(4)
    %         subplot(5,5,c)
    % 
    %         
    %         errobarPatch([1:2], base(1:end,keep)', coRegionColors(c,:), 0.1);
    % 
    %         xlim([-1.5 3.5])
    %         ax = gca;
    %         ax.XTick = [1:2];
    %         ax.XTickLabel = {'load3', 'load1'};
    %         fig = gcf;
    %         fig.Color = 'w';

    %         [stats, df, P, surrog] = statcond({base(1,keep), base(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
    %         Pval(iRa,iRb,1) = P;
    %         [stats, df, P, surrog] = statcond({noR(1,keep), noR(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
    %         Pval(iRa,iRb,2) = P;
    %         [stats, df, P, surrog] = statcond({coR(1,keep), coR(2,keep)},'paired','on','tail', 'one', 'method', 'perm', 'naccu', 10000, 'verbose','off');
    %         Pval(iRa,iRb,3) = P;
    %         
            [H,P,CI,STATS] = ttest(base(1,keep), base(2,keep),'tail', 'right');
            Pval(iRa,iRb,2,t-1) = P;
            [H,P,CI,STATS] = ttest(noR(1,keep), noR(2,keep),'tail', 'right');
            Pval(iRa,iRb,3,t-1) = P;
            [H,P,CI,STATS] = ttest(coR(1,keep), coR(2,keep), 'tail', 'right');
            Pval(iRa,iRb,1,t-1) = P;
            
            modulAll(iRa,iRb,1,t-1) = YY(1);
            modulAll(iRa,iRb,3,t-1) = YY(3);
            modulAll(iRa,iRb,2,t-1) = YY(2);
% 
%             [H,P,CI,STATS] = ttest(baseModul,zeros(1,length(baseModul)),'tail', 'right');
%             Pval(iRa,iRb,1,t-1) = P;
%             [H,P,CI,STATS] = ttest(noRmodul,zeros(1,length(noRmodul)),'tail', 'right');
%             Pval(iRa,iRb,3,t-1) = P;
%             [H,P,CI,STATS] = ttest(coRmodul,zeros(1,length(coRmodul)), 'tail', 'right');
%             Pval(iRa,iRb,2,t-1) = P;


    %         [P,H] = signrank(base(1,keep), base(2,keep),'tail', 'right');
    %         Pval(iRa,iRb,1) = P;
    %         [P,H] = signrank(noR(1,keep), noR(2,keep),'tail', 'right');
    %         Pval(iRa,iRb,2) = P;
    %         [P,H] = signrank(coR(1,keep), coR(2,keep), 'tail', 'right');
    %         Pval(iRa,iRb,3) = P;




             fprintf('\n')
            c = c+1;

        end
    end
end


adj_p = nan(size(Pval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(Pval))]=fdr_bh(Pval(~isnan(Pval)),0.05,'pdep','yes');


% Create color map for p-values with more distinct colors
pValueColorMap = [
    0.8, 0.8, 0.8;    % p > 0.05 (not significant) - light gray
    0.1, 0.7, 0.9;    % p < 0.05 - sky blue
    0.4, 0.0, 0.7;    % p < 0.01 - indigo
    0.8, 0.0, 0.0     % p < 0.001 - dark red
];

axRng = nan(length(regionsPlot), length(regionsPlot));
for t = 3:4
    for iRa = 1:length(regionsPlot)
        for iRb = iRa:length(regionsPlot)
            c = (iRb - 1) * 5 + iRa;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);
            figure(1)
            subplot(5,5,c)
            ax = gca;
            if t == 3
                axRng(iRa, iRb) = range(ax.YLim);
                ax.YLim(2) = ax.YLim(2) + 5/60*axRng(iRa, iRb);
                ax.YLim(1) = ax.YLim(1) - 5/30*axRng(iRa, iRb);
                ax.YTick = round(min(ax.YLim)/10)*10:10:max(ax.YLim);
                ax.YTick = unique([ax.YTick(1), 0, ax.YTick(end)]);
%                 ax.XTick = [2:3:9];
                ax.XTickLabel = {'M','P'};
                ax.FontSize = 9;
                ax.LineWidth = 1;
        
            end
            
            if t ==4
                vl = vline([7.1]);
                for v = 1:length(vl); vl(v).LineStyle = '--'; vl(v).Color = 'r';vl(v).LineWidth = 1; end
                
%                 vl = vline([1:9]); hold on;
%                 for v = 1:length(vl); vl(v).LineStyle = '-'; vl(v).Color = [0.5 0.5 0.5];vl(v).LineWidth = 0.4; end
%                 uistack(vl, 'bottom');
            end
            
            pReg = adj_p(iRa, iRb, :, t-1);
            for ii = 1:3
                if pReg(ii) < 0.05 && pReg(ii) >= 0.01
%                     bf_all{iRa,iRb,t-2,ii}.FaceColor = pValueColorMap(2,:);
%                     bl_all{iRa,iRb,t-2,ii}.Color = pValueColorMap(2, :);
                    tx = text((t-2)*3+ii+((t-2)*spacing), y_max(iRa, iRb, ii, t-1) + (1/9)*axRng(iRa, iRb) ,'*', 'FontSize', 13, ...
                              'HorizontalAlignment', 'center', ...
                              'Interpreter', 'none'); hold on
%                     tx.Color = pValueColorMap(2,:);
                elseif pReg(ii) < 0.01 && pReg(ii) >= 0.001
%                     bf_all{iRa,iRb,t-2,ii}.FaceColor = pValueColorMap(3,:);
%                     bl_all{iRa,iRb,t-2,ii}.Color = pValueColorMap(3, :);
                    tx = text((t-2)*3+ii+((t-2)*spacing), y_max(iRa, iRb, ii, t-1) + (1.5/9)*axRng(iRa, iRb) ,char(167), 'FontSize', 11, ...
                              'HorizontalAlignment', 'center', ...
                              'Interpreter', 'none'); hold on
%                     tx.Color = pValueColorMap(3,:);
                elseif pReg(ii) < 0.001 
%                     bf_all{iRa,iRb,t-2,ii}.FaceColor = pValueColorMap(4,:);
%                     bl_all{iRa,iRb,t-2,ii}.Color = pValueColorMap(4, :);
                    tx = text((t-2)*3+ii+((t-2)*spacing), y_max(iRa, iRb, ii, t-1) + (1/9)*axRng(iRa, iRb) ,char(9660), 'FontSize', 9, ...
                              'HorizontalAlignment', 'center', ...
                              'Interpreter', 'none'); hold on
%                     tx.Color = pValueColorMap(4,:);
                else
                    bf_all{iRa,iRb,t-2,ii}.FaceColor = pValueColorMap(1,:);
%                     bl_all{iRa,iRb,t-2,ii}.Color = pValueColorMap(1, :);
                    if modulAll(iRa,iRb,ii,t-1) >= 2
%                         tx = text((t-2)*3+ii+((t-2)*spacing), y_max(iRa, iRb, ii, t-1) + (1/9)*axRng(iRa, iRb) ,'ns', 'FontSize', 8, ...
%                                   'HorizontalAlignment', 'center', ...
%                                   'Interpreter', 'none'); hold on       
                    else
%                         tx = text((t-2)*3+ii+((t-2)*spacing), y_min(iRa, iRb, ii, t-1) - (1/9)*axRng(iRa, iRb) ,'ns', 'FontSize', 8, ...
%                                   'HorizontalAlignment', 'center', ...
%                                   'Interpreter', 'none'); hold on 
                    end
                end
            end


        end
    end
end

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskLoadAcrossRegions_%s.pdf', tag)))



conds = {'all', 'co-rip', 'no-rip'};
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
    ax.XTickLabel = regionsPlot;
    ax.YTickLabel = regionsPlot;
    
    title(conds{ii})



end
% sgtitle('co firing load 3 vs load 1')
fig= gcf;
fig.Color = 'w';

% figure('Position',[296 1084 400 132]); 
figure('Position',[1 842 299 289])
cBox = [];
cBox(1,:) = mean(regionColors(1:5,:));
cBox(2,:) = mean(regionColors(1:5,:))*0.6;
cBox(3,:) = mean(regionColors(1:5,:))*0.3;
cBox(1,:) = [1 0 0];
cBox(2,:) = [1 0 0]*0.3;
cBox(3,:) = [1 0 0]*0.6;
cnt = 1;
for t = 2:4
    for ii = [1,2,3]
        N = sum(adj_p(:,:,ii,t-1) < 0.05 & modulAll(:,:,ii,t-1) > 0,'all');
        b = bar(2*(t-1)+0.5*ii, N,0.35); hold on;
        b.FaceColor = EdgeColors(ii,:);
        b.EdgeColor = 'k';
        b.LineWidth = 1;
%         N = sum(adj_p(:,:,ii,t-1) < 0.05 & modulAll(:,:,ii,t-1) < 0,'all');
%         b = bar(cnt, -N,0.7); hold on;
%         b.FaceColor = cBox(ii,:);
%         b.EdgeColor = 'b';
%         b.LineWidth = 1;

        cnt = cnt + 1;

    end
end
ylim([0 15])
xlim([3.8 8.0])
box off

ax = gca;
ax.XTick = [2.5:0.5:3.5, 4.5:0.5:5.5, 6.5:0.5:7.5];
ax.XTickLabel = {'','',''};
ax.LineWidth = 1.0;
ax.FontSize = 11;
ax.YTickLabel = abs(ax.YTick);
fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskLoadAcrossRegionsBar_%s.pdf', tag)))

% close all
plot_condition = load_trials_all == 3;
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(~plot_condition, 1, 'first');
thresh = 0; 
coRegionColors = nan(10,3);

EdgeColors = [clr(8,:);
            mean([clr(1,:); clr(8,:)]);
              clr(1,:)];
          
% figure('Position',[1 891 round(525*1085/883) round(670*1334/1037)])
% figure('Position',[1 891 round(670*1085/883) round(670*1334/1037)])
Pval = nan(length(regionsPlot), length(regionsPlot), 3, 3);
y_max = nan(length(regionsPlot), length(regionsPlot),3,3);
y_min = nan(length(regionsPlot), length(regionsPlot),3,3);
modulAll = nan(length(regionsPlot), length(regionsPlot),3,3);
modulAllPairs = cell(1,2);
bl_all = cell(length(regionsPlot), length(regionsPlot), 3, 3);
bf_all = cell(length(regionsPlot), length(regionsPlot), 3, 3);
spacing = 0.4;
strng = {'M', 'P'};
for t = 3:4
    figure('Position',[1 847 410 296])
    xx_jitt = [];
    YY = [];
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
                    if iiA == iiB; continue; end
                    if isempty(whole_trials_unit_cross_region_all{iiA, iiB,3}); continue; end
                    types = UnitPairType{iiA, iiB};
    %                 iiType = strcmp(types(1,:), 'pyr') &  strcmp(types(2,:), 'pyr');
                    iiType = strcmp(types(1,:), 'pyr') |  strcmp(types(1,:), 'int');
                    fprintf('%s %s \n', regions{iiA}, regions{iiB})
                    tempBase= [];
                    tempCoRcount= [];
                    tempNoR= [];
                    for ii = 1:2

                        tempBase(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,7+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,7+ii}(iiType,t);
    %                     tempAll = whole_trials_unit_cross_region_all{iiA, iiB,3}(:,ii) ./ whole_trials_cross_region_dur{iiA, iiB,3}(:,4) * 1e3;                
                        tempCoRcount(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,3+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,3+ii}(iiType,t);                
                        tempNoR(ii,:) = whole_trials_unit_cross_region_all{iiA, iiB,5+ii}(iiType,t) ./ whole_trials_cross_region_dur{iiA, iiB,5+ii}(iiType,t);

                    end 

    %                 all = [all tempAll(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 coR = [coR tempCoR(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 noR = [noR tempNoR(tempBase > 0.00 & ~isnan(tempCoR))'];
    %                 base = [base tempBase(tempBase > 0.00 & ~isnan(tempCoR))];

    %                 all = [all tempAll];
                    coR = [coR tempCoRcount];
                    noR = [noR tempNoR];
                    base = [base tempBase];


                end
            end

            c = (iRb - 1) * 5 + iRa;
            coRegionColors(c,:) = mean([regionColors(iRa,:);regionColors(iRb,:)]);

            figure(t+1)
%             subplot(5,5,c)
            

%             keep = coR(1,:) > 0/70 | coR(2,:) > 0/70;
            keep = keepAll{iRa,iRb};
            
            
            baseModul = [base(1,keep)-base(2,keep)]./mean([base(1,keep) base(2,keep)]) * 100;
            coRmodul = [coR(1,keep)-coR(2,keep)]./mean([coR(1,keep) coR(2,keep)]) * 100;
            noRmodul = [noR(1,keep)-noR(2,keep)]./mean([noR(1,keep) noR(2,keep)]) * 100;
            modul = [coRmodul; baseModul; noRmodul];
            
            modulAllPairs{t-2} = [modulAllPairs{t-2}, modul];
    %         vp = violinplot([baseModul, coRmodul, noRmodul], groups); 
    %         for iV = 1:length(vp)
    %             vp(iV).ViolinPlot.Visible = 'off';
    %             vp(iV).ScatterPlot.Visible = 'off';
    %             vp(iV).WhiskerPlot.Visible = 'off';
    %             vp(iV).BoxWidth = 0.1;
    %             
    %         end
            YY(:,iRa,iRb) = mean(modul, 2, 'omitnan');
             % Compute MAD

            YYdev = std(modul,[], 2, 'omitnan') / sqrt(sum(keep));
            xx = [(t-2)*3+1:(t-2)*3+3]+((t-2)*spacing);
            w = 0.10;
            cBox = [];
            cBox(1,:) = coRegionColors(c,:);
            cBox(2,:) = coRegionColors(c,:)*0.6;
            cBox(3,:) = coRegionColors(c,:)*0.3;
            
            xx_jitt(:,iRa,iRb) = (1:3) + ((rand(1) - 0.5)*0.25);
            pl = plot(xx_jitt(:,iRa,iRb),YY(:,iRa,iRb), '-'); hold on;
            pl.Color = coRegionColors(c,:);
            pl.LineWidth = 1.5;
            
            
            
   

            xlim([0 3.5])
    %         ylim([-20 40])
    %         ax = gca;
    %         ax.XTick = [1:2];
    %         ax.XTickLabel = {'load3', 'load1'};
            fig = gcf;
            fig.Color = 'w';
    %         title(num2str(sum(keep)))
            hl = hline(0);
            hl.LineStyle = '-';
            hl.Color = 'k';
            hl.LineWidth = 1;
            uistack(hl, 'bottom')
            ax = gca;
            ax.XTick = [(3-2)*3+2+((3-2)*spacing) (4-2)*3+2+((4-2)*spacing)];
   box off
            




             fprintf('\n')
            c = c+1;

        end
    end
    
    axRng = nan(length(regionsPlot), length(regionsPlot));
    for iRa = 1:length(regionsPlot)
        for iRb = iRa:length(regionsPlot)
            c = (iRb - 1) * 5 + iRa;

            

            pReg = adj_p(iRa, iRb, :, t-1);
            for ii = 1:3
                pl = plot(xx_jitt(ii,iRa,iRb),YY(ii, iRa, iRb), 'o'); hold on;
                pl.MarkerSize = 7;
                pl.LineWidth = 1.5;
                pl.Color = coRegionColors(c,:);
                if pReg(ii) < 0.05 && pReg(ii) >= 0.01
                    pl.MarkerFaceColor = pValueColorMap(2,:);

                elseif pReg(ii) < 0.01 && pReg(ii) >= 0.001
                    pl.MarkerFaceColor = pValueColorMap(3,:);

                elseif pReg(ii) < 0.001 
                    pl.MarkerFaceColor = pValueColorMap(4,:);

                else
                    pl.MarkerFaceColor = pValueColorMap(1,:);

                end
            end


        end
    end
    ax = gca;
    ax.XTick = 1:3;
    ax.XTickLabel = {'coR', 'all', 'noR'};
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
    
    savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipCoFireTaskLoadAcrossRegionsSingle%s_%s.pdf', strng{t-2},    tag)))

    

end

%%

mean(modulAllPairs{1}(1,:))
[H,P,CI,STATS] = ttest(modulAllPairs{1}(1,:))








