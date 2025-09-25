

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

regionsPlot = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  


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

    trialsProbe = zeros(length(regions),length(regions),1e3);

    coR_replay = zeros(2,length(regions),length(regions));
    coR_replay_shuff = zeros(2,length(regions),length(regions));
    noR_replay = zeros(2,length(regions),length(regions));
    noR_replay_shuff = zeros(2,length(regions),length(regions));
    full_replay = zeros(2,length(regions),length(regions));
    full_replay_mismatch = zeros(2,5);
    fullCoF = nan(5e5, 2, length(regions),length(regions));
    fullCoF_mm = nan(5e5, 2, length(regions),length(regions));
    
    replay_event = cell(1,8);

    coRdur = zeros(2,2,length(regions),length(regions));
    noRdur = zeros(2,2,length(regions),length(regions));
    coRap = zeros(2,2,length(regions),length(regions));
    noRap = zeros(2,2,length(regions),length(regions));
    loads = nan(5e5, length(regions),length(regions));
    respTimes = nan(5e5, length(regions),length(regions));
    coRcoF = nan(5e5, 2, length(regions),length(regions));
    noRcoF = nan(5e5, 2, length(regions),length(regions), 25);
    coRcoFShuff = nan(5e5, 2, length(regions),length(regions));
    
    coR_replay_mm = zeros(2,length(regions),length(regions));
    noR_replay_mm = zeros(2,length(regions),length(regions));
    
    coRdur_mm = zeros(2,2,length(regions),length(regions));
    noRdur_mm = zeros(2,2,length(regions),length(regions));
    coRap_mm = zeros(2,2,length(regions),length(regions));
    noRap_mm = zeros(2,2,length(regions),length(regions));
    coRcoF_mm = nan(5e5, 2, length(regions),length(regions));
    noRcoF_mm = nan(5e5, 2, length(regions),length(regions), 25);
    
    noRcoFCount = ones(length(regions),length(regions), 25);
    noRcoF_mmCount = ones(length(regions),length(regions), 25);
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
        


         %across region co-rippling
        for iRseed = 1 %:length(regionsPlot)
            iiSeed = find(contains(regions, regions{iRseed}));
            iiB = find(~contains(regions, regions{iRseed}));

            for iRa = 1:length(regions)
                for iRb = iRa+1:length(regions) %iiB
%                     iRa = iRa;
                    datArip =  whole_trials_region_all(cTrial,:,iRa);
                    datBrip =  whole_trials_region_all(cTrial,:,iRb);

                    if iRa == iRb

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


        %                 
                    datArip = sum(rippMask(contains(locations, regions(iRa)),:), 'omitnan') > 0;
                    datBrip = sum(rippMask(contains(locations, regions(iRb)),:), 'omitnan') > 0;

                    if sum(datArip) == 0 || sum(datBrip) == 0; continue; end


                    nUnits = length(uLFPregionA) + length(uLFPregionB);
                    probeIm = []; encodIm =[];
                    for iT = 1:length(trials.start_time)
                        trialLoad = trials.loads(iT);
                        trialAcc = trials.response_accuracy;

                        imageTime(1) =  round(trials.timestamps_Encoding1(iT)*1e3);
                        imageTime(2) =  round(trials.timestamps_Encoding2(iT)*1e3);
                        imageTime(3) =  round(trials.timestamps_Encoding3(iT)*1e3);
                        imageTime(4) =  round(trials.timestamps_Maintenance(iT)*1e3);
                        probeTime =  round(trials.timestamps_Probe(iT)*1e3);
                        respTime =  round(trials.timestamps_Response(iT)*1e3);



                        probeIm(iT) = trials.PicIDs_Probe(iT);
                        encodIm(1, iT) = trials.PicIDs_Encoding1(iT);
                        encodIm(2, iT) = trials.PicIDs_Encoding2(iT);
                        encodIm(3, iT) = trials.PicIDs_Encoding3(iT);

                        if (trials.probe_in_out(iT)) 

                            for uA = 1:length(uRegionA)
                                for uB = 1:length(uRegionB)
                                    typeA = units{uRegionA(uA),3};
                                    typeB = units{uRegionB(uB),3};
                                    spikeMaskP = spikeMask(:,probeTime+1:probeTime+1e3);
                                    datABcoripP = datArip(probeTime+1:probeTime+1e3) > 0 & datBrip(probeTime+1:probeTime+1e3) > 0;
                                    datABnoripP = datArip(probeTime+1:probeTime+1e3) == 0 & datBrip(probeTime+1:probeTime+1e3) == 0;
                                    if trials.loads(iT) == 3                                                      
                                        iE = find(encodIm(:,iT) == probeIm(iT));
                                        iEmm = find(encodIm(:,iT) ~= probeIm(iT));
                                        iEmm = iEmm(randi(length(iEmm)));
        %                                 spikeMaskE = spikeMask(:,imageTime(iE):imageTime(iE+1));
                                        spikeMaskE = spikeMask(:,imageTime(iE):imageTime(iE)+2e3);
                                        spikeMaskEmm = spikeMask(:,imageTime(iEmm):imageTime(iEmm)+2e3);
                                        datABcoripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & datBrip(imageTime(iE):imageTime(iE)+2e3);
                                        datABnoripE = datArip(imageTime(iE):imageTime(iE)+2e3) == 0 & datBrip(imageTime(iE):imageTime(iE)+2e3) == 0;
                                    else
                                        spikeMaskE = spikeMask(:,imageTime(1):imageTime(1)+2e3);
                                        datABcoripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & datBrip(imageTime(1):imageTime(1)+2e3) > 0;
                                        datABnoripE = datArip(imageTime(1):imageTime(1)+2e3) == 0 & datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                                    end

                                    if sum(datABcoripP) == 0 || sum(datABcoripE) == 0; continue; end
%                                     if sum(datABnoripP) == 0 || sum(datABnoripE) == 0; continue; end

                                    %co-ripples
                                    datAP =  find(spikeMaskP(uRegionA(uA),datABcoripP)); datBP =  find(spikeMaskP(uRegionB(uB),datABcoripP));                                  
                                    sP = datAP - datBP'; 
                                    coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);% & sP(:) <= coFireWindow*2);
                                    
                                    
                                    datAE =  find(spikeMaskE(uRegionA(uA),datABcoripE)); datBE =  find(spikeMaskE(uRegionB(uB),datABcoripE)); 
                                    sE = datAE - datBE';
                                    coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);% & sE(:) <= coFireWindow*2);
                                   
                                    
                                    ii = find(isnan(coRcoF(:,1,iRa, iRb)), 1, 'first');
                                    coRcoF(ii,1,iRa, iRb) = coFp;
                                    coRcoF(ii,2,iRa, iRb) = coFe;
                                    loads(ii,iRa, iRb) = trialLoad;
                                    respTimes(ii,iRa, iRb) = respTime - probeTime;
                                    
                                    if coFp > 0 && coFe > 0
                                        coR_replay(1,iRa, iRb) = coR_replay(1,iRa, iRb) + 1;
                                        coRdur(1,1,iRa, iRb) = coRdur(1,1,iRa, iRb) + sum(datABcoripE);
                                        coRdur(1,2,iRa, iRb) = coRdur(1,2,iRa, iRb) + sum(datABcoripP);
                                        coRap(1,1,iRa, iRb) = coRap(1,1,iRa, iRb) + sum(coFe);
                                        coRap(1,2,iRa, iRb) = coRap(1,2,iRa, iRb) + sum(coFp);
                                        replay_event{1} = [replay_event{1} subj];
                                        replay_event{2} = [replay_event{2} uRegionA(uA)];
                                        replay_event{3} = [replay_event{3} uRegionB(uB)];  
                                        replay_event{4} = [replay_event{4} iT];
                                        replay_event{5} = [replay_event{5} {datAP}];
                                        replay_event{6} = [replay_event{6} {datBP}];
                                        replay_event{7} = [replay_event{7} {datAE}];
                                        replay_event{8} = [replay_event{8} {datBE}];
                                    else
                                        coR_replay(2,iRa, iRb) = coR_replay(2,iRa, iRb) + 1;
                                        coRdur(2,1,iRa, iRb) = coRdur(2,1,iRa, iRb) + sum(datABcoripE);
                                        coRdur(2,2,iRa, iRb) = coRdur(2,2,iRa, iRb) + sum(datABcoripP);
                                        coRap(2,1,iRa, iRb) = coRap(2,1,iRa, iRb) + sum(coFe);
                                        coRap(2,2,iRa, iRb) = coRap(2,2,iRa, iRb) + sum(coFp);
                                    end
                                    
                                    %co-ripple shuffle

                                    datAshuff = randi(sum(datABcoripP)+sum(datABcoripE), [1, length(datAP)+length(datAE)]); 
                                    datBshuff = randi(sum(datABcoripP)+sum(datABcoripE), [1, length(datBP)+length(datBE)]);
                                    A = datAshuff(datAshuff<=sum(datABcoripP));
                                    B = datBshuff(datBshuff<=sum(datABcoripP));
                                    if isempty(A) || isempty(B)
                                        sPshuff = []; 
                                    else
                                        sPshuff = A - B'; 
                                    end
                                    coFpshuff = sum(sPshuff(:) >= -coFireWindow*2 & sPshuff(:) <= coFireWindow*2);
                                    
                                    A = datAshuff(datAshuff>sum(datABcoripP));
                                    B = datBshuff(datBshuff>sum(datABcoripP));
                                    if isempty(A) || isempty(B)
                                        sEshuff = []; 
                                    else
                                        sEshuff = A - B'; 
                                    end
                                    coFeshuff = sum(sEshuff(:) >= -coFireWindow*2 & sEshuff(:) <= coFireWindow*2);
                                    
                                    coRcoFShuff(ii,1,iRa, iRb) = coFpshuff;
                                    coRcoFShuff(ii,2,iRa, iRb) = coFeshuff;
        
                                    if coFpshuff > 0 && coFeshuff > 0
                                        coR_replay_shuff(1,iRa, iRb) = coR_replay_shuff(1,iRa, iRb) + 1;

                                    else
                                        coR_replay_shuff(2,iRa, iRb) = coR_replay_shuff(2,iRa, iRb) + 1;
                                    end
                                    
                                    for iter = 1:25
                                        %no-ripples
                                        ripDur = sum(datABcoripP);
                                        b = mask2bounds(datABnoripP);
                                        bDur = b(:,2) - b(:,1);
                                        if isempty(b) 
                                            b = [1,1]; ripDur = 0;
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
                                        datA =  find(spikeMaskP(uRegionA(uA),datABnoripPiter)); datB =  find(spikeMaskP(uRegionB(uB),datABnoripPiter));
    %                                     datAshuff = randi(sum(datABnoripP), [1, length(datA)]); datBshuff = randi(sum(datABnoripP), [1, length(datB)]);
                                        sP = datA - datB'; 
                                        noFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);% & sP(:) <= coFireWindow*2);
    %                                     sPshuff = datAshuff - datBshuff'; 
    %                                     coFpshuff = sum(sPshuff(:) >= -coFireWindow*2 & sPshuff(:) <= coFireWindow*2);

                                        ripDur = sum(datABcoripE);
                                        b = mask2bounds(datABnoripE);
                                        bDur = b(:,2) - b(:,1);
                                        if isempty(b) 
                                            b = [1,1]; ripDur = 0;
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
                                        datA =  find(spikeMaskE(uRegionA(uA),datABnoripEiter)); datB =  find(spikeMaskE(uRegionB(uB),datABnoripEiter));
    %                                     datAshuff = randi(sum(datABnoripE), [1, length(datA)]); datBshuff = randi(sum(datABnoripE), [1, length(datB)]);
                                        sE = datA - datB';
                                        noFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);% & sE(:) <= coFireWindow*2);
    %                                     sEshuff = datAshuff - datBshuff'; 
    %                                     coFeshuff = sum(sEshuff(:) >= -coFireWindow*2 & sEshuff(:) <= coFireWindow*2);
                                    
                                        ii = noRcoFCount(iRa, iRb,iter);
                                        noRcoF(ii,1,iRa, iRb, iter) = noFp;
                                        noRcoF(ii,2,iRa, iRb, iter) = noFe;
                                        noRcoFCount(iRa, iRb,iter) = noRcoFCount(iRa, iRb,iter) + 1;
                                        if noFp > 0 && noFe > 0
                                            noRdur(1,1,iRa, iRb) = noRdur(1,1,iRa, iRb) + sum(datABnoripEiter);
                                            noRdur(1,2,iRa, iRb) = noRdur(1,2,iRa, iRb) + sum(datABnoripPiter);
                                            noRap(1,1,iRa, iRb) = noRap(1,1,iRa, iRb) + sum(noFe);
                                            noRap(1,2,iRa, iRb) = noRap(1,2,iRa, iRb) + sum(noFp);
                                        else
                                            noRdur(2,1,iRa, iRb) = noRdur(2,1,iRa, iRb) + sum(datABnoripEiter);
                                            noRdur(2,2,iRa, iRb) = noRdur(2,2,iRa, iRb) + sum(datABnoripPiter);
                                            noRap(2,1,iRa, iRb) = noRap(2,1,iRa, iRb) + sum(noFe);
                                            noRap(2,2,iRa, iRb) = noRap(2,2,iRa, iRb) + sum(noFp);
                                        end
                                    end
                                    
                                    if noFp > 0 && noFe > 0
                                        noR_replay(1,iRa, iRb) = noR_replay(1,iRa, iRb) + 1;
                                    else
                                        noR_replay(2,iRa, iRb) = noR_replay(2, iRa, iRb) + 1;
                                    end
                                    
%co-ripples
                                    datAP =  find(spikeMaskP(uRegionA(uA),:)); datBP =  find(spikeMaskP(uRegionB(uB),:));                                  
                                    sP = datAP - datBP'; 
                                    coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);%);% & sP(:) <= coFireWindow*2);
                                    
                                    
                                    datAE =  find(spikeMaskE(uRegionA(uA),:)); datBE =  find(spikeMaskE(uRegionB(uB),:)); 
                                    sE = datAE - datBE';
                                    coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2); % & sE(:) <= coFireWindow*2);
                                   
                                    
                                    ii = find(isnan(fullCoF(:,1,iRa, iRb)), 1, 'first');
                                    fullCoF(ii,1,iRa, iRb) = coFp;
                                    fullCoF(ii,2,iRa, iRb) = coFe;

                                end
                            end

                        else
                            for uA = 1:length(uRegionA)
                                for uB = 1:length(uRegionB)
                                    typeA = units{uRegionA(uA),3};
                                    typeB = units{uRegionB(uB),3};
                                    spikeMaskP = spikeMask(:,probeTime+1:probeTime+1e3);
                                    datABcoripP = datArip(probeTime+1:probeTime+1e3) > 0 & datBrip(probeTime+1:probeTime+1e3) > 0;
                                    datABnoripP = datArip(probeTime+1:probeTime+1e3) == 0 & datBrip(probeTime+1:probeTime+1e3) == 0;
                                    if trials.loads(iT) == 3                                                      
                                        iEmm = find(encodIm(:,iT) ~= probeIm(iT));
                                        iEmm = iEmm(randi(length(iEmm)));
                                        spikeMaskEmm = spikeMask(:,imageTime(iEmm):imageTime(iEmm)+2e3);
                                        datABcoripEmm = datArip(imageTime(iEmm):imageTime(iEmm)+2e3) > 0 & datBrip(imageTime(iEmm):imageTime(iEmm)+2e3);
                                        datABnoripEmm = datArip(imageTime(iEmm):imageTime(iEmm)+2e3) == 0 & datBrip(imageTime(iEmm):imageTime(iEmm)+2e3) == 0;
                                    else
                                        spikeMaskEmm  = spikeMask(:,imageTime(1):imageTime(1)+2e3);
                                        datABcoripEmm = datArip(imageTime(1):imageTime(1)+2e3) > 0 & datBrip(imageTime(1):imageTime(1)+2e3) > 0;
                                        datABnoripEmm = datArip(imageTime(1):imageTime(1)+2e3) == 0 & datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                                    end

                                    if sum(datABcoripP) == 0 || sum(datABcoripE) == 0; continue; end
%                                     if sum(datABnoripP) == 0 || sum(datABnoripE) == 0; continue; end

                                    %co-ripples
                                    datA =  find(spikeMaskP(uRegionA(uA),datABcoripP)); datB =  find(spikeMaskP(uRegionB(uB),datABcoripP));
                                    sP = datA - datB'; 
                                    coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);% & sP(:) <= coFireWindow*2);

                                    datA =  find(spikeMaskEmm(uRegionA(uA),datABcoripEmm)); datB =  find(spikeMaskEmm(uRegionB(uB),datABcoripEmm));
                                    sE = datA - datB';
                                    coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);% & sE(:) <= coFireWindow*2);
                                    
                                    ii = find(isnan(coRcoF_mm(:,1,iRa, iRb)), 1, 'first');
                                    coRcoF_mm(ii,1,iRa, iRb) = coFp;
                                    coRcoF_mm(ii,2,iRa, iRb) = coFe;
                                    if coFp > 0 && coFe > 0
                                        coR_replay_mm(1,iRa, iRb) = coR_replay_mm(1,iRa, iRb) + 1;
                                        coRdur_mm(1,1,iRa, iRb) = coRdur_mm(1,1,iRa, iRb) + sum(datABcoripEmm);
                                        coRdur_mm(1,2,iRa, iRb) = coRdur_mm(1,2,iRa, iRb) + sum(datABcoripP);
                                        coRap_mm(1,1,iRa, iRb) = coRap_mm(1,1,iRa, iRb) + sum(coFe);
                                        coRap_mm(1,2,iRa, iRb) = coRap_mm(1,2,iRa, iRb) + sum(coFp);

                                    else
                                        coR_replay_mm(2,iRa, iRb) = coR_replay_mm(2,iRa, iRb) + 1;
                                        coRdur_mm(2,1,iRa, iRb) = coRdur_mm(2,1,iRa, iRb) + sum(datABcoripEmm);
                                        coRdur_mm(2,2,iRa, iRb) = coRdur_mm(2,2,iRa, iRb) + sum(datABcoripP);
                                        coRap_mm(2,1,iRa, iRb) = coRap_mm(2,1,iRa, iRb) + sum(coFe);
                                        coRap_mm(2,2,iRa, iRb) = coRap_mm(2,2,iRa, iRb) + sum(coFp);
                                    end


        %                             end

                                    %no-ripples
                                    for iter = 1:25
                                        ripDur = sum(datABcoripP);
                                        b = mask2bounds(datABnoripP);
                                        bDur = b(:,2) - b(:,1);
                                        if isempty(b) 
                                            b = [1,1]; ripDur = 0;
                                        elseif sum(bDur >= ripDur) == 0
                                            b = b(bDur == max(bDur),:);
                                            ripDur = max(bDur);
                                        elseif sum(bDur >= ripDur) > 0
                                            b(bDur < ripDur,:) = [];
                                            b = b(randi(size(b,1),1),:);
                                        end
                                        
                                        datABnoripPiter = datABnoripP;
                                        b = b(randi(size(b,1),1),:);
                                        randStart = randi([b(1), b(2) - ripDur]);
                                        datABnoripPiter(1:randStart) = false;
                                        datABnoripPiter(randStart+ripDur+1:end) = false;
                                        datA =  find(spikeMaskP(uRegionA(uA),datABnoripPiter)); datB =  find(spikeMaskP(uRegionB(uB),datABnoripPiter));
                                        sP = datA - datB'; 
                                        noFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);% & sP(:) <= coFireWindow*2);

                                        ripDur = sum(datABcoripEmm);
                                        b = mask2bounds(datABnoripEmm);
                                        bDur = b(:,2) - b(:,1);
                                        if isempty(b) 
                                            b = [1,1]; ripDur = 0;
                                        elseif sum(bDur >= ripDur) == 0
                                            b = b(bDur == max(bDur),:);
                                            ripDur = max(bDur);
                                        elseif sum(bDur >= ripDur) > 0
                                            b(bDur < ripDur,:) = [];
                                            b = b(randi(size(b,1),1),:);
                                        end

                                        datABnoripEmmiter = datABnoripEmm;
                                        randStart = randi([b(1), b(2) - ripDur]);
                                        datABnoripEmmiter(1:randStart) = false;
                                        datABnoripEmmiter(randStart+ripDur+1:end) = false;
                                        datA =  find(spikeMaskEmm(uRegionA(uA),datABnoripEmmiter)); datB =  find(spikeMaskEmm(uRegionB(uB),datABnoripEmmiter));
                                        sE = datA - datB';
                                        noFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);% & sE(:) <= coFireWindow*2);

                                        ii = noRcoF_mmCount(iRa, iRb,iter);
                                        noRcoF_mm(ii,1,iRa, iRb, iter) = noFp;
                                        noRcoF_mm(ii,2,iRa, iRb, iter) = noFe;
                                        noRcoF_mmCount(iRa, iRb,iter) = noRcoF_mmCount(iRa, iRb,iter) +1;
                                        if noFp > 0 && noFe > 0
                                            noR_replay_mm(1,iRa, iRb) = noR_replay_mm(1,iRa, iRb) + 1;
                                            noRdur_mm(1,1,iRa, iRb) = noRdur_mm(1,1,iRa, iRb) + sum(datABnoripEmmiter);
                                            noRdur_mm(1,2,iRa, iRb) = noRdur_mm(1,2,iRa, iRb) + sum(datABnoripPiter);
                                            noRap_mm(1,1,iRa, iRb) = noRap_mm(1,1,iRa, iRb) + sum(noFe);
                                            noRap_mm(1,2,iRa, iRb) = noRap_mm(1,2,iRa, iRb) + sum(noFp);
                                        else
                                            noRdur_mm(2,1,iRa, iRb) = noRdur_mm(2,1,iRa, iRb) + sum(datABnoripEmmiter);
                                            noRdur_mm(2,2,iRa, iRb) = noRdur_mm(2,2,iRa, iRb) + sum(datABnoripPiter);
                                            noRap_mm(2,1,iRa, iRb) = noRap_mm(2,1,iRa, iRb) + sum(noFe);
                                            noRap_mm(2,2,iRa, iRb) = noRap_mm(2,2,iRa, iRb) + sum(noFp);
                                        end
                                    end
                                    if noFp > 0 && noFe > 0
                                        noR_replay_mm(1,iRa, iRb) = noR_replay_mm(1,iRa, iRb) + 1;
                                    else
                                        noR_replay_mm(2,iRa, iRb) = noR_replay_mm(2,iRa, iRb) + 1;
                                    end   
                                    
                                    datAP =  find(spikeMaskP(uRegionA(uA),:)); datBP =  find(spikeMaskP(uRegionB(uB),:));                                  
                                    sP = datAP - datBP'; 
                                    coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);% & sP(:) <= coFireWindow*2);
                                    
                                    
                                    datAE =  find(spikeMaskEmm(uRegionA(uA),:)); datBE =  find(spikeMaskEmm(uRegionB(uB),:)); 
                                    sE = datAE - datBE';
                                    coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);% & sE(:) <= coFireWindow*2);
                                   
                                    
                                    ii = find(isnan(fullCoF_mm(:,1,iRa, iRb)), 1, 'first');
                                    fullCoF_mm(ii,1,iRa, iRb) = coFp;
                                    fullCoF_mm(ii,2,iRa, iRb) = coFe;
                                    
                                end %uB
                            end %uA
                        end
                    end %trials
                end%iRb
            end %iRa
        end %iRseed
    end %subjct

    
    correct_trials_all(cTrial:end) = [];
    load_trials_all(cTrial:end) =[];

    probe_in_out(cTrial:end) = [];
    subjID(cTrial:end) = [];

    fprintf('done processing .. \n')
    toc
%     save(sprintf('Sternbger_ripple_results_z25_%icoF_splitHem_template_v2.mat', coFireWindow*2), '-v7.3')
end    
% %% pairwise split hem
% 
% close all
% 
% load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
% figure;
% figure;
% figure;
% pValue = nan(10,10);
% for iRa = 1:10
%     figure(3)
%     regionCoordsAll{iRa} = double(regionCoordsAll{iRa});
% %     plot(mean(regionCoordsAll{iRa}(:,1)),mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,3)), 'o'); hold on;
%     if (strcmp(regionsAll{iRa}(1), 'L') && ~contains(regionsAll{iRa}, 'OFC')) || (strcmp(regionsAll{iRa}(1), 'R') && contains(regionsAll{iRa}, 'OFC'))
%         plot3(mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3)), 'ko'); hold on;
%         text(mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3)),regionsAll{iRa})
%         xA = [mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3))];
%     else
%         reflectedPoints = 2*40 - regionCoordsAll{iRa}(:,2);
%         plot3(mean(reflectedPoints), mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3)), 'go'); hold on;
%         xA = [mean(mean(reflectedPoints)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3))];
%         text(mean(reflectedPoints), mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3)),regionsAll{iRa})
% 
%     end
%     
%     
%     for iRb = iRa+1:10
%         
%          if (strcmp(regionsAll{iRb}(1), 'L') && ~contains(regionsAll{iRb}, 'OFC')) || (strcmp(regionsAll{iRb}(1), 'R') && contains(regionsAll{iRb}, 'OFC'))
%             xB = [mean(regionCoordsAll{iRb}(:,2)),mean(regionCoordsAll{iRb}(:,1)), mean(regionCoordsAll{iRb}(:,3))];
%         else
%             reflectedPoints = 2*40 - regionCoordsAll{iRb}(:,2);
%             xB = [mean(mean(reflectedPoints)),mean(regionCoordsAll{iRb}(:,1)),mean(regionCoordsAll{iRb}(:,3))];
% 
%         end
%         
%         
%         
%         coR_fr = sum(coRap(:,:,iRa,iRb), 'all')/sum(coRdur(:,:,iRa,iRb), 'all') * 1e3;
%         noR_fr = sum(noRap(:,:,iRa,iRb), 'all')/sum(noRdur(:,:,iRa,iRb), 'all') * 1e3;
%         % coR_fr = sum(coR_replay(1))/sum(coRdur(:)) * 1e3;
%         % noR_fr = sum(noR_replay(1))/sum(noRdur(:)) * 1e3;
% %         X = [4*(iR-1)+1, 4*(iR-1)+2, 4*(iR-1)+3, 4*(iR-1)+4];
%         c = (iRb - 1) * 10 + iRa;
% 
%         X = [1:4];
%         Y = [coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb)), coR_replay_mm(1,iRa,iRb)/sum(coR_replay_mm(:,iRa,iRb)), noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb)), noR_replay_mm(1,iRa,iRb)/sum(noR_replay_mm(:,iRa,iRb))]*100;
%         % Y = [coR_replay(1)/sum(coR_replay)/coR_fr, noR_replay(1)/sum(noR_replay)/noR_fr, coR_replay(1)/sum(coR_replay)/coR_fr, noR_replay(1)/sum(noR_replay)/noR_fr]*100;
% 
%         % Y = [Y, coR_fr, noR_fr];
%         figure(1);
%         subplot(10,10,c)
%         b = bar(X,Y,0.5); hold on;
% %         imagesc(Y); 
%         ax = gca;
% %         ax.XTickLabels = {'coR E & P coF m [%]',  'coR E & P coF mm [%]', 'noR E & P coF m [%]',  'noR E & P coF mm [%]'};
%         ax.XTickLabelRotation = 50;
%         title(sprintf('%s %s', regions{iRa}, regions{iRb}));
%         
%         figure(3)
% %         [chi2stat, pValue(iRa, iRb), expected] = chiSquaredTest([[coR_replay(1,iRa,iRb),  sum(coR_replay(:,iRa,iRb))]; [coR_replay_mm(1,iRa,iRb) , sum(coR_replay_mm(:,iRa,iRb))]]);
%         [chi2stat, pValue(iRa, iRb), expected] = chiSquaredTest([[noR_replay(1,iRa,iRb),  sum(noR_replay(:,iRa,iRb))]; [noR_replay_mm(1,iRa,iRb) , sum(noR_replay_mm(:,iRa,iRb))]]);
%         if pValue(iRa, iRb) < 0.05
% %             if (coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb))) > (coR_replay_mm(1,iRa,iRb)/sum(coR_replay_mm(:,iRa,iRb)))
%             if (noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb))) > (noR_replay_mm(1,iRa,iRb)/sum(noR_replay_mm(:,iRa,iRb)))
%                 plot3([xA(1) xB(1)], [xA(2) xB(2)], [xA(3) xB(3)], 'r-'); hold on;
%             else
%                 plot3([xA(1) xB(1)], [xA(2) xB(2)], [xA(3) xB(3)],  'b-'); hold on;
%             end
%         end
% 
% 
% %         X = [2*(iR-1)+1, 2*(iR-1)+2];
%         X = [1,2];
%         Y = [(coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb))) / (noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb))), (coR_fr / noR_fr)];
% 
%         figure(2);
%         subplot(10,10,c)
%         b = bar(X,Y,0.5); hold on;
% %         title(sprintf('ratio match %s', regionsPlot{iR}));
%         ylim([0 2.5])
% 
%         coR_fr = sum(coRap_mm(:,:,iRa,iRb), 'all')/sum(coRdur_mm(:,:,iRa,iRb), 'all') * 1e3;
%         noR_fr = sum(noRap_mm(:,:,iRa,iRb), 'all')/sum(noRdur_mm(:,:,iRa,iRb), 'all') * 1e3;
%         % coR_fr = sum(coR_replay(1))/sum(coRdur(:)) * 1e3;
%         % noR_fr = sum(noR_replay(1))/sum(noRdur(:)) * 1e3;
%         X = 1:4;
%         Y = [coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb)), noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb)), coR_fr/1e2, noR_fr/1e2]*100;
%         % Y = [coR_replay_mm(1)/sum(coR_replay_mm)/coR_fr, noR_replay_mm(1)/sum(noR_replay_mm)/noR_fr, coR_replay_mm(1)/sum(coR_replay_mm)/coR_fr, noR_replay_mm(1)/sum(noR_replay_mm)/noR_fr]*100;
%         % Y = [Y, coR_fr, noR_fr];
%     %     figure;
%     %     b = bar(X,Y,0.5);
%     %     ax = gca;
%     %     ax.XTickLabels = {'coR E & P coF [%]',  'noR E % P coF [%]', 'CoR CoF Rate [Hz]', 'NoR CoF Rate [Hz]'};
%     %     title(sprintf('E and P mismatch %s', regionsPlot{iR}));
%     % 
%     %     X = 1:2;
%     %     Y = [(coR_replay_mm(1,iR)/sum(coR_replay_mm(:,iR))) / (noR_replay_mm(1,iR)/sum(noR_replay_mm(:,iR))), (coR_fr / noR_fr)];
%     %     figure;
%     %     b = bar(X,Y,0.5);
%     %     title(sprintf('ration mismatch %s', regionsPlot{iR}));
%     %     ylim([0 2.5])
%     end
% 
% 
% end
% 
% %% pairwise combine hem
% close all
% load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
% figure;
% figure;
% figure;
% pValue = [];
% for iRa = 4 %:length(regionsPlot)
%     figure(6)
%     plot3(mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3)), 'ko'); hold on;
%     xA = [mean(regionCoordsAll{iRa}(:,2)),mean(regionCoordsAll{iRa}(:,1)), mean(regionCoordsAll{iRa}(:,3))];
%     
%        
%     
%     
%     for iRb = iRa:length(regionsPlot)
%         
%         xB = [mean(regionCoordsAll{iRb}(:,2)),mean(regionCoordsAll{iRb}(:,1)), mean(regionCoordsAll{iRb}(:,3))];
%        
% %         iRapl = find(contains(regions, regionsPlot{iRa}));
%         iRapl = find(contains(regions, regionsPlot([1:3])));
%         iRbpl = find(contains(regions, regionsPlot([1:3])));
%         coR_replay_pl = zeros(2,1);
%         coR_replay_mm_pl = zeros(2,1);
%         noR_replay_pl = zeros(2,1);
%         noR_replay_mm_pl = zeros(2,1);
%         
%         for Aloop = iRapl
%                 for Bloop = iRbpl
%                     iiA = min([Aloop, Bloop]);
%                     iiB = max([Aloop, Bloop]);
%                     coR_replay_pl = coR_replay_pl + coR_replay(:,iiA,iiB);
%                     coR_replay_mm_pl = coR_replay_mm_pl + coR_replay_mm(:,iiA,iiB);
%                     noR_replay_pl = noR_replay_pl + noR_replay(:,iiA,iiB);
%                     noR_replay_mm_pl = noR_replay_mm_pl + noR_replay_mm(:,iiA,iiB);
%                     
%                 end
%         end
%         coR_fr = sum(coRap(:,:,iRa,iRb), 'all')/sum(coRdur(:,:,iRa,iRb), 'all') * 1e3;
%         noR_fr = sum(noRap(:,:,iRa,iRb), 'all')/sum(noRdur(:,:,iRa,iRb), 'all') * 1e3;
%         % coR_fr = sum(coR_replay(1))/sum(coRdur(:)) * 1e3;
%         % noR_fr = sum(noR_replay(1))/sum(noRdur(:)) * 1e3;
% %         X = [4*(iR-1)+1, 4*(iR-1)+2, 4*(iR-1)+3, 4*(iR-1)+4];
%         c = (iRb - 1) * 10 + iRa;
% 
%         X = [1:4];
%         Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)]*100;
%         % Y = [coR_replay(1)/sum(coR_replay)/coR_fr, noR_replay(1)/sum(noR_replay)/noR_fr, coR_replay(1)/sum(coR_replay)/coR_fr, noR_replay(1)/sum(noR_replay)/noR_fr]*100;
% 
%         % Y = [Y, coR_fr, noR_fr];
%         figure(4);
%         subplot(10,10,c)
%         b = bar(X,Y,0.5); hold on;
% %         imagesc(Y); 
%         ax = gca;
% %         ax.XTickLabels = {'coR E & P coF m [%]',  'coR E & P coF mm [%]', 'noR E & P coF m [%]',  'noR E & P coF mm [%]'};
%         ax.XTickLabelRotation = 50;
% %         title(sprintf('E and P match %s', regionsPlot{iR}));
%         
%         figure(6)
%         [chi2stat, pValue(iRa, iRb), expected] = chiSquaredTest([[coR_replay_pl(1),  sum(coR_replay_pl)]; [coR_replay_mm_pl(1) , sum(coR_replay_mm_pl)]]);
%         [chi2stat, pValue(iRa, iRb), expected] = chiSquaredTest([[noR_replay_pl(1),  sum(noR_replay_pl)]; [noR_replay_mm_pl(1) , sum(noR_replay_mm_pl)]]);
%         if pValue(iRa, iRb) < 0.05
% %             if (coR_replay_pl(1)/sum(coR_replay_pl)) > (coR_replay_mm_pl(1)/sum(coR_replay_mm_pl))
%             if (noR_replay_pl(1)/sum(noR_replay_pl)) > (noR_replay_mm_pl(1)/sum(noR_replay_mm_pl))
%                 plot3([xA(1) xB(1)], [xA(2) xB(2)], [xA(3) xB(3)], 'r-'); hold on;
%             else
%                 plot3([xA(1) xB(1)], [xA(2) xB(2)], [xA(3) xB(3)],  'b-'); hold on;
%             end
%         end
% 
% 
% %         X = [2*(iR-1)+1, 2*(iR-1)+2];
%         X = [1,2];
%         Y = [(coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb))) / (noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb))), (coR_fr / noR_fr)];
% 
%         figure(5);
%         subplot(10,10,c)
%         b = bar(X,Y,0.5); hold on;
% %         title(sprintf('ratio match %s', regionsPlot{iR}));
%         ylim([0 2.5])
% 
%         coR_fr = sum(coRap_mm(:,:,iRa,iRb), 'all')/sum(coRdur_mm(:,:,iRa,iRb), 'all') * 1e3;
%         noR_fr = sum(noRap_mm(:,:,iRa,iRb), 'all')/sum(noRdur_mm(:,:,iRa,iRb), 'all') * 1e3;
%         % coR_fr = sum(coR_replay(1))/sum(coRdur(:)) * 1e3;
%         % noR_fr = sum(noR_replay(1))/sum(noRdur(:)) * 1e3;
%         X = 1:4;
%         Y = [coR_replay(1,iRa,iRb)/sum(coR_replay(:,iRa,iRb)), noR_replay(1,iRa,iRb)/sum(noR_replay(:,iRa,iRb)), coR_fr/1e2, noR_fr/1e2]*100;
%         % Y = [coR_replay_mm(1)/sum(coR_replay_mm)/coR_fr, noR_replay_mm(1)/sum(noR_replay_mm)/noR_fr, coR_replay_mm(1)/sum(coR_replay_mm)/coR_fr, noR_replay_mm(1)/sum(noR_replay_mm)/noR_fr]*100;
%         % Y = [Y, coR_fr, noR_fr];
%     %     figure;
%     %     b = bar(X,Y,0.5);
%     %     ax = gca;
%     %     ax.XTickLabels = {'coR E & P coF [%]',  'noR E % P coF [%]', 'CoR CoF Rate [Hz]', 'NoR CoF Rate [Hz]'};
%     %     title(sprintf('E and P mismatch %s', regionsPlot{iR}));
%     % 
%     %     X = 1:2;
%     %     Y = [(coR_replay_mm(1,iR)/sum(coR_replay_mm(:,iR))) / (noR_replay_mm(1,iR)/sum(noR_replay_mm(:,iR))), (coR_fr / noR_fr)];
%     %     figure;
%     %     b = bar(X,Y,0.5);
%     %     title(sprintf('ration mismatch %s', regionsPlot{iR}));
%     %     ylim([0 2.5])
%     end
% 
% 
% end
% 
% 
% 
 %%
regionColors =  brewermap(12, 'Dark2');
contraColor = brewermap(12, 'Accent');
% close all
load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
pValue = [];

bW = 1.1;
iRa = 4; 
iRapl = find(contains(regions, regionsPlot{iRa}));
iRbpl = find(contains(regions, regionsPlot([1:3])));
coR_replay_pl = zeros(2,1);
coR_replay_mm_pl = zeros(2,1);
coR_replay_shuff = zeros(2,1);
noR_replay_pl = zeros(2,1);
noR_replay_mm_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end

            dat = [sum(coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0); sum(coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0)];
            coR_replay_pl = coR_replay_pl + dat;
            dat = [sum(coRcoF_mm(:,1,iiA, iiB) > 0 & coRcoF_mm(:,2,iiA, iiB) > 0); sum(coRcoF_mm(:,1,iiA, iiB) == 0 | coRcoF_mm(:,2,iiA, iiB) == 0)];
            coR_replay_mm_pl = coR_replay_mm_pl + dat;
            
            dat = [];
            for iter = 1:25
                dat = [dat,[sum(noRcoF(:,1,iiA, iiB, iter) > 0 & noRcoF(:,2,iiA, iiB, iter) > 0); sum(noRcoF(:,1,iiA, iiB,iter) == 0 | noRcoF(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_pl = noR_replay_pl + round(mean(dat,2));
            dat = [];
            for iter = 1:25
                dat = [dat, [sum(noRcoF_mm(:,1,iiA, iiB, iter) > 0 & noRcoF_mm(:,2,iiA, iiB,iter) > 0); sum(noRcoF_mm(:,1,iiA, iiB,iter) == 0 | noRcoF_mm(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_mm_pl = noR_replay_mm_pl + round(mean(dat,2));

            coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];

            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');
            
            
        end
end

coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;

% figure('Position',[1191 647 608 638]);


X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)]*100;

% figure('Position',[1201 947 781 362]);
figure('Position',[1191 647 608 638]);
ha = tight_subplot(2,2,[0.25 ,0.15],[.05 .05],[.15 .15]);
axes(ha(1))
yyaxis left
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
% b = bar_striped(X(2),Y(2),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(iRa,:);
% b.FaceAlpha = 1;
% b.LineWidth = 1;
b = bar(X(3),Y(3),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;
% b = bar_striped(X(4),Y(4),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(iRa,:);
% b.FaceAlpha = 0.5;
% b.LineWidth = 1;

ax = gca;
ax.XTick = 1:4;
ax.XTickLabelRotation = 50;
box off

X = [0.65:3.65;1.35:4.35]';
Y = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3] %1:length(Y)
%     pl = plot([X(iL,1) X(iL,2)], [Y(iL) Y(iL)], '-'); hold on;
%     pl.LineWidth = 2.5;
%     pl.Color = regionColors(3,:);
    
    pl = plot(mean([X(iL,1) X(iL,2)]), Y(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3,:);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;

end
% b = bar(X,Y,0.4); hold on;
ax = gca;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;
box off

ax.YAxis(1).Label.String = '% coF in E & P';
% ax.YAxis(2).Label.String = 'coFire rate [Hz]';
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(2).Limits = [0 max(Y)+0.05*max(Y)];
ax.YAxis(1).FontSize = 11;
ax.YAxis(2).FontSize = 11;
ax.LineWidth = 1;

[chi2stat, pValue(1,1), expected] = chiSquaredTest([[coR_replay_pl(1),  coR_replay_pl(2)]; [coR_replay_mm_pl(1) , coR_replay_mm_pl(2)]]);
[chi2stat, pValue(1,2), expected] = chiSquaredTest([[noR_replay_pl(1),  noR_replay_pl(2)]; [noR_replay_mm_pl(1) , noR_replay_mm_pl(2)]]);
[chi2stat, pValue(1,3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]),  sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]),  sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

iRa = 5; 
iRapl = find(contains(regions, regionsPlot{iRa}));
iRbpl = find(contains(regions, regionsPlot([1:3])));
coR_replay_pl = zeros(2,1);
coR_replay_mm_pl = zeros(2,1);
coR_replay_shuff = zeros(2,1);
noR_replay_pl = zeros(2,1);
noR_replay_mm_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end

            dat = [sum(coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0); sum(coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0)];
            coR_replay_pl = coR_replay_pl + dat;
            dat = [sum(coRcoF_mm(:,1,iiA, iiB) > 0 & coRcoF_mm(:,2,iiA, iiB) > 0); sum(coRcoF_mm(:,1,iiA, iiB) == 0 | coRcoF_mm(:,2,iiA, iiB) == 0)];
            coR_replay_mm_pl = coR_replay_mm_pl + dat;
            dat = [sum(coRcoFShuff(:,1,iiA, iiB) > 0 & coRcoFShuff(:,2,iiA, iiB) > 0); sum(coRcoFShuff(:,1,iiA, iiB) == 0 | coRcoFShuff(:,2,iiA, iiB) == 0)];
            coR_replay_shuff = coR_replay_shuff + dat;
            dat = [];
            for iter = 1:25
                dat = [dat,[sum(noRcoF(:,1,iiA, iiB, iter) > 0 & noRcoF(:,2,iiA, iiB, iter) > 0); sum(noRcoF(:,1,iiA, iiB,iter) == 0 | noRcoF(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_pl = noR_replay_pl + round(mean(dat,2));
            dat = [];
            for iter = 1:25
                dat = [dat, [sum(noRcoF_mm(:,1,iiA, iiB, iter) > 0 & noRcoF_mm(:,2,iiA, iiB,iter) > 0); sum(noRcoF_mm(:,1,iiA, iiB,iter) == 0 | noRcoF_mm(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_mm_pl = noR_replay_mm_pl + round(mean(dat,2));

            coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');
            
          
            

        end
end
    




coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;




X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)]*100;

% figure('Position',[1216 1008 495 565]);
axes(ha(2))
yyaxis left
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
% b = bar_striped(X(2),Y(2),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(iRa,:);
% b.FaceAlpha = 1;
% b.LineWidth = 1;
b = bar(X(3),Y(3),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;
% b = bar_striped(X(4),Y(4),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(iRa,:);
% b.FaceAlpha = 0.5;
% b.LineWidth = 1;
box off

X = [0.65:3.65;1.35:4.35]';
Y = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

yyaxis right
for iL = [1 3] %1:length(Y)
%     pl = plot([X(iL,1) X(iL,2)], [Y(iL) Y(iL)], '-'); hold on;
%     pl.LineWidth = 2.5;
%     pl.Color = regionColors(3,:);

    pl = plot(mean([X(iL,1) X(iL,2)]), Y(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3,:);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;

end
ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

ax.YAxis(1).Label.String = '% coF in E & P';
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 0.25];
ax.YAxis(2).Limits = [0 max(Y)+0.05*max(Y)];
ax.YAxis(1).FontSize = 11; 
ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pValue(2,1), expected] = chiSquaredTest([[coR_replay_pl(1),  coR_replay_pl(2)]; [coR_replay_mm_pl(1) , coR_replay_mm_pl(2)]]);
[chi2stat, pValue(2,2), expected] = chiSquaredTest([[noR_replay_pl(1),  noR_replay_pl(2)]; [noR_replay_mm_pl(1) , noR_replay_mm_pl(2)]]);
[chi2stat, pValue(2,3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]),  sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]),  sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);

iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_pl = zeros(2,1);
coR_replay_mm_pl = zeros(2,1);
coR_replay_shuff = zeros(2,1);
noR_replay_pl = zeros(2,1);
noR_replay_mm_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if ~strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end

            dat = [sum(coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0); sum(coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0)];
            coR_replay_pl = coR_replay_pl + dat;
            dat = [sum(coRcoF_mm(:,1,iiA, iiB) > 0 & coRcoF_mm(:,2,iiA, iiB) > 0); sum(coRcoF_mm(:,1,iiA, iiB) == 0 | coRcoF_mm(:,2,iiA, iiB) == 0)];
            coR_replay_mm_pl = coR_replay_mm_pl + dat;
            dat = [];
            for iter = 1:25
                dat = [dat,[sum(noRcoF(:,1,iiA, iiB, iter) > 0 & noRcoF(:,2,iiA, iiB, iter) > 0); sum(noRcoF(:,1,iiA, iiB,iter) == 0 | noRcoF(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_pl = noR_replay_pl + round(mean(dat,2));
            dat = [];
            for iter = 1:25
                dat = [dat, [sum(noRcoF_mm(:,1,iiA, iiB, iter) > 0 & noRcoF_mm(:,2,iiA, iiB,iter) > 0); sum(noRcoF_mm(:,1,iiA, iiB,iter) == 0 | noRcoF_mm(:,2,iiA, iiB,iter) == 0)]];
            end
            noR_replay_mm_pl = noR_replay_mm_pl + round(mean(dat,2));

            

            coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end
    
coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;




X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)]*100;

% figure('Position',[1216 1008 495 565]);
axes(ha(3))
yyaxis left
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(6,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
% b = bar_striped(X(2),Y(2),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(6,:);
% b.FaceAlpha = 1;
% b.LineWidth = 1;
b = bar(X(3),Y(3),bW); hold on;
b.FaceColor = regionColors(6,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;
% b = bar_striped(X(4),Y(4),bW, 0.02, 5); hold on;
% b.FaceColor = regionColors(6,:);
% b.FaceAlpha = 0.5;
% b.LineWidth = 1;



X = [0.65:3.65;1.35:4.35]';
Y = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3] %1:length(Y)
%     pl = plot([X(iL,1) X(iL,2)], [Y(iL) Y(iL)], '-'); hold on;
%     pl.LineWidth = 2.5;
%     pl.Color = regionColors(3,:);

    pl = plot(mean([X(iL,1) X(iL,2)]), Y(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3,:);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;


end
ax = gca;
ax.XTick = 1:4;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8*regionColors(3,:);
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;
box off

ax.YAxis(1).Label.String = '% coF in E & P';
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(2).Limits = [0 max(Y)+0.05*max(Y)];
ax.YAxis(1).FontSize = 11; 
ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pValue(3,1), expected] = chiSquaredTest([[coR_replay_pl(1),  coR_replay_pl(2)]; [coR_replay_mm_pl(1) , coR_replay_mm_pl(2)]]);
[chi2stat, pValue(3,2), expected] = chiSquaredTest([[noR_replay_pl(1),  noR_replay_pl(2)]; [noR_replay_mm_pl(1) , noR_replay_mm_pl(2)]]);
[chi2stat, pValue(3,3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]),  sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]),  sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);


iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_pl = zeros(2,1);
coR_replay_mm_pl = zeros(2,1);
full_replay_pl = zeros(2,1);
full_replay_mm_pl = zeros(2,1); 
coR_replay_shuff = zeros(2,1);
noR_replay_pl = zeros(2,1);
noR_replay_mm_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
            if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
%             dat = [sum(coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0); sum(coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0)];
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0) & (respTimes(:,iiA,iiB) > median(respTimes, 'all', 'omitnan'));
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) > median(respTimes, 'all', 'omitnan'));
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF(ii1,1,iiA, iiB)) + sum(coRcoF(ii1,2,iiA, iiB)); sum(coRcoF(ii2,1,iiA, iiB)) + sum(coRcoF(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_pl = coR_replay_pl + dat;
            
            ii1 = coRcoF_mm(:,1,iiA, iiB) > 0 & coRcoF_mm(:,2,iiA, iiB) > 0;
            ii2 = coRcoF_mm(:,1,iiA, iiB) == 0 | coRcoF_mm(:,2,iiA, iiB) == 0;
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF_mm(ii1,1,iiA, iiB)) + sum(coRcoF_mm(ii1,2,iiA, iiB)); sum(coRcoF_mm(ii2,1,iiA, iiB)) + sum(coRcoF_mm(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_mm_pl = coR_replay_mm_pl + dat;
            dat = [];
%             
%             dat = [sum(fullCoF(:,1,iiA, iiB) > 0 & fullCoF(:,2,iiA, iiB) > 0); sum(fullCoF(:,1,iiA, iiB) == 0 | fullCoF(:,2,iiA, iiB) == 0)];
%             full_replay_pl = full_replay_pl + dat;
            
%             dat = [sum(fullCoF_mm(:,1,iiA, iiB) > 0 & fullCoF_mm(:,2,iiA, iiB) > 0); sum(fullCoF_mm(:,1,iiA, iiB) == 0 | fullCoF_mm(:,2,iiA, iiB) == 0)];
%             full_replay_mm_pl = full_replay_mm_pl + dat;
            dat = [];
            
            
            for iter = 1:25
                ii1 = (noRcoF(:,1,iiA, iiB, iter) > 0 & noRcoF(:,2,iiA, iiB, iter) > 0) & (respTimes(:,iiA,iiB) > median(respTimes, 'all', 'omitnan'));
                ii2 = (noRcoF(:,1,iiA, iiB, iter) == 0 | noRcoF(:,2,iiA, iiB, iter) == 0) & (respTimes(:,iiA,iiB) > median(respTimes, 'all', 'omitnan'));
%                 dat = [dat,[sum(noRcoF(ii1,1,iiA, iiB, iter)) + sum(noRcoF(ii1,2,iiA, iiB, iter)); sum(noRcoF(ii2,1,iiA, iiB,iter)) +  sum(noRcoF(ii2,2,iiA, iiB,iter))]];
                dat = [dat,[sum(ii1); sum(ii2)]];
            end
            noR_replay_pl = noR_replay_pl + round(mean(dat,2));
            dat = [];
            for iter = 1:25
                ii1 = noRcoF_mm(:,1,iiA, iiB, iter) > 0 & noRcoF_mm(:,2,iiA, iiB, iter) > 0;
                ii2 = noRcoF_mm(:,1,iiA, iiB, iter) == 0 | noRcoF_mm(:,2,iiA, iiB, iter) == 0;
%                 dat = [dat,[sum(noRcoF_mm(ii1,1,iiA, iiB, iter)) + sum(noRcoF_mm(ii1,2,iiA, iiB, iter)); sum(noRcoF_mm(ii2,1,iiA, iiB,iter)) +  sum(noRcoF_mm(ii2,2,iiA, iiB,iter))]];
                dat = [dat,[sum(ii1); sum(ii2)]];

            end
            noR_replay_mm_pl = noR_replay_mm_pl + round(mean(dat,2));

            coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end
    
fprintf('\n')



coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;


% fullCoF

X = [1:4];
Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), noR_replay_pl(1)/sum(noR_replay_pl), noR_replay_mm_pl(1)/sum(noR_replay_mm_pl)]*100;

% figure('Position',[1216 1008 495 565]);
axes(ha(4))
yyaxis left
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = contraColor(7,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
% b = bar_striped(X(2),Y(2),bW, 0.02, 5); hold on;
% b.FaceColor = contraColor(7,:);
% b.FaceAlpha = 1;
% b.LineWidth = 1;
b = bar(X(3),Y(3),bW); hold on;
b.FaceColor = contraColor(7,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;
% b = bar_striped(X(4),Y(4),bW, 0.02, 5); hold on;
% b.FaceColor = contraColor(7,:);
% b.FaceAlpha = 0.5;
% b.LineWidth = 1;

X = [0.65:3.65;1.35:4.35]';
Y = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];
yyaxis right
for iL = [1 3] %1:length(Y)
%     pl = plot([X(iL,1) X(iL,2)], [Y(iL) Y(iL)], '-'); hold on;
%     pl.LineWidth = 2.5;
%     pl.Color = regionColors(3,:);

    pl = plot(mean([X(iL,1) X(iL,2)]), Y(iL), '^'); hold on;
    pl.LineWidth = 1.5;
    pl.MarkerEdgeColor = regionColors(3,:);
    pl.MarkerFaceColor = [0.7 0.7 0.7];
    pl.MarkerSize = 9;


end

ax = gca;
ax.YAxis(1).Label.String = '% coF in E & P';
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(2).Limits = [0 max(Y)+0.05*max(Y)]; 
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).FontSize = 11; 
ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;
[chi2stat1, pValue(4,1), expected1] = chiSquaredTest([[coR_replay_pl(1),  coR_replay_pl(2)]; [coR_replay_mm_pl(1) , coR_replay_mm_pl(2)]]);
[chi2stat2, pValue(4,2), expected2] = chiSquaredTest([[noR_replay_pl(1),  noR_replay_pl(2)]; [noR_replay_mm_pl(1) , noR_replay_mm_pl(2)]]);
[chi2stat, pValue(4,3), expected] = chiSquaredTest([[sum([coR_replay_pl(1), coR_replay_mm_pl(1)]),  sum([coR_replay_pl(2), coR_replay_mm_pl(2)])]; [sum([noR_replay_pl(1), noR_replay_mm_pl(1)]),  sum([noR_replay_pl(2), noR_replay_mm_pl(2)])]]);
% [chi2stat2, pV, expected2] = chiSquaredTest([[full_replay_pl(1),  full_replay_pl(2)]; [full_replay_mm_pl(1) , full_replay_mm_pl(2)]]);

fig = gcf;
fig.Color = 'w'; 
box off

adj_p = nan(size(pValue));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pValue))]=fdr_bh(pValue(~isnan(pValue)),0.05,'pdep','yes');



savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplay_%s.pdf', tag)))





%%


results = compare_chi_squared_effect_sizes(chi2stat1, sum(expected1(:)), 2, 2, chi2stat2, sum(expected2(:)), 2,2)

%% Compare co-ripple and no-ripple


close all
iRa = 5; bW = 0.8;
iRapl = find(contains(regions, regionsPlot{iRa}));
iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_pl = zeros(2,1);
noR_replay_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
            if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            
            dat = [sum(ii1); sum(ii2)];
            coR_replay_pl = coR_replay_pl + dat;
            
            ii1 = (noRcoF(:,1,iiA, iiB) > 0 & noRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (noRcoF(:,1,iiA, iiB) == 0 | noRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            noR_replay_pl = noR_replay_pl + dat;

            
            
            
                

%             coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end

coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;
    
fprintf('\n')

Y = [coR_replay_pl(1)/sum(coR_replay_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), coR_replay_slow_pl(1)/sum(coR_replay_slow_pl), noR_replay_pl(1)/sum(noR_replay_pl)]*100;
X = [1,2,3,3.5];

figure('Position',[1 627 256 292]);
% ha = tight_subplot(1,4,[0.25 ,0.05],[.05 .05],[.075 .05]);
% axes(ha(1))
b = bar(X(1),Y(1),bW); hold on;
% b.FaceColor = regionColors(iRa,:);
b.FaceColor = regionColors(7,:);
b.FaceAlpha = 1;
b.LineWidth = 1;


b = bar(X(4),Y(4),bW); hold on;
% b.FaceColor = regionColors(iRa,:);
b.FaceColor = regionColors(7,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

box off

X = [0.65:3.65;1.35:4.35]';
Yfr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

yyaxis right

pl = plot(1, Yfr(1), '^'); hold on;
pl.LineWidth = 1.5;
pl.MarkerEdgeColor = regionColors(3,:);
pl.MarkerFaceColor = [0.7 0.7 0.7];
pl.MarkerSize = 9;

pl = plot(3.5, Yfr(4), '^'); hold on;
pl.LineWidth = 1.5;
pl.MarkerEdgeColor = regionColors(3,:);
pl.MarkerFaceColor = [0.7 0.7 0.7];
pl.MarkerSize = 9;



ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

ax.YAxis(1).Label.String = '% coF in E & P'; 
ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 max(Y)+0.05*max(Y)];
ax.YAxis(2).Limits = [0 max(Yfr)+0.05*max(Yfr)];
ax.YAxis(1).FontSize = 11; 
ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pV, expected] = chiSquaredTest([[coR_replay_fast_pl(1),  coR_replay_fast_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);

fig = gcf;
fig.Color = 'w'; 
box off


savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplay_cono_contra_%s.pdf', tag)))

%% Compare fast RT to slow RT
% close all
iRa = 4; bW = 0.8;
iRapl = find(contains(regions, regionsPlot{iRa}));
% iRapl = find(contains(regions, regionsPlot([1:3])));
iRbpl = find(contains(regions, regionsPlot([1:3])));
coR_replay_slow_pl = zeros(2,1);
coR_replay_fast_pl = zeros(2,1);
noR_replay_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
%             if ~strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            
            dat = [sum(ii1); sum(ii2)];
            coR_replay_fast_pl = coR_replay_fast_pl + dat;
            
            ii1 = (noRcoF(:,1,iiA, iiB) > 0 & noRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (noRcoF(:,1,iiA, iiB) == 0 | noRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            noR_replay_pl = noR_replay_pl + dat;

            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);     
%             ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
%             ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF(ii1,1,iiA, iiB)) + sum(coRcoF(ii1,2,iiA, iiB)); sum(coRcoF(ii2,1,iiA, iiB)) + sum(coRcoF(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_slow_pl = coR_replay_slow_pl + dat;
            
            
                

%             coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end

coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;
    
fprintf('\n')

Y = [coR_replay_fast_pl(1)/sum(coR_replay_fast_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), coR_replay_slow_pl(1)/sum(coR_replay_slow_pl), noR_replay_pl(1)/sum(noR_replay_pl)]*100;
X = [1,2,3,3.5];

figure('Position',[2 632 755 299]);
ha = tight_subplot(1,4,[0.25 ,0.05],[.05 .05],[.075 .05]);
axes(ha(1))
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
b = bar_striped(X(2),Y(3),bW, 0.02, 5); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(4),Y(4),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

box off

X = [0.65:3.65;1.35:4.35]';
Yfr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

% yyaxis right
% 
% pl = plot(mean([X(1:2,1) X(1:2,2)],'all'), Yfr(1), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;
% 
% pl = plot(3.5, Yfr(4), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;



ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

%ax.YAxis(1).Label.String = '% coF in E & P'; 
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
% ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 max(Y)+0.05*max(Y)];
% ax.YAxis(2).Limits = [0 max(Yfr)+0.05*max(Yfr)];
ax.YAxis(1).FontSize = 11; 
% ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pV, expected] = chiSquaredTest([[coR_replay_fast_pl(1),  coR_replay_fast_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);
% [chi2stat, pV, expected] = chiSquaredTest([[noR_replay_pl(1),  noR_replay_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);

iRa = 5; bW = 0.8;
iRapl = find(contains(regions, regionsPlot{iRa}));
% iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:3])));
coR_replay_slow_pl = zeros(2,1);
coR_replay_fast_pl = zeros(2,1);
noR_replay_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
%             if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            coR_replay_fast_pl = coR_replay_fast_pl + dat;
            
            ii1 = (noRcoF(:,1,iiA, iiB) > 0 & noRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (noRcoF(:,1,iiA, iiB) == 0 | noRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            noR_replay_pl = noR_replay_pl + dat;

            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF(ii1,1,iiA, iiB)) + sum(coRcoF(ii1,2,iiA, iiB)); sum(coRcoF(ii2,1,iiA, iiB)) + sum(coRcoF(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_slow_pl = coR_replay_slow_pl + dat;
            
            
                

%             coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end
  coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3; 
fprintf('\n')

Y = [coR_replay_fast_pl(1)/sum(coR_replay_fast_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), coR_replay_slow_pl(1)/sum(coR_replay_slow_pl), noR_replay_pl(1)/sum(noR_replay_pl)]*100;
X = [1,2,3,3.5];

axes(ha(2))
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
b = bar_striped(X(2),Y(3),bW, 0.02, 5); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(4),Y(4),bW); hold on;
b.FaceColor = regionColors(iRa,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

box off

X = [0.65:3.65;1.35:4.35]';
Yfr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

% yyaxis right
% 
% pl = plot(mean([X(1:2,1) X(1:2,2)],'all'), Yfr(1), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;
% 
% pl = plot(3.5, Yfr(4), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;



ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

%ax.YAxis(1).Label.String = '% coF in E & P'; 
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
% ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 max(Y)+0.05*max(Y)];
% ax.YAxis(2).Limits = [0 max(Yfr)+0.05*max(Yfr)];
ax.YAxis(1).FontSize = 11; 
% ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pV, expected] = chiSquaredTest([[coR_replay_fast_pl(1),  coR_replay_fast_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);

iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_slow_pl = zeros(2,1);
coR_replay_fast_pl = zeros(2,1);
noR_replay_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
            if ~strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            coR_replay_fast_pl = coR_replay_fast_pl + dat;
            
            ii1 = (noRcoF(:,1,iiA, iiB) > 0 & noRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (noRcoF(:,1,iiA, iiB) == 0 | noRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            noR_replay_pl = noR_replay_pl + dat;

            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF(ii1,1,iiA, iiB)) + sum(coRcoF(ii1,2,iiA, iiB)); sum(coRcoF(ii2,1,iiA, iiB)) + sum(coRcoF(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_slow_pl = coR_replay_slow_pl + dat;
            
            
                

%             coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end
  coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3; 
fprintf('\n')

Y = [coR_replay_fast_pl(1)/sum(coR_replay_fast_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), coR_replay_slow_pl(1)/sum(coR_replay_slow_pl), noR_replay_pl(1)/sum(noR_replay_pl)]*100;
X = [1,2,3,3.5];

axes(ha(3))
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(6,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
b = bar_striped(X(2),Y(3),bW, 0.02, 5); hold on;
b.FaceColor = regionColors(6,:);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(4),Y(4),bW); hold on;
b.FaceColor = regionColors(6,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

box off

X = [0.65:3.65;1.35:4.35]';
Yfr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

% yyaxis right
% 
% pl = plot(mean([X(1:2,1) X(1:2,2)],'all'), Yfr(1), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;
% 
% pl = plot(3.5, Yfr(4), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;



ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

%ax.YAxis(1).Label.String = '% coF in E & P'; 
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
% ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 max(Y)+0.05*max(Y)];
% ax.YAxis(2).Limits = [0 max(Yfr)+0.05*max(Yfr)];
ax.YAxis(1).FontSize = 11; 
% ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pV, expected] = chiSquaredTest([[coR_replay_fast_pl(1),  coR_replay_fast_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);

iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_slow_pl = zeros(2,1);
coR_replay_fast_pl = zeros(2,1);
noR_replay_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
            if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})
                
            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) <= median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            coR_replay_fast_pl = coR_replay_fast_pl + dat;
            
            ii1 = (noRcoF(:,1,iiA, iiB) > 0 & noRcoF(:,2,iiA, iiB) > 0)   & (loads(:,iiA,iiB) == 3);
            ii2 = (noRcoF(:,1,iiA, iiB) == 0 | noRcoF(:,2,iiA, iiB) == 0) & (loads(:,iiA,iiB) == 3);
            dat = [sum(ii1); sum(ii2)];
            noR_replay_pl = noR_replay_pl + dat;

            ii1 = (coRcoF(:,1,iiA, iiB) > 0 & coRcoF(:,2,iiA, iiB) > 0)   & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
            ii2 = (coRcoF(:,1,iiA, iiB) == 0 | coRcoF(:,2,iiA, iiB) == 0) & (respTimes(:,iiA,iiB) > median(respTimes(loads ==3), 'all', 'omitnan')) & (loads(:,iiA,iiB) == 3);
%             dat = [sum(ii1); sum(ii2)];
%             dat = [sum(coRcoF(ii1,1,iiA, iiB)) + sum(coRcoF(ii1,2,iiA, iiB)); sum(coRcoF(ii2,1,iiA, iiB)) + sum(coRcoF(ii2,2,iiA, iiB))];
            dat = [sum(ii1); sum(ii2)];
            coR_replay_slow_pl = coR_replay_slow_pl + dat;
            
            
                

%             coR_all = [coR_all; coRcoF(~isnan(coRcoF(:,1,iiA, iiB)),:,iiA, iiB); coRcoF_mm(~isnan(coRcoF_mm(:,1,iiA, iiB)),:,iiA, iiB)];
            
            coRap_pl = coRap_pl + sum(coRap(:,:,iiA, iiB), 'all');
            coRap_mm_pl = coRap_mm_pl + sum(coRap_mm(:,:,iiA, iiB), 'all');
            noRap_pl = noRap_pl + sum(noRap(:,:,iiA, iiB), 'all');
            noRap_mm_pl = noRap_mm_pl + sum(noRap_mm(:,:,iiA, iiB), 'all');
            
            coRdur_pl = coRdur_pl + sum(coRdur(:,:,iiA, iiB), 'all');
            coRdur_mm_pl = coRdur_mm_pl + sum(coRdur_mm(:,:,iiA, iiB), 'all');
            noRdur_pl = noRdur_pl + sum(noRdur(:,:,iiA, iiB), 'all');
            noRdur_mm_pl = noRdur_mm_pl + sum(noRdur_mm(:,:,iiA, iiB), 'all');


        end
end
  coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3; 
fprintf('\n')

Y = [coR_replay_fast_pl(1)/sum(coR_replay_fast_pl), coR_replay_mm_pl(1)/sum(coR_replay_mm_pl), coR_replay_slow_pl(1)/sum(coR_replay_slow_pl), noR_replay_pl(1)/sum(noR_replay_pl)]*100;
X = [1,2,3,3.5];

axes(ha(4))
b = bar(X(1),Y(1),bW); hold on;
b.FaceColor = regionColors(7,:);
b.FaceAlpha = 1;
b.LineWidth = 1;
b = bar_striped(X(2),Y(3),bW, 0.02, 5); hold on;
b.FaceColor = regionColors(7,:);
b.FaceAlpha = 1;
b.LineWidth = 1;

b = bar(X(4),Y(4),bW); hold on;
b.FaceColor = regionColors(7,:);
b.FaceAlpha = 0.5;
b.LineWidth = 1;

box off

X = [0.65:3.65;1.35:4.35]';
Yfr = [coR_fr, coR_mm_fr, noR_fr, noR_mm_fr];

% yyaxis right
% 
% pl = plot(mean([X(1:2,1) X(1:2,2)],'all'), Yfr(1), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;
% 
% pl = plot(3.5, Yfr(4), '^'); hold on;
% pl.LineWidth = 1.5;
% pl.MarkerEdgeColor = regionColors(3,:);
% pl.MarkerFaceColor = [0.7 0.7 0.7];
% pl.MarkerSize = 9;



ax = gca;
ax.XTick = 1:4;
% ax.XTickLabels = {'coR m',  'coR mm', 'noR m',  'noR mm'};
ax.XTickLabelRotation = 50;

%ax.YAxis(1).Label.String = '% coF in E & P'; 
%ax.YAxis(2).Label.String = 'coFire rate [Hz]'; 
ax.YAxis(1).Color = [0 0 0];
% ax.YAxis(2).Color = 0.8*regionColors(3,:);
ax.YAxis(1).Limits = [0 max(Y)+0.05*max(Y)];
% ax.YAxis(2).Limits = [0 max(Yfr)+0.05*max(Yfr)];
ax.YAxis(1).FontSize = 11; 
% ax.YAxis(2).FontSize = 11; 
ax.LineWidth = 1;

[chi2stat, pV, expected] = chiSquaredTest([[coR_replay_fast_pl(1),  coR_replay_fast_pl(2)]; [coR_replay_slow_pl(1),  coR_replay_slow_pl(2)]]);

fig = gcf;
fig.Color = 'w'; 
box off


savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplay_RT_%s.pdf', tag)))

%%
iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));
coR_replay_pl = zeros(2,1);
coR_replay_mm_pl = zeros(2,1);
full_replay_pl = zeros(2,1);
full_replay_mm_pl = zeros(2,1); 
coR_replay_shuff = zeros(2,1);
noR_replay_pl = zeros(2,1);
noR_replay_mm_pl = zeros(2,1);

coRap_pl = 0;
coRap_mm_pl = 0;
noRap_pl = 0;
noRap_mm_pl = 0;
coRdur_pl = 0;
coRdur_mm_pl = 0;
noRdur_pl = 0;
noRdur_mm_pl = 0;
coR_all = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(coR_replay(:,iiA,iiB)) == 0; continue; end
%             if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            fprintf('%s %s \n', regions{iiA}, regions{iiB})

            ii = fullCoF(:,1,iiA, iiB) > 0 & fullCoF(:,2,iiA, iiB) > 0;
            dat = sum(fullCoF(ii,1,iiA, iiB))+sum(fullCoF(ii,2,iiA, iiB));
            full_replay_pl = full_replay_pl/sum(~isnan(fullCoF(ii,1,iiA, iiB))) + dat;
            ii = fullCoF_mm(:,1,iiA, iiB) > 0 & fullCoF_mm(:,2,iiA, iiB) > 0;
            dat = sum(fullCoF_mm(ii,1,iiA, iiB))+sum(fullCoF_mm(ii,2,iiA, iiB));
            full_replay_mm_pl = full_replay_mm_pl/sum(~isnan(fullCoF_mm(ii,1,iiA, iiB))) + dat;
            dat = [];
            
         


        end
end
    
fprintf('\n')



coR_fr = coRap_pl/coRdur_pl * 1e3;
coR_mm_fr = coRap_mm_pl/coRdur_mm_pl * 1e3;

noR_fr = noRap_pl/noRdur_pl * 1e3;
noR_mm_fr = noRap_mm_pl/noRdur_mm_pl * 1e3;








