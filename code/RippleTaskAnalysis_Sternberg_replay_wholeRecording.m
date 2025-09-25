

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

    
    load_trials_all_unitPair = cell(length(regions), length(regions));
    unitPairID = cell(length(regions), length(regions));
    testCoFire = cell(length(regions), length(regions));
    UnitPairType = cell(length(regions), length(regions));

    trialsProbe = zeros(length(regions),length(regions),1e3);

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
%     full_replay = zeros(2,length(regions),length(regions));
%     full_replay_mismatch = zeros(2,5);
    
    coRdur = zeros(2,2,length(regions),length(regions));
    noRdur = zeros(2,2,length(regions),length(regions));
    coRap = zeros(2,2,length(regions),length(regions));
    noRap = zeros(2,2,length(regions),length(regions));
    coRcoF = nan(5e5, 2, length(regions),length(regions));
    noRcoF = nan(5e5, 2, length(regions),length(regions), 25);
    coRcoFShuff = nan(5e5, 2, length(regions),length(regions));
    
    coR_replay_mm = zeros(2,length(regions),length(regions));
    noR_replay_mm = zeros(2,length(regions),length(regions));
    
    coR_replay_PRTHe = nan(5e5, length(regions),length(regions));
    coR_replay_PRTHp = nan(5e5, length(regions),length(regions));
    coR_replay_PRTHe_shuff = nan(5e5, length(regions),length(regions));
    coR_replay_PRTHp_shuff = nan(5e5, length(regions),length(regions));
    rA_replay_PRTHe = nan(5e5, length(regions),length(regions));
    rA_replay_PRTHp = nan(5e5, length(regions),length(regions));
    
    countPRTHp = ones(length(regions),length(regions));
    countPRTHe = ones(length(regions),length(regions));
    countCoRp = ones(length(regions),length(regions));
    countCoRe = ones(length(regions),length(regions));
    
    countPRTHp_shuff = ones(length(regions),length(regions));
    countPRTHe_shuff = ones(length(regions),length(regions));
    countCoRp_shuff = ones(length(regions),length(regions));
    countCoRe_shuff = ones(length(regions),length(regions));
    
    countPRTHpA = ones(length(regions),length(regions));
    countPRTHeA = ones(length(regions),length(regions));
    countRAe = ones(length(regions),length(regions));
    countRAp = ones(length(regions),length(regions));

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
        for iRseed = 1 %1:length(regionsPlot)
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
                                    datAripP = datArip(probeTime+1:probeTime+1e3) > 0  & datBrip(probeTime+1:probeTime+1e3) == 0;
                                    datABnoripP = datArip(probeTime+1:probeTime+1e3) == 0 & datBrip(probeTime+1:probeTime+1e3) == 0;
                                    if trials.loads(iT) == 3                                                      
                                        iE = find(encodIm(:,iT) == probeIm(iT));
                                        iEmm = find(encodIm(:,iT) ~= probeIm(iT));
                                        iEmm = iEmm(randi(length(iEmm)));
        %                                 spikeMaskE = spikeMask(:,imageTime(iE):imageTime(iE+1));
                                        spikeMaskE = spikeMask(:,imageTime(iE):imageTime(iE)+2e3);
                                        spikeMaskEmm = spikeMask(:,imageTime(iEmm):imageTime(iEmm)+2e3);
                                        datABcoripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & datBrip(imageTime(iE):imageTime(iE)+2e3);
                                        datAripE = datArip(imageTime(iE):imageTime(iE)+2e3) > 0 & datBrip(imageTime(iE):imageTime(iE)+2e3) == 0;
                                        datABnoripE = datArip(imageTime(iE):imageTime(iE)+2e3) == 0 & datBrip(imageTime(iE):imageTime(iE)+2e3) == 0;
                                    else
                                        spikeMaskE = spikeMask(:,imageTime(1):imageTime(1)+2e3);
                                        datABcoripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & datBrip(imageTime(1):imageTime(1)+2e3) > 0;
                                        datAripE = datArip(imageTime(1):imageTime(1)+2e3) > 0 & datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                                        datABnoripE = datArip(imageTime(1):imageTime(1)+2e3) == 0 & datBrip(imageTime(1):imageTime(1)+2e3) == 0;
                                    end

%                                     if sum(datABcoripP) == 0 || sum(datABcoripE) == 0; continue; end
%                                     if sum(datABnoripP) == 0 || sum(datABnoripE) == 0; continue; end

                                    %co-ripples
                                    datAP =  find(spikeMaskP(uRegionA(uA),:)); datBP =  find(spikeMaskP(uRegionB(uB),:));                                  
                                    sP = datAP - datBP'; 
                                    [coFpB, coFpA] = find(sP >= -coFireWindow*2 & sP <= coFireWindow*2);
                                    
                                    
                                    datAE =  find(spikeMaskE(uRegionA(uA),:)); datBE =  find(spikeMaskE(uRegionB(uB),:)); 
                                    sE = datAE - datBE';
                                    [coFeB, coFeA] = find(sE >= -coFireWindow*2 & sE <= coFireWindow*2);
                                   
                                    if ~isempty(coFeA) > 0 && ~isempty(coFpA)
                                        bCoRipE = mean(mask2bounds(datABcoripE),2)';
                                        bCoRipP = mean(mask2bounds(datABcoripP),2)';
                                        bRipAE = mean(mask2bounds(datAripE),2)';
                                        bRipAP = mean(mask2bounds(datAripP),2)';
                                        
                                        mat = nan(2,length(coFpA)); mat(1,:) = datAP(coFpA); mat(2,:) = datBP(coFpB);
                                        coFp = mean(mat, 1);
                                        mat = nan(2,length(coFeA)); mat(1,:) = datAE(coFeA); mat(2,:) = datBE(coFeB);
                                        coFe = mean(mat, 1);
                                        
                                        cP = bCoRipP-coFp';
                                        cE = bCoRipE-coFe';
                                        
                                        coR_replay_PRTHe(countPRTHe(iRa,iRb):countPRTHe(iRa,iRb)+length(cE(:))-1, iRa,iRb) = cE(:);
                                        countPRTHe(iRa,iRb) = countPRTHe(iRa,iRb) + length(cE(:));
                                        countCoRe(iRa,iRb) = countCoRe(iRa,iRb) + length(bCoRipE);

                                        coR_replay_PRTHp(countPRTHp(iRa,iRb):countPRTHp(iRa,iRb)+length(cP(:))-1, iRa,iRb) = cP(:);
                                        countPRTHp(iRa,iRb) = countPRTHp(iRa,iRb) + length(cP(:));
                                        countCoRp(iRa,iRb) = countCoRp(iRa,iRb) + length(bCoRipP);

                                       
                                        
                                        cP = bRipAP-coFp';
                                        cE = bRipAE-coFe';
                                        
                                        rA_replay_PRTHe(countPRTHeA(iRa,iRb):countPRTHeA(iRa,iRb)+length(cE(:))-1, iRa,iRb) = cE(:);
                                        countPRTHeA(iRa,iRb) = countPRTHeA(iRa,iRb) + length(cE(:));
                                        countRAe(iRa,iRb) = countRAe(iRa,iRb) + length(bRipAE);

                                        rA_replay_PRTHp(countPRTHpA(iRa,iRb):countPRTHpA(iRa,iRb)+length(cP(:))-1, iRa,iRb) = cP(:);
                                        countPRTHpA(iRa,iRb) = countPRTHpA(iRa,iRb) + length(cP(:));
                                        countRAp(iRa,iRb) = countRAp(iRa,iRb) + length(bRipAP);

                                    end
                                    
                                                                   
%                                     %co-ripple shuffle
                                    if ~isempty(datAP) && ~isempty(datAE) && ~isempty(datBP) && ~isempty(datBE)
                                        datAP = randi(length(spikeMaskP), [1,length(datAP)]); datBP = randi(length(spikeMaskP), [1,length(datBP)]);                                
                                        sP = datAP - datBP'; 
                                        [coFpB, coFpA] = find(sP >= -coFireWindow*2 & sP <= coFireWindow*2);


                                        datAE = randi(length(spikeMaskE), [1,length(datAE)]); datBE = randi(length(spikeMaskE), [1,length(datBE)]);                                  
                                        sE = datAE - datBE';
                                        [coFeB, coFeA] = find(sE >= -coFireWindow*2 & sE <= coFireWindow*2);

                                        if ~isempty(coFeA) && ~isempty(coFpA)
                                            bCoRipE = mean(mask2bounds(datABcoripE),2)';
                                            bCoRipP = mean(mask2bounds(datABcoripP),2)';
                                            mat = nan(2,length(coFpA)); mat(1,:) = datAP(coFpA); mat(2,:) = datBP(coFpB);
                                            coFp = mean(mat, 1);
                                            mat = nan(2,length(coFeA)); mat(1,:) = datAE(coFeA); mat(2,:) = datBE(coFeB);
                                            coFe = mean(mat, 1);

                                            cP = bCoRipP-coFp';
                                            cE = bCoRipE-coFe';

                                            coR_replay_PRTHe_shuff(countPRTHe_shuff(iRa,iRb):countPRTHe_shuff(iRa,iRb)+length(cE(:))-1, iRa,iRb) = cE(:);
                                            countPRTHe_shuff(iRa,iRb) = countPRTHe_shuff(iRa,iRb) + length(cE(:));
                                            countCoRe_shuff(iRa,iRb) = countCoRe_shuff(iRa,iRb) + length(bCoRipE);

                                            coR_replay_PRTHp_shuff(countPRTHp_shuff(iRa,iRb):countPRTHp_shuff(iRa,iRb)+length(cP(:))-1, iRa,iRb) = cP(:);
                                            countPRTHp_shuff(iRa,iRb) = countPRTHp_shuff(iRa,iRb) + length(cP(:));
                                            countCoRp_shuff(iRa,iRb) = countCoRp_shuff(iRa,iRb) + length(bCoRipP);

                                        end
                                    end
                                end
                            end

                        else
%                             for uA = 1:length(uRegionA)
%                                 for uB = 1:length(uRegionB)
%                                     typeA = units{uRegionA(uA),3};
%                                     typeB = units{uRegionB(uB),3};
%                                     spikeMaskP = spikeMask(:,probeTime+1:probeTime+1e3);
%                                     datABcoripP = datArip(probeTime+1:probeTime+1e3) > 0 & datBrip(probeTime+1:probeTime+1e3) > 0;
%                                     datABnoripP = datArip(probeTime+1:probeTime+1e3) == 0 & datBrip(probeTime+1:probeTime+1e3) == 0;
%                                     if trials.loads(iT) == 3                                                      
%                                         iEmm = find(encodIm(:,iT) ~= probeIm(iT));
%                                         iEmm = iEmm(randi(length(iEmm)));
%                                         spikeMaskEmm = spikeMask(:,imageTime(iEmm):imageTime(iEmm)+2e3);
%                                         datABcoripEmm = datArip(imageTime(iEmm):imageTime(iEmm)+2e3) > 0 & datBrip(imageTime(iEmm):imageTime(iEmm)+2e3);
%                                         datABnoripEmm = datArip(imageTime(iEmm):imageTime(iEmm)+2e3) == 0 & datBrip(imageTime(iEmm):imageTime(iEmm)+2e3) == 0;
%                                     else
%                                         spikeMaskEmm  = spikeMask(:,imageTime(1):imageTime(1)+2e3);
%                                         datABcoripEmm = datArip(imageTime(1):imageTime(1)+2e3) > 0 & datBrip(imageTime(1):imageTime(1)+2e3) > 0;
%                                         datABnoripEmm = datArip(imageTime(1):imageTime(1)+2e3) == 0 & datBrip(imageTime(1):imageTime(1)+2e3) == 0;
%                                     end
% 
%                                     if sum(datABcoripP) == 0 || sum(datABcoripE) == 0; continue; end
% %                                     if sum(datABnoripP) == 0 || sum(datABnoripE) == 0; continue; end
% 
%                                     %co-ripples
%                                     datA =  find(spikeMaskP(uRegionA(uA),datABcoripP)); datB =  find(spikeMaskP(uRegionB(uB),datABcoripP));
%                                     sP = datA - datB'; 
%                                     coFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);
% 
%                                     datA =  find(spikeMaskEmm(uRegionA(uA),datABcoripEmm)); datB =  find(spikeMaskEmm(uRegionB(uB),datABcoripEmm));
%                                     sE = datA - datB';
%                                     coFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);
%                                     
%                                     ii = find(isnan(coRcoF_mm(:,1,iRa, iRb)), 1, 'first');
%                                     coRcoF_mm(ii,1,iRa, iRb) = coFp;
%                                     coRcoF_mm(ii,2,iRa, iRb) = coFe;
%                                     if coFp > 0 && coFe > 0
%                                         coR_replay_mm(1,iRa, iRb) = coR_replay_mm(1,iRa, iRb) + 1;
%                                         coRdur_mm(1,1,iRa, iRb) = coRdur_mm(1,1,iRa, iRb) + sum(datABcoripEmm);
%                                         coRdur_mm(1,2,iRa, iRb) = coRdur_mm(1,2,iRa, iRb) + sum(datABcoripP);
%                                         coRap_mm(1,1,iRa, iRb) = coRap_mm(1,1,iRa, iRb) + sum(coFe);
%                                         coRap_mm(1,2,iRa, iRb) = coRap_mm(1,2,iRa, iRb) + sum(coFp);
% 
%                                     else
%                                         coR_replay_mm(2,iRa, iRb) = coR_replay_mm(2,iRa, iRb) + 1;
%                                         coRdur_mm(2,1,iRa, iRb) = coRdur_mm(2,1,iRa, iRb) + sum(datABcoripEmm);
%                                         coRdur_mm(2,2,iRa, iRb) = coRdur_mm(2,2,iRa, iRb) + sum(datABcoripP);
%                                         coRap_mm(2,1,iRa, iRb) = coRap_mm(2,1,iRa, iRb) + sum(coFe);
%                                         coRap_mm(2,2,iRa, iRb) = coRap_mm(2,2,iRa, iRb) + sum(coFp);
%                                     end
% 
% 
%         %                             end
% 
%                                     %no-ripples
%                                     for iter = 1:25
%                                         ripDur = sum(datABcoripP);
%                                         b = mask2bounds(datABnoripP);
%                                         bDur = b(:,2) - b(:,1);
%                                         if isempty(b) 
%                                             b = [1,1]; ripDur = 0;
%                                         elseif sum(bDur >= ripDur) == 0
%                                             b = b(bDur == max(bDur),:);
%                                             ripDur = max(bDur);
%                                         elseif sum(bDur >= ripDur) > 0
%                                             b(bDur < ripDur,:) = [];
%                                             b = b(randi(size(b,1),1),:);
%                                         end
%                                         
%                                         datABnoripPiter = datABnoripP;
%                                         b = b(randi(size(b,1),1),:);
%                                         randStart = randi([b(1), b(2) - ripDur]);
%                                         datABnoripPiter(1:randStart) = false;
%                                         datABnoripPiter(randStart+ripDur+1:end) = false;
%                                         datA =  find(spikeMaskP(uRegionA(uA),datABnoripPiter)); datB =  find(spikeMaskP(uRegionB(uB),datABnoripPiter));
%                                         sP = datA - datB'; 
%                                         noFp = sum(sP(:) >= -coFireWindow*2 & sP(:) <= coFireWindow*2);
% 
%                                         ripDur = sum(datABcoripEmm);
%                                         b = mask2bounds(datABnoripEmm);
%                                         bDur = b(:,2) - b(:,1);
%                                         if isempty(b) 
%                                             b = [1,1]; ripDur = 0;
%                                         elseif sum(bDur >= ripDur) == 0
%                                             b = b(bDur == max(bDur),:);
%                                             ripDur = max(bDur);
%                                         elseif sum(bDur >= ripDur) > 0
%                                             b(bDur < ripDur,:) = [];
%                                             b = b(randi(size(b,1),1),:);
%                                         end
% 
%                                         datABnoripEmmiter = datABnoripEmm;
%                                         randStart = randi([b(1), b(2) - ripDur]);
%                                         datABnoripEmmiter(1:randStart) = false;
%                                         datABnoripEmmiter(randStart+ripDur+1:end) = false;
%                                         datA =  find(spikeMaskEmm(uRegionA(uA),datABnoripEmmiter)); datB =  find(spikeMaskEmm(uRegionB(uB),datABnoripEmmiter));
%                                         sE = datA - datB';
%                                         noFe = sum(sE(:) >= -coFireWindow*2 & sE(:) <= coFireWindow*2);
% 
%                                         ii = noRcoF_mmCount(iRa, iRb,iter);
%                                         noRcoF_mm(ii,1,iRa, iRb, iter) = noFp;
%                                         noRcoF_mm(ii,2,iRa, iRb, iter) = noFe;
%                                         noRcoF_mmCount(iRa, iRb,iter) = noRcoF_mmCount(iRa, iRb,iter) +1;
%                                         if noFp > 0 && noFe > 0
%                                             noR_replay_mm(1,iRa, iRb) = noR_replay_mm(1,iRa, iRb) + 1;
%                                             noRdur_mm(1,1,iRa, iRb) = noRdur_mm(1,1,iRa, iRb) + sum(datABnoripEmmiter);
%                                             noRdur_mm(1,2,iRa, iRb) = noRdur_mm(1,2,iRa, iRb) + sum(datABnoripPiter);
%                                             noRap_mm(1,1,iRa, iRb) = noRap_mm(1,1,iRa, iRb) + sum(noFe);
%                                             noRap_mm(1,2,iRa, iRb) = noRap_mm(1,2,iRa, iRb) + sum(noFp);
%                                         else
%                                             noRdur_mm(2,1,iRa, iRb) = noRdur_mm(2,1,iRa, iRb) + sum(datABnoripEmmiter);
%                                             noRdur_mm(2,2,iRa, iRb) = noRdur_mm(2,2,iRa, iRb) + sum(datABnoripPiter);
%                                             noRap_mm(2,1,iRa, iRb) = noRap_mm(2,1,iRa, iRb) + sum(noFe);
%                                             noRap_mm(2,2,iRa, iRb) = noRap_mm(2,2,iRa, iRb) + sum(noFp);
%                                         end
%                                     end
%                                     if noFp > 0 && noFe > 0
%                                         noR_replay_mm(1,iRa, iRb) = noR_replay_mm(1,iRa, iRb) + 1;
%                                     else
%                                         noR_replay_mm(2,iRa, iRb) = noR_replay_mm(2,iRa, iRb) + 1;
%                                     end                             
%                                 end %uB
%                             end %uA
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

% 
 %%
regionColors =  brewermap(12, 'Dark2');
contraColor = brewermap(12, 'Accent');
close all
load('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/regionLocations.mat')
pValue = [];
win = 130;
winPlot = 100;
smWin = 11;
binSz = 1;
iRa = 4; 
iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));

close all

% figure('Position',[1133 782 527 105]);
% figure('Position',[1133 782 527 105]);
% ha1 = tight_subplot(1,2,[0.2 ,0.2],[.05 .05],[.15 .15]);

figure('Position',[1324 421 364 136]);
figure('Position',[1324 421 364 136]);
% ha2 = tight_subplot(1,2,[0.2 ,0.2],[.05 .05],[.15 .15]);





iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));

Np = [];
Ne = [];
NpA = [];
NeA = [];
NpShuff = [];
NeShuff = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(countPRTHp(iiA,iiB)) == 1; continue; end
            if ~strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end

            
            times = coR_replay_PRTHe(:,iiA,iiB);
            temp = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRe(iiA,iiB), 'gaussian', smWin);
            
            times = rA_replay_PRTHe(:,iiA,iiB);
            tempA = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countRAe(iiA,iiB), 'gaussian', smWin);
            
            times = coR_replay_PRTHe_shuff(:,iiA,iiB);
            tempShuff = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRe_shuff(iiA,iiB), 'gaussian', smWin);
                        
            z = zscore([temp, tempA, tempShuff]);
            Ne = [Ne; z(1:length(z)/3)];
            NeA = [NeA; z((length(z)/3)+1:(2*length(z)/3))];
            NeShuff = [NeShuff;  z((2*length(z)/3)+1:end)];

            times = coR_replay_PRTHp(:,iiA,iiB);
            temp = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRp(iiA,iiB), 'gaussian', smWin);

            times = rA_replay_PRTHp(:,iiA,iiB);
            tempA = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countRAp(iiA,iiB), 'gaussian', smWin);
            
            times = coR_replay_PRTHp_shuff(:,iiA,iiB);
            tempShuff = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRp_shuff(iiA,iiB), 'gaussian', smWin);
                  
            z = zscore([temp, tempA, tempShuff]);
            Np = [Np; z(1:length(z)/3)];
            NpA = [NpA; z((length(z)/3)+1:(2*length(z)/3))];
            NpShuff = [NpShuff;  z((2*length(z)/3)+1:end)];


        end
end


% figure('Position',[1216 1008 495 565]);
% axes(ha1(1))
figure(1)
subplot(1,2,1)
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255); hold on;
f.LineStyle = '--';
f.LineWidth = 1;
[bl1,bf1] = boundedline(-win:binSz:win, mean(Ne), std(Ne)/sqrt(size(Ne,1)), 'b'); hold on; 
bl1.Color = regionColors(6,:);
bl1.LineWidth = 1.5;
bf1.FaceColor = regionColors(6,:);
bf1.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NeA), std(NeA)/sqrt(size(NeA,1))); hold on; 
bl3.Color = 0.8*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NeShuff), std(NeShuff)/sqrt(size(NeShuff,1))); hold on; 
bl3.Color = 0.5*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;
xlim([-winPlot winPlot])
vl = vline(0); vl.LineWidth =1.5; 
times = -win:binSz:win;
box off

xx = [];
xx(:,1) = mean(Ne(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NeA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NeShuff(:, times >= -25 & times <= 25), 2);
p  = anova1(xx);
fprintf('encode ipsi p: %.4f\n', p)
ax = gca;
ax.LineWidth = 1;
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode ipsi coR vs rA p: %.4f\n', p_perm)
% xx = mean(Ne(:, times >= -25 & times <= 25), 2);
% yy = mean(NeShuff(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode ipsi coR vs shuff p: %.4f\n', p_perm)
% xx = mean(NeShuff(:, times >= -25 & times <= 25), 2);
% yy = mean(NeA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode ipsi rA vs shuff p: %.4f\n', p_perm)

% axes(ha2(1))
figure(2)
subplot(1,2,1)
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255); hold on;
f.LineStyle = '--';
f.LineWidth = 1;
[bl2,bf2] = boundedline(-win:binSz:win, mean(Np), std(Np)/sqrt(size(Np,1)), 'r'); hold on; 
bl2.Color = regionColors(6,:);
bl2.LineWidth = 1.5;
bf2.FaceColor = regionColors(6,:);
bf2.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NpA), std(NpA)/sqrt(size(NpA,1))); hold on; 
bl3.Color = 0.8*bl2.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

[bl4,bf4] = boundedline(-win:binSz:win, mean(NpShuff), std(NpShuff)/sqrt(size(NpShuff,1))); hold on; 
bl4.Color = 0.5*bl2.Color;
bl4.LineWidth = 1.5;
bf4.FaceColor = bl3.Color;
bf4.FaceAlpha = 0.3;

xlim([-winPlot winPlot])
vl = vline(0); vl.LineWidth =1.5; 
box off

xx = [];
xx(:,1) = mean(Np(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NpA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NpShuff(:, times >= -25 & times <= 25), 2);
p  = anova1(xx);
fprintf('probe ipsi p: %.4f\n', p)
ax = gca;
ax.LineWidth = 1;
% xx = mean(Np(:, times >= -25 & times <= 25), 2);
% yy = mean(NpA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe ipsi coR vs rA p: %.4f\n', p_perm)
% xx = mean(Np(:, times >= -25 & times <= 25), 2);
% yy = mean(NpShuff(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe ipsi coR vs shuff p: %.4f\n', p_perm)
% xx = mean(NpShuff(:, times >= -25 & times <= 25), 2);
% yy = mean(NpA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe ipsi rA vs shuff p: %.4f\n', p_perm)

iRapl = find(contains(regions, regionsPlot([1:5])));
iRbpl = find(contains(regions, regionsPlot([1:5])));

Np = [];
Ne = [];
NpA = [];
NeA = [];
NpShuff = [];
NeShuff = [];
for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl,iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end
            if sum(countPRTHp(iiA,iiB)) == 1; continue; end
            if strcmp(regions{iiA}(1), regions{iiB}(1)); continue; end
            
            times = coR_replay_PRTHe(:,iiA,iiB);
            temp = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRe(iiA,iiB), 'gaussian', smWin);
            
            times = rA_replay_PRTHe(:,iiA,iiB);
            tempA = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countRAe(iiA,iiB), 'gaussian', smWin);
            
            times = coR_replay_PRTHe_shuff(:,iiA,iiB);
            tempShuff = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRe_shuff(iiA,iiB), 'gaussian', smWin);
                        
            z = zscore([temp, tempA, tempShuff]);
            Ne = [Ne; z(1:length(z)/3)];
            NeA = [NeA; z((length(z)/3)+1:(2*length(z)/3))];
            NeShuff = [NeShuff;  z((2*length(z)/3)+1:end)];

            times = coR_replay_PRTHp(:,iiA,iiB);
            temp = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRp(iiA,iiB), 'gaussian', smWin);

            times = rA_replay_PRTHp(:,iiA,iiB);
            tempA = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countRAp(iiA,iiB), 'gaussian', smWin);
            
            times = coR_replay_PRTHp_shuff(:,iiA,iiB);
            tempShuff = smoothdata(histcounts(times(~isnan(times)),-win-(binSz/2):binSz:win+(binSz/2))/countCoRp_shuff(iiA,iiB), 'gaussian', smWin);
                  
            z = zscore([temp, tempA, tempShuff]);
            Np = [Np; z(1:length(z)/3)];
            NpA = [NpA; z((length(z)/3)+1:(2*length(z)/3))];
            NpShuff = [NpShuff;  z((2*length(z)/3)+1:end)];
            


        end
end
    
fprintf('\n')




figure(1)
subplot(1,2,2)
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255); hold on;
f.LineStyle = '--';
f.LineWidth = 1;
% axes(ha1(2))
[bl1,bf1] = boundedline(-win:binSz:win, mean(Ne), std(Ne)/sqrt(size(Ne,1)), 'b'); hold on; 
bl1.Color = contraColor(7,:);
bl1.LineWidth = 1.5;
bf1.FaceColor = contraColor(7,:);
bf1.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NeA), std(NeA)/sqrt(size(NeA,1))); hold on; 
bl3.Color = 0.8*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NeShuff), std(NeShuff)/sqrt(size(NeShuff,1))); hold on; 
bl3.Color = 0.5*bl1.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;
xlim([-winPlot winPlot])
vl = vline(0); vl.LineWidth =1.5; 
fig = gcf;
fig.Color = 'w';
box off

savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplayPRTH_%s_E.pdf', tag)))
times = -win:binSz:win;
xx = [];
xx(:,1) = mean(Ne(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NeA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NeShuff(:, times >= -25 & times <= 25), 2);
p  = anova1(xx);
fprintf('encode contra p: %.4f\n', p)
ax = gca;
ax.LineWidth = 1;
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode contra coR vs rA p: %.4f\n', p_perm)
% xx = mean(Ne(:, times >= -25 & times <= 25), 2);
% yy = mean(NeShuff(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode contra coR vs shuff p: %.4f\n', p_perm)
% xx = mean(NeShuff(:, times >= -25 & times <= 25), 2);
% yy = mean(NeA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('encode contra rA vs shuff p: %.4f\n', p_perm)
figure(2)
subplot(1,2,2)
f = fill([-25 25 25 -25], [-1 -1 2 2], [233 247 244]/255); hold on;
f.LineStyle = '--';
f.LineWidth = 1;
% axes(ha2(2))
[bl2,bf2] = boundedline(-win:binSz:win, mean(Np), std(Np)/sqrt(size(Np,1)), 'r'); hold on; 
bl2.Color = contraColor(7,:);
bl2.LineWidth = 1.5;
bf2.FaceColor = contraColor(7,:);
bf2.FaceAlpha = 0.3;

[bl3,bf3] = boundedline(-win:binSz:win, mean(NpA), std(NpA)/sqrt(size(NpA,1))); hold on; 
bl3.Color = 0.8*bl2.Color;
bl3.LineWidth = 1.5;
bf3.FaceColor = bl3.Color;
bf3.FaceAlpha = 0.3;

[bl4,bf4] = boundedline(-win:binSz:win, mean(NpShuff), std(NpShuff)/sqrt(size(NpShuff,1))); hold on; 
bl4.Color = 0.5*bl2.Color;
bl4.LineWidth = 1.5;
bf4.FaceColor = bl3.Color;
bf4.FaceAlpha = 0.3;
xlim([-winPlot winPlot])
vl = vline(0); vl.LineWidth =1.5;
fig = gcf;
fig.Color = 'w';
box off

xx = [];
xx(:,1) = mean(Np(:, times >= -25 & times <= 25), 2);
xx(:,2) = mean(NpA(:, times >= -25 & times <= 25), 2);
xx(:,3) = mean(NpShuff(:, times >= -25 & times <= 25), 2);
p  = anova1(xx);
fprintf('probe contra p: %.4f\n', p)
ax = gca;
ax.LineWidth = 1;
% xx = mean(Np(:, times >= -25 & times <= 25), 2);
% yy = mean(NpA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe contra coR vs rA p: %.4f\n', p_perm)
% xx = mean(Np(:, times >= -25 & times <= 25), 2);
% yy = mean(NpShuff(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe contra coR vs shuff p: %.4f\n', p_perm)
% xx = mean(NpShuff(:, times >= -25 & times <= 25), 2);
% yy = mean(NpA(:, times >= -25 & times <= 25), 2);
% [~,~,p_perm] = statcond({xx',yy'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% fprintf('probe contra rA vs shuff p: %.4f\n', p_perm)

adj_p = nan(size(pValue));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(pValue))]=fdr_bh(pValue(~isnan(pValue)),0.05,'pdep','yes');


figure(2)
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipReplayPRTH_%s_P.pdf', tag)))




















