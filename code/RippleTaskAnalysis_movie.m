

close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
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
win = [-1 3] * 1e3;
winResp = [-2 1] * 1e3;
recog_trials_all_region = [];

recog_trials_units_all_region = [];


computeCoRipple = false;
rippleFR = nan(length(subj_list_full), 100);
CoRippleRatesRegionAll = nan(length(regions), length(regions), length(subj_list_full));

%prealocate variables
recog_trials_all = nan(6.5e3, length(win(1):win(2)));
recog_trials_units_all = nan(6.5e3, length(win(1):win(2)));

recog_trials_resp_all = nan(6.5e3, length(winResp(1):winResp(2)));
% recog_trials_all_region = [];

whole_trial_all = nan(6.5e3, length(win(1):win(2)));
whole_trials_region_all = nan(6.5e3, length(win(1):win(2)), length(regions));
whole_trials_units_region_all  = nan(6.5e3, length(win(1):win(2)), length(regions));
recog_trials_cross_region_all = cell(length(regions));
whole_trials_unit_cross_region_all = cell(length(regions), length(regions), 3);
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        recog_trials_cross_region_all{iRa, iRb} = nan(6.5e3, length(win(1):win(2)));       
        for ii = 1:size(whole_trials_unit_cross_region_all,3)
            whole_trials_unit_cross_region_all{iRa, iRb,ii} = nan(6.5e3, length(win(1):win(2)));
        end
    end
end

subjID = nan(6.5e3, 1);
correct_trials_all = nan(6.5e3, 1);
novel_trials_all = nan(6.5e3, 1);
resp_confidence = nan(6.5e3, 1);
respLatency = nan(6.5e3, 1);
densityAll = [];
cTrial = 1;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};

    fprintf('processing %s ... \n', subject)

    f = contains(LFPfilesMacr, subject);
    macrObj = load(fullfile(dataDirectory, LFPfilesMacr{f}));
    LFPtimeMacr = macrObj.times;

    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    LFPtimeMicr = micrObj.times;

    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, taskfiles{f}));


    modifier = 'IISremv';
    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'macro_bp', tag);
    load(fullfile(matExportFolder, filename))
    
    modifier = 'IISremv_singleWire_v3_z25';
    modifier = '1kHz_template_z25';
%     modifier = 'IISremv_v3';
    tag = [recordingState,'_',location,'_',modifier];
%     filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
        
    rippleStats.recordingType = [repmat({'macro'}, [1 length(rippleStats.locs)]) repmat({'micro'}, [1 length(microObj.rippleStats.locs)])];
    rippleStats.locs = [cellfun(@(X) LFPtimeMacr(X), rippleStats.locs, 'UniformOutput', false), ...
                        cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false)];
    rippleStats.window(cellfun(@(X) isempty(X), rippleStats.window)) = {[1 1]};
    microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = [cellfun(@(X) [LFPtimeMacr(X(:,1)); LFPtimeMacr(X(:,2))]', rippleStats.window, 'UniformOutput', false), ...
                          cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', microObj.rippleStats.window, 'UniformOutput', false)];

    rippleStats.chanLabels = [rippleStats.chanLabels; microObj.rippleStats.chanLabels];
    rippleStats.density = [rippleStats.density, microObj.rippleStats.density];
    rippleStats.macroTimes = LFPtimeMacr;
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

    rippAll = sum(rippMask, 'omitnan');
    
    units = LoadSpikeTimes(subject,'RutishauserLab', 'bmovie');

    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3)); 
    uLocations = units(U,5);

    unitsAll = find(U);
    nNeuron = sum(U);
    binWidth = 0.001; %seconds
    
    binEdges = 0:binWidth:(rippleStats.recordingLength/1000);
    spikeMatrix = nan(nNeuron, length(binEdges)-1);
    spikeMatrix_sm = nan(nNeuron, length(binEdges)-1);
    for iU = 1:nNeuron 
        X = units{unitsAll(iU),2};
    %     ISIx = [0 diff(X')];
    %     xx = randperm(length(ISIx), length(ISIx));
    %     Xshuff =[X(1)];
    %     for ii = 2:length(ISIx); Xshuff(ii) = Xshuff(end) + ISIx(xx(ii)); end
        [N,EDGES] = histcounts(X,binEdges);
        FR = sum(N) / rippleStats.recordingLength * 1000;
    %     if FR > 0
        spikeMatrix(iU,:) = N;
        spikeMatrix_sm(iU,:) = smoothdata(N, 'gaussian',200);    %     end
    end
    unitAll = smoothdata(sum(spikeMatrix, 'omitnan'), 'gaussian',200);
    unitAll = zscore(unitAll); % - mean(unitAll);

    recog_trials = nan(length(trials.start_time)-1, diff(win)+1);
    recog_trials_resp = nan(length(trials.start_time)-1, diff(winResp)+1);
    recog_trials_region = nan(length(trials.start_time)-1, diff(win)+1, length(regions));
    recog_trials_cross_region = nan(length(trials.start_time)-1, diff(win)+1, length(regions), length(regions));
    recog_trials_units = nan(length(trials.start_time)-1, diff(win)+1);

    dat_region = cell(1,length(regions));
    u_dat_region_sm = cell(1,length(regions));
    u_dat_region = cell(1,length(regions));
    for iT = 1:length(trials.start_time)-1
        respTime = round(trials.response_time(iT+1)*1e3);
        trialTime = round(trials.start_time(iT+1)*1e3);
        respLatency(cTrial)= respTime-trialTime;


        %overal co-rippling
        trialData = nan(1, win(2)-win(1)+1);
        dat = rippAll(trialTime+win(1):trialTime+win(2));
        if length(dat) > length(trialData)
            trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;
        end
        recog_trials_all(cTrial, :) = trialData;

        dat = rippAll(respTime+winResp(1):respTime+winResp(2));
        recog_trials_resp_all(cTrial, :) = dat;

%         if sum(dat) == 0; continue; end


        
        %within region co-rippling
        for iR = 1:length(regions)
            chRegion = contains(rippleStats.chanLabels, regions(iR));
            uRegion = contains(uLocations, regions(iR));

            if isempty(dat_region{iR}) && sum(chRegion) > 0
                 dat_region{iR} = sum(rippMask(chRegion,:),1, 'omitnan');
            elseif sum(chRegion) == 0
                 dat_region{iR} = nan(1, length(rippMask));
            end
            
             if isempty(u_dat_region_sm{iR}) && sum(uRegion) > 0
                 u_dat_region_sm{iR} = zscore(sum(spikeMatrix_sm(uRegion,:),1, 'omitnan'));
                 u_dat_region{iR} = zscore(sum(spikeMatrix(uRegion,:),1, 'omitnan'));
            elseif sum(uRegion) == 0
                 u_dat_region_sm{iR} = nan(1, length(rippMask));
                 u_dat_region{iR} = nan(1, length(rippMask));
             end
            
            
            
            
            
            
            trialData  = dat_region{iR}(trialTime+win(1):trialTime+win(2));
            utrialData = u_dat_region{iR}(trialTime+win(1):trialTime+win(2));
            utrialDatasm = u_dat_region_sm{iR}(trialTime+win(1):trialTime+win(2));            
            
            if length(trialData) <= 16e3
                whole_trials_region_all(cTrial, 1:length(trialData), iR) = trialData;
                whole_trials_units_region_all(cTrial, 1:length(trialData), iR) = utrialDatasm;
            else
               whole_trials_region_all(cTrial, :, iR) = trialData(1:size(whole_trials_region_all,2));
               whole_trials_units_region_all(cTrial, :, iR) = utrialDatasm(1:size(whole_trials_region_all,2));

            end       
            
            
                        
        end
        
        
        
         %across region co-rippling
        for iRa = 1:length(regions)
            for iRb = iRa:length(regions)
                datArip =  whole_trials_region_all(cTrial,:,iRa);
                datBrip =  whole_trials_region_all(cTrial,:,iRb);
                
                
                
                if iRa == iRb
                    datAB = datArip;
                else

                    datAB = datArip + datBrip;
                    datAB(datArip <= 0 | datBrip <= 0) = 0;
                end
                
                recog_trials_cross_region_all{iRa, iRb}(cTrial, :) = datAB;
            end
        end

        trialData = nan(1, win(2)-win(1)+1);
        dat = unitAll(trialTime+win(1):trialTime+win(2));
        if length(dat) > 16e3
                trialData = dat(1:length(trialData));
        else
            trialData(1:length(dat)) = dat;

        end
        recog_trials_units_all(cTrial,:)= trialData;
        
        
        correct_trials_all(cTrial) = trials.response_correct(iT+1);
        novel_trials_all(cTrial) = contains(trials.stimulus_file(iT+1), 'new');

        resp_confidence(cTrial) =  trials.response_confidence(iT+1);
        subjID(cTrial) = subj;

        cTrial = cTrial + 1;
    end

    


    

%     save(fullfile(exportDir, 'ripDetections', sprintf('%s_rippleStats.mat', subject)), 'rippleStats' ,'-v7.3')

end

fprintf('done processing .. \n')
confPerSubj = arrayfun(@(X) sum(resp_confidence(subjID == X) == 3) / sum(subjID == X), unique(subjID(~isnan(subjID))));
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
% figure('Position', [904 414 700 522]);
% 
% % plot_condition = true(1, length(correct_trials_all));
% plot_condition = resp_confidence ==3;
% plot_data = smoothdata(recog_trials_resp_all(plot_condition,:), 2, 'gaussian',1);
% 
% [bl1, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
%                                 std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% 
% plot_data = smoothdata(recog_trials_resp_all(~plot_condition,:), 2, 'gaussian',1);
% 
% [bl2, bf] = boundedline(winResp(1):winResp(2), mean(plot_data, 'omitnan'), ...
%                               std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% 
% figure('Position', [904 414 700 522]);
% subplot(2,3,1)
% plot_condition = correct_trials_all == 1;
% plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',200);
% [bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                                 std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',200);
% 
% [bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                               std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% legend([bl1, bl2], 'correct', 'incorrect', 'location', 'best')
% xlim([-1000 3000])
% ylim([0.7 1.9])
% ylabel('co-ripples')
% xlabel('time from trial start [ms]')
% 
% 
% subplot(2,3,4)  
% % figure;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% ylim([-0.3 0.3])
% xlim([-1000 3000])
% ylabel('unit firing [z-score]')
% xlabel('time from trial start [ms]')
% fig = gcf;
% fig.Color = 'w';
% 
% subplot(2,3,2)
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));
% plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',200);
% [bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                                 std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',200);
% 
% [bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                               std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% legend([bl1, bl2], 'novel', 'familiar', 'location', 'best')
% xlim([-1000 3000])
% ylim([0.7 1.9])
% ylabel('co-ripples')
% xlabel('time from trial start [ms]')
% 
% subplot(2,3,5)  
% % figure;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% ylim([-0.3 0.3])
% xlim([-1000 3000])
% ylabel('unit firing [z-score]')
% xlabel('time from trial start [ms]')
% fig = gcf;
% fig.Color = 'w';
% 
% subplot(2,3,3)
% plot_condition = resp_confidence == 3;
% plot_data = smoothdata(recog_trials_all(plot_condition,:), 2, 'gaussian',300);
% [bl1, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                                 std(plot_data, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% plot_data = smoothdata(recog_trials_all(~plot_condition,:), 2, 'gaussian',300);
% 
% [bl2, bf] = boundedline(win(1):win(2), mean(plot_data, 'omitnan'), ...
%                               std(plot_data, 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% legend([bl1, bl2], 'confident', 'unsure', 'location', 'best')
% xlim([-1000 3000])
% ylim([0.7 1.9])
% % ylim([0 1.5])
% ylabel('co-ripples')
% xlabel('time from trial start [ms]')
% 
% subplot(2,3,6)  
% % figure;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(plot_condition,:), 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
% bf.FaceAlpha = 0.7;
% [bl, bf] = boundedline(win(1):win(2), mean(recog_trials_units_all(~plot_condition,:), 'omitnan'), ...
%     std(recog_trials_units_all(~plot_condition,:), 'omitnan')/sqrt(sum(~plot_condition)), 'b'); hold on;
% bf.FaceAlpha = 0.7;
% vline(0)
% ylim([-0.3 0.3])
% xlim([-1000 3000])
% ylabel('unit firing [z-score]')
% xlabel('time from trial start [ms]')
% fig = gcf;
% fig.Color = 'w';
% savepdf(gcf, fullfile(exportDirFigs, 'coRipTask.pdf'))
% 
regionColors =  brewermap(12, 'Dark2');

figure('Position',[669 805 168 729]);
plot_condition = resp_confidence ==3;
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));
% plot_condition = logical(correct_trials_all(~isnan(correct_trials_all)));

blall = [];
binSz = 1000;
bins = -1000:binSz:5000;
times = win(1):win(2);
for iR = 1:length(regions)
    subplot(5,1,iR)

    plot_data_conf = smoothdata(whole_trials_region_all(plot_condition,:, iR), 2, 'gaussian',100);    
    [bl1, bf] = boundedline(win(1):win(2), mean(plot_data_conf, 'omitnan'), ...
        std(plot_data_conf, 'omitnan')/sqrt(sum(plot_condition)), 'r'); hold on;
     bl1.Color = regionColors(iR,:);
%     
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;


    plot_data_uns = smoothdata(whole_trials_region_all(~plot_condition,:, iR), 2, 'gaussian',100);
    [bl2, bf] = boundedline(win(1):win(2), mean(plot_data_uns, 'omitnan'), ...
            std(plot_data_uns, 'omitnan')/sqrt(sum(~plot_condition)), 'k'); hold on;
    bf.FaceAlpha = 0.7;

    
    hl = hline(mean(whole_trials_region_all(~plot_condition,1:1e3, iR), 'all','omitnan'));
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

    title(regions{iR})

%     legend([bl1, bl2], 'confident', 'unsure', 'location', 'best')
%     legend([bl1, bl2], 'novel', 'familiar', 'location', 'best')
%     legend([bl1, bl2], 'correct', 'incorrect', 'location', 'best')


%     blall(iR) =bl;
    xlim([-1000 3000])
%     ylim([0.3 1])
    ylabel('co-ripples')
    xlabel('time from trial start [ms]')
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'coRipTaskWithinRegions.pdf'))

figure('Position',[669 805 168 729]);

for iR = 1:length(regions)
    subplot(5,1,iR)



    plot_data_1 = whole_trials_units_region_all(~plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));

    [bl2, bf] = boundedline(win(1):win(2), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    
    plot_data_3    = whole_trials_units_region_all(plot_condition,:,iR);
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    

    [bl1, bf] = boundedline(win(1):win(2), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
    
    vl = vline(0); 
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
    
    probeResp_1 = trapz(plot_data_1(:,abs(win(1))+1:abs(win(1))+1.5e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,abs(win(1))+1:abs(win(1))+1.5e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s firing p = %.2f\n', regions{iR}, p_perm)

    title(regions{iR})

%     legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


%     blall(iR) =bl;
    xlim([-1000 3000])
%     ylim([0.3 1])
    ylabel('unit firing')
    xlabel('time from trial start [ms]')
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('FiringTaskWithinRegions_%s.pdf', 'Movie')))
%%
figure('Position',[1277 763 688 810]);

confSubj = find(confPerSubj > 0);
% c = 1;
tabCoRip = cell(5);
lmeCoRip = cell(5);
lmePval = nan(5);
for iRa = 1:length(regions)
    for iRb = iRa:length(regions)
        c = (iRb - 1) * 5 + iRa;
        subplot(5,5,c)


        plot_data_1 = smoothdata(recog_trials_cross_region_all{iRa,iRb}(~plot_condition & ismember(subjID, confSubj),:), 2, 'gaussian',250);
%         plot_data(isnan(whole_trials_cross_region_load1_all{iRa,iRb})) = nan;
        plot_data_mu_1 = mean(plot_data_1, 'omitnan');
        plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
%         plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
%         plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
        [bl2, bf] = boundedline(win(1):win(2), plot_data_mu_1, plot_data_sem, 'k'); hold on;
        bf.FaceAlpha = 0.7;

        plot_data_3    = smoothdata(recog_trials_cross_region_all{iRa,iRb}(plot_condition & ismember(subjID, confSubj),:), 2, 'gaussian',250);

        Arate    = smoothdata(whole_trials_region_all(plot_condition,:,iRa), 2, 'gaussian',500);
%         Arate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
        Arate = mean(Arate, 'all' , 'omitnan');
        Brate = smoothdata(whole_trials_region_all(plot_condition,:,iRb), 2, 'gaussian',500);
%         Brate(:,isnan((whole_trial_all(condTrial1,:))))  = nan;
        Brate = mean(Brate, 'all' ,'omitnan');
        plot_data_3_AB    = Arate.*Brate;
%         plot_data(isnan(whole_trials_cross_region_load3_all{iRa,iRb})) = nan;
        plot_data_mu_3 = mean(plot_data_3, 'omitnan');
        plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
%         plot_data_mu_3(isnan(whole_trial_all(condTrial1,:)))  = nan;  
%         plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

        [bl1, bf] = boundedline(win(1):win(2), plot_data_mu_3, plot_data_sem, 'r'); hold on;
        bf.FaceAlpha = 0.7;
        bl1.Color = mean([regionColors(iRa,:);regionColors(iRb,:)]) ;
        bf.FaceColor = mean([regionColors(iRa,:);regionColors(iRb,:)]);
        bl1.MarkerFaceColor = bl1.Color;
        bf.FaceAlpha = 0.3;
%         pl = hline(plot_data_3_AB, '--'); hold on;
%         pl.Color = bl1.Color;


        vl = vline(0); 
        vl.LineWidth = 1.5;

        
        hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));   
        hl.Color = 'k';
        hl.LineStyle = '-';
        hl.LineWidth = 1.0;
        
%         ylim([0 0.6])
        xlim([-1000 3500])
        
        
     
        

         coRresp = mean(recog_trials_cross_region_all{iRa,iRb}(:,abs(win(1))+1:abs(win(1))+1.5e3), 2, 'omitnan');
         coRbaseline = mean(recog_trials_cross_region_all{iRa,iRb}(:,1:1e3), 2, 'omitnan');
%          coRresp(cTrial:end) = [];
%          coRbaseline(cTrial:end) = [];
        
         tabCoRip{iRa,iRb} = table(subjID(~isnan(coRresp) & ismember(subjID, confSubj)), plot_condition(~isnan(coRresp) & ismember(subjID, confSubj)), coRresp(~isnan(coRresp) & ismember(subjID, confSubj)),coRbaseline(~isnan(coRresp) & ismember(subjID, confSubj)), ...
                 'VariableNames', {'subject',               'condition',                     'rates_resp',            'rates_baseline'});
         lmeCoRip{iRa,iRb} = fitlme(tabCoRip{iRa,iRb}, 'rates_resp ~ condition + (1|subject)');   
%          lmeCoRip{iRa,iRb} = fitlme(tabCoRip{iRa,iRb}, 'rates_baseline ~ rates_maintenance + (1|subject)');   
         lmePval(iRa,iRb) = lmeCoRip{iRa,iRb}.Coefficients(2,6);


%         title([regions{iRa}, ' <--> ', regions{iRb}])
%         if iRa == iRb
%             legend([bl1, bl2], 'confident', 'unsure', 'location', 'best')
%         end
        probeResp_1 = trapz(plot_data_1(:,abs(win(1))+1:abs(win(1))+1.5e3), 2); probeResp_1(isnan(probeResp_1)) = [];
        probeResp_3 = trapz(plot_data_3(:,abs(win(1))+1:abs(win(1))+1.5e3), 2); probeResp_3(isnan(probeResp_3)) = [];
        [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
% %         [pRank,H] = ranksum(probeResp_3, probeResp_1);
% 
        fprintf('%s <--> %s ripples response p = %.2f\n', regions{iRa}, regions{iRb}, p_perm)
% 
%         probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_1(isnan(probeResp_1)) = [];
%         probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
%         probeResp_3(isnan(probeResp_3)) = [];
%         [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
%         [pRank,H] = ranksum(probeResp_3, probeResp_1);

%         fprintf('%s <--> %s ripples maintenance p = %.2f\n', regions{iRa}, regions{iRb},  p)


%         blall(iR) =bl;
    %     xlim(Lim)
    %     ylim([0.3 1])
        ylabel('co-ripples')
        xlabel('time from trial start [ms]')
        
%         c = c+1;
    end
end
fig= gcf;
fig.Color = 'w';
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'coRipTaskAcrossRegions_movie.pdf'))
adj_p = nan(size(lmePval));
[h, crit_p, adj_ci_cvrg, adj_p(~isnan(lmePval))]=fdr_bh(lmePval(~isnan(lmePval)),0.05,'pdep','yes');

figure('Position', [1210 1274 420 299]); 
imagesc(adj_p', [min(adj_p(:))-1e-2 1]); hold on;
colorbar;
hl = hline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
vl = vline((1/2):size(adj_p,1)+(1/2), 'k-'); hold on;
cmap = flipud(slanCM('viridis'));
cmap = [1 1 1; cmap];
colormap(cmap)
[yy xx] = find(adj_p' < 0.05);
plot(xx, yy, 'r*'); hold on;

% subplot(1,2,2)
% for iR = 1:length(regions)
%     [bl] = plot(0,0); hold on;
%     bl.Color = regionColors(iR,:);
% 
%     bf.FaceColor = regionColors(iR,:);
%     bf.FaceAlpha = 0.7;
% 
%     blall(iR) =bl;
% end
% legend(blall, regions, 'location', 'northwest');
% axis off


figure('Position', [899 818 710 106]);
h = boxplot(respLatency, 'Orientation', 'horizontal', 'Symbol', '');
xlabel('Response Latency')
whiskerLimits = [quantile(respLatency, 0.01), quantile(respLatency, 0.95)];
xlim([0 whiskerLimits(2)]);
fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, 'respLatency.pdf'))

% ylim([0.9 1.1])
% savepdf(gcf, 'MemoryResponses_microOnly.pdf')
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
savepdf(gcf, fullfile(exportDirFigs, 'density.pdf'))


























