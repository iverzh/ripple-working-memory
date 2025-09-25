
close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%

subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
dataDirectory =     '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/preprocess';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end

state = 'wake';
rippleContactType = 'micro';

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfiles = dir(fullfile(dataDirectory, sprintf('*%s*', rippleContactType)));
LFPfiles = {LFPfiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../out', '*micro*'));
micrFiles = {micrFiles.name}';
timeShift = 0; %shift unit times to match LFP;


recordingState = 'wake';
location = 'NC';
c = 0;
% figure('Position', [1219 1160 1552 413]); 
% figure;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    units = LoadSpikeTimes(subject,'RutishauserLab');
    c = c + length(units);
%     for iU = 1:size(units,1)
%         subplot(1,2,1)
%         plot(cell2mat(units(iU,4))');
%         subplot(1,2,2)
%         times = units{iU,2}*1e3;
%         acg = times - times';
%         acg(acg == 0) = [];
%         histogram(acg(:), -100.5:100.5)
% 
%         waitforbuttonpress; clf;
%     end
    ccg = nan(size(units,1));
    for a = 1:length(ccg)
        for b = a+1:length(ccg)
            if a == b; continue; end
            if strcmp(units{a,end}(1:4), units{b,end}(1:4)); continue; end
            timesA = units{a,2}*1e3;
            timesB = units{b,2}*1e3;
            ccgAB = timesA - timesB';
            [N, times] = histcounts(ccgAB(:), -100.5:100.5);
            times = movmean(times,2); times(1) = [];
            baseline = mean(N(abs(times) > 50));
            middle = mean(N(abs(times) < 5));
            if ( length(timesA) > 200 ) && ( length(timesB) > 200)
                ccg(a,b) = middle/baseline;
            else
                ccg(a,b) = nan;
            end

        end
    end
    figure;
    imagesc(ccg, [0 10]); colorbar;
    title(subject)
%     waitforbuttonpress; clf;
    
    fldr =     '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/artifactSearch';
    savepdf(gcf, sprintf('%s/%s_ccgPeaks.pdf', fldr, subject))
    save(sprintf('%s/%s_ccgPeaks.mat', fldr, subject), 'ccg', '-v7.3')
%     clf;
end

%%

nRipTotal = 0;
Npyr = [];
Nint = [];
ISI = [];
INphase = [];
PYphase = [];
rippleLFP_py = {};
rippleLFP_in = {};
uc = 1;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    units = LoadSpikeTimes(subject,'RutishauserLab');

    tag = [recordingState,'_',location,'_',rippleContactType];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    load(fullfile(matExportFolder, filename));
    chan_labels = rippleStats.chanLabels;
%     f = contains(micrFiles, subject);
%     micrObj = load(fullfile(dataDirectory, '../out', micrFiles{f}));

    f = contains(LFPfiles, subject);
    micrObj = load(fullfile(dataDirectory, LFPfiles{f}));
    
    f = contains(bpFiles, subject);
    bpObj = load(fullfile(dataDirectory, '../out', bpFiles{f}));
    
    
    
    
    window = 1.500; %seconds
    % window = window * rippMedial.rippleStats.fs;
    window = window * rippleStats.fs;
    
    
    
    % RippleSpikeTrain = nan(length(rippMedial.rippleStats.chanLabels),1e1);
    RippleSpikeTrain = nan(length(rippleStats.chanLabels),1e1);
    for ui = 1:size(units,1)
        switch subject 
            case 'MG67'
                fileName = units{ui,1};
                chanLabel = fileName(18:19);
                if isnan(chanLabel); chanLabel = str2double(fileName(18)); end
    
                chanLabel = find(ismember(rippleStats.chanLabels,chanLabel));
            otherwise
    
    %             chanLabel = units{ch,1} - 96;
                chanLabel = units{ui,5};% - 96;
                chanLabel = find(ismember(chan_labels,chanLabel));
                if isempty(chanLabel); chanLabel = 1; end

                if strcmp(rippleContactType, 'macro')
                    nL = numel(chan_labels);
                    [alph,numr] = deal(cell(nL,1));
                    for iL = 1:nL
                      nMsk = isstrprop(chan_labels{iL},'digit');
                      if ~nMsk(end) % if label does not end in a number
                        alph{iL} = chan_labels{iL};
                        numr{iL} = NaN;
                      else% if label ends in a number
                        % consider everything before the number at the end the name aka alph
                        % (no matter how many digits that number is)
                        idx = find(nMsk,1,'first');
                        alph{iL} = chan_labels{iL}(1:idx-1);
                        numr{iL} = chan_labels{iL}(idx:end);
                      end
                    end
%                     numr = cell2mat(numr)

                    unitLoc = units{ui, 5};
                    nMsk = isstrprop(unitLoc,'digit');
                    idx = find(~nMsk,1,'last');
                    unitLoc = unitLoc(1:idx);

                    chanLabel = find(strcmp(alph, unitLoc) & strcmp(numr, '02-01'));


                end
        end
    
        if strcmp(rippleContactType, 'micro'); dataChan = micrObj.lfp_data(chanLabel, :);
        elseif strcmp(rippleContactType, 'macro'); dataChan = bpObj.data(:,chanLabel)';
        end
%         
    
    %     spikeTimes = units{ch,2} * rippMedial.rippleStats.fs;
        spikeTimes = (units{ui,2} + timeShift) * rippleStats.fs;
        spikeTimes(spikeTimes > length(dataChan)) = [];
        if strcmp(subject,'MG63')
            spikeTimes(spikeTimes > 4229990 & spikeTimes < 15900000) = [];
        end
         
        nSpk  = length(spikeTimes);
    %     nSpk/secs
        chLFP = units{ui,1};
        if nSpk > 0 && ~strcmp(units{ui,3},'mult') && ~isempty(units{ui,3})
            %compute spike PRTH
    %         tempCh = find(ismember(newLabels,str2double(a.rippleStats.chanLabels{chanLabel})));
    %         ripples = rippLateral.rippleStats.locs{chanLabel};
            ripples = rippleStats.locs{chanLabel};
            
            nRip = length(ripples);
            rippleMat = zeros(nRip,2*window+1);
            PRTH = [];
            c = 1;
            phasesUnit = [];
            for rip = 1:length(ripples)
                ripLoc = ripples(rip);
                ripWin = micrObj.times(rippleStats.window{chanLabel}(rip,:));
                if ripLoc + window < length(dataChan)
                    rippleMat(rip,:) = dataChan(ripLoc-window:ripLoc+window);
                    ripLoc = micrObj.times(ripLoc) * rippleStats.fs;
                    localSpikes = spikeTimes(spikeTimes > ripLoc-window & spikeTimes < ripLoc+window);
                    localSpikesRip = spikeTimes(spikeTimes > ripWin(1) & spikeTimes < ripWin(2));
    %                 phasesUnit = [phasesUnit RBphaseAll(chanLabel, round(localSpikesRip))];
                    ISItemp = diff(localSpikes);
                    ISI = [ISI, ISItemp];
                    localSpikesShift = localSpikes - ripLoc;
                    PRTH = [PRTH, localSpikesShift];

%                     if isempty(localSpikesRip); rippleMat(rip,:) = []; end

    
    %                 if ~isempty(localSpikes)
    %                     phases  = ripple_phases(chanLabel,round(localSpikes));
    % 
    %                     if strcmp(units{ch,3},'pyr')
    %                         PYphase = [PYphase phases];
    %                     elseif strcmp(units{ch,3},'int')
    %                         INphase = [INphase phases];
    %                     end
    % 
    % 
    %                 end
                end
                
            end
            
            if strcmp(units{ui,3},'pyr')
                [Npyr(uc,:),edges] = histcounts(PRTH,-window:window);
                PYphase{uc} = phasesUnit;
                rippleLFP_py{uc} = rippleMat;
                uc = uc + 1;
            elseif strcmp(units{ui,3},'int')
                [Nint(ui,:),edges] = histcounts(PRTH,-window:window);
                INphase{ui} = phasesUnit;
                rippleLFP_in{ui} = rippleMat;
    
            end
    %         N(chanLabel,:) = zscore(N(chanLabel,:));
            edges = smoothdata(edges,'movmean',2);
            edges(1) = [];
            nRipTotal = nRipTotal + nRip; 
            
            %find spike trains (looking for conditions that underlie plateau
            % potentials)
    %         for rip = 1:length(ripples)
    %             ripLoc = ripples(rip);
    %             rippleStart = rippleStats.window{chanLabel}(rip,1);
    %             rippleEnd = rippleStats.window{chanLabel}(rip,2);
    %             rippleSpikes = spikeTimes > rippleStart & spikeTimes < rippleEnd;
    % 
    % %             if sum(rippleSpikes)>=3
    % %                 spkTemp = spikeTimes(rippleSpikes) - ripLoc + 101;
    % %                 RippleSpikeTrain(chanLabel,c) = ripLoc;
    % %                 c = c + 1;
    % %                 
    % %                 [b,a] = butter(3,rippleStats.RB/(rippleStats.fs/2));
    % %                 rippleband = filtfilt(b,a,dataChan);
    % %                 figure('Position',[356 276 1201 413]);
    % %                 yyaxis left
    % %                 plot(dataChan(ripLoc-100:ripLoc+100)); hold on;
    % %                 yyaxis right
    % %                 plot(rippleband(ripLoc-100:ripLoc+100)); hold on;
    % %                 ln2 = plot(101-(ripLoc-rippleStart):101+(rippleEnd-ripLoc),rippleband(rippleStart:rippleEnd),'g-'); hold on;
    % %                 ln2.LineWidth =2;
    % %                 
    % %                 vl = vline(spkTemp,'r-');
    % %                 for l = 1:length(vl); vl(l).LineWidth= 2.5; end
    % % %                 vl.Color = 'r';
    % % %                 vl.LineWidth = 2.5;
    % %                 
    % %                 RB = rippleband(ripLoc-100:ripLoc+100);
    % %                 mark = plot(spkTemp, repmat(max(RB)+10,[1 length(spkTemp)]), 'o');
    % %                 pause; close all;
    % %             end
    %         end
        end
    
        ui
    end
end
%%
exportDirec =     '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/rippleUnitFiring';
edges = (-window+0.5):(window-0.5);

win = 100;
% close all
[A,I] = sort(max(Npyr,[],2),'descend');
Nsort = Npyr(I,:);
Nsort(A==0,:) = [];
% Nsort = Nsort / 0.5 / nRipTotal * rippleStats.fs;
% N_Hz = Nsort / nRipTotal * rippMedial.rippleStats.fs;
nRipTotal = cellfun(@(X) length(X), rippleLFP_py);
nRipTotal = sum(nRipTotal);
N_Hz = Nsort / nRipTotal * rippleStats.fs;
figure;
subplot(4,1,[2 4])
spikingWindow = N_Hz(:, 1501-win:1501+win);
mu = mean(spikingWindow(:), 'omitnan'); s = std(spikingWindow(:), 'omitnan');
imagesc(N_Hz,[mu-s mu+5*s] ); hold on;% colorbar;
xlim([1501-win 1501+win])
vline(1001)
% colorbar
% rippleMat = cell2mat(rippleLFP_py);
rippleMat = rippleLFP_py{1};
subplot(4,1,1)
plot(-window:window, -nanmean(rippleMat)); hold on;
box off
vline(0)
xlim([-win win])

[ISIcounts,ISIedges] = histcounts(ISI,0:1:50);

ISIedges = smoothdata(ISIedges,'movmean',2);
ISIedges(1) = [];
figure; plot(ISIedges,smoothdata(ISIcounts,'gaussian',5));
vline(11)
grid minor

% close all
figure('Position',[661 327 429 542]);
plotMat = cellfun(@(X) X', rippleLFP_py, 'UniformOutput', false);
plotMat = cell2mat(plotMat)'*1e6;
% rippleMat = cell2mat(rippleLFP_py);
subplot(2,1,1)
boundedline(-window:window, -mean(plotMat, 'omitnan'), std(plotMat, 'omitnan') / sqrt(size(plotMat,1)), 'k'); hold on;
ylabel('microV')
xlim([-win win])

subplot(2,1,2)
yyaxis left
s = cellfun(@(X) size(X,1), rippleLFP_py);
A = Npyr./s' * 1e3;
% sc = 1/max(A);

ln = boundedline(edges,mean(A, 'omitnan'),std(A, 'omitnan') / sqrt(size(A,1)), 'b'); hold on;
ln.LineWidth = 1;
ln.Color = [67 151 186]/255;

ax = gca;Stealth Neurotech
ax.YAxis(1).Color = ln.Color;
ylabel('Firing Rate [Hz]')

yyaxis right
s = cellfun(@(X) size(X,1), rippleLFP_in);
A = Nint./s' * 1e3;

if ~isempty(A)
    % sc = 1/max(A);
    ln = plot(edges,mean(A, 'omitnan'),'r-'); hold on;
    ln.LineWidth = 1;
    ln.Color = [233 48 34]/255;
    ylabel('Firing Rate [Hz]')
    % ylim([0.3 1.1])
    ax = gca;
    ax.YAxis(2).Color = ln.Color;
    fig = gcf;
    fig.Color = [1 1 1];
    box off
    filename = sprintf('%s_%s_ripple-unit.pdf', subject, state);
    savepdf(gcf,fullfile(exportDirec, filename))
end
xlim([-win win])

% subplot(2,1,1)
% plot(-window:window, -nanmean(rippleMat), 'k'); hold on;
% ylabel('microV')
% xlim([-win win])
axis off
filename = sprintf('%s_%s_%s_ripple-unit_noAxis.pdf', subject, state, rippleContactType);
savepdf(gcf,fullfile(exportDirec, filename))

figure('Position',[661 327 429 542]);
subplot(2,1,1)
ln = plot(-window:window, -nanmean(rippleMat), 'k'); hold on;
ln.LineWidth = 2;

ylabel('microV')
xlim([-50 50])

subplot(2,1,2)
yyaxis left
s = cellfun(@(X) size(X,1), rippleLFP_py);
A = Npyr./s' * 1e3;
% sc = 1/max(A);

ln = plot(edges,mean(A, 'omitnan'),'b-'); hold on;
ln.LineWidth = 1.5;
ln.Color = [67 151 186]/255;
ylim([floor(min(mean(A, 'omitnan'))) ceil(max(mean(A, 'omitnan')))])
ax = gca;
ax.YAxis(1).Color = ln.Color;
ylabel('Firing Rate [Hz]')

yyaxis right
s = cellfun(@(X) size(X,1), rippleLFP_in);
A = Nint./s' * 1e3;
if ~isempty(A)
    % sc = 1/max(A);
    ln = plot(edges,mean(A, 'omitnan'),'r-'); hold on;
    ln.LineWidth = 1.5;
    ln.Color = [233 48 34]/255;
    ylabel('Firing Rate [Hz]')
    xlim([-50 50])
    % ylim([0.3 1.1])
    ax = gca;
    ax.YAxis(2).Color = ln.Color;
    fig = gcf;
    fig.Color = [1 1 1];
    box off
    filename = sprintf('%s_%s_ripple-unit_zoom.pdf', subject, state);
    savepdf(gcf,fullfile(exportDirec, filename))
end
subplot(2,1,1)
plot(-window:window, -nanmean(rippleMat), 'k'); hold on;
ylabel('microV')
xlim([-50 50])
axis off
filename = sprintf('%s_%s_ripple-unit_zoom_noAxis.pdf', subject, state);
savepdf(gcf,fullfile(exportDirec, filename))

PYmeanPhase  = [];
for u = 1:length(PYphase)
    if ~isempty(PYphase{u})
        PYmeanPhase = [PYmeanPhase PYphase{u}];
    end
end

INmeanPhase  = [];
for u = 1:length(INphase)
    if ~isempty(INphase{u})
        INmeanPhase = [INmeanPhase INphase{u}];
    end
end

phaseEdges = linspace(-pi,pi,51);



%%

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
cd(exportDir)
cUnit = 0; coRipPad = 0;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    units = LoadSpikeTimes(subject,'RutishauserLab');

    tag = [recordingState,'_',location,'_',rippleContactType];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    load(fullfile(matExportFolder, filename));
    chan_labels = rippleStats.chanLabels;
% %     f = contains(micrFiles, subject);
% %     micrObj = load(fullfile(dataDirectory, '../out', micrFiles{f}));
% 
    f = contains(LFPfiles, subject);
    micrObj = load(fullfile(dataDirectory, LFPfiles{f}));
    times = micrObj.times * 1e3;
%     
%     f = contains(bpFiles, subject);
%     bpObj = load(fullfile(dataDirectory, '../out', bpFiles{f}));
    
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    for chRipp = 1:size(rippMask,1) 
        if ~isempty(rippleStats.window{chRipp}) 
            iS = round(times(rippleStats.window{chRipp}(:,1)));
            iE = round(times(rippleStats.window{chRipp}(:,2)));

        
            for ii = 1:length(iE)            
                if iS(ii) <= 0 || iE(ii) <= 0; continue; end

                rippMask(chRipp,iS(ii):iE(ii)) = 1;
            end
        end
        
    end

    [b,a] = butter(3,[70 100]/(rippleStats.fs/2));
    data = micrObj.lfp_data;
    rippleBand = nan(size(data));
    rippleBandPhase = nan(size(data));
    for ch = 1:size(data,1)
        rippleBand(ch,:) = filtfilt(b,a,data(ch,:));
        rippleBandPhase(ch,:) = angle(hilbert(rippleBand(ch,:))); 
        clc% zero phase
        fprintf('\n\n\n ... calculating ripple phase ...\n')
        fprintf('     %s [ %02i / %02i ]\n', subject, ch, size(data, 1))
    end


    coRipSpike = cell(length(units),length(units));
    coRipDur = cell(length(units),length(units));
    noRipDur = cell(length(units),length(units));
    noRipSpike = cell(length(units),length(units));
    interactionType = cell(length(units),length(units));
    coFireProbNoR = cell(length(units),length(units));
    ccgNoR = cell(length(units),length(units));
    ccgCoR = cell(length(units),length(units));
    coSpikeTimes = cell(length(units),length(units));
    coSpikePhase = cell(length(units),length(units));
    
    for uAi = 1:length(units)
        for uBi = 1:length(units)
    %         uA = unitsA(uAi); uB = unitsB(uBi);
            uA = uAi; uB = uBi;
            chA = units{uA,1}; chB = units{uB,1};
            typeA = units{uA,3}; typeB = units{uB,3};
            
            if chA ~= chB && any(strcmp(typeA,{'pyr', 'int'})) && any(strcmp(typeB,{'pyr', 'int'}))
                uTa = units{uA,2} * 1e3; uTa = uTa(uTa >= 0.5);
                uTb = units{uB,2} * 1e3; uTb = uTb(uTb >= 0.5);
    
                coRipMask = rippMask(chA, :) & rippMask(chB, :);
                noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :);
    
                
    
                uTcoRa = uTa(coRipMask(round(uTa)));
                uTcoRb = uTb(coRipMask(round(uTb)));
    
                uTnoRa = uTa(noRipMask(round(uTa)));
                uTnoRb = uTb(noRipMask(round(uTb)));
    
                controlMask = zeros(size(coRipMask));
    
                if ~isempty(uTcoRa) && ~isempty(uTcoRb)
                    
                    
                    bnds = mask2bounds(coRipMask);
                    dur = bnds(:,2) - bnds(:,1);
                    bnds(dur < 25, :) = [];
                    dur(dur < 25) = [];

                    if size(bnds,1) < 2; continue; end
                    
                    coSpikeCount = nan(2,length(bnds));
                    coSpikeTimesAB = cell(1,length(bnds));
                    coSpikePhaseAB = cell(1,length(bnds));
                    noSpikeCount = nan(2,length(bnds));
                    bndsBaselineAll = nan(size(bnds)); 
                    for ibnd  = 1:length(bnds)
                        aSpikesCo = sum(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                        bSpikesCo = sum(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));
    
                        shift = 2e3;
                        bndsBaseline = bnds(ibnd,:) - shift;
                        if any(bndsBaseline < 0); continue; end
    
                        loopCount = 0;
                        while ~all(noRipMask(bndsBaseline(1):bndsBaseline(2))) && loopCount < 1000  
                            bndsBaseline = bndsBaseline - 100;
                            loopCount = loopCount + 1;
                            if any(bndsBaseline < 1)
                                loopCount = 1000;
                                break
                            end
                        end
    
                        if loopCount > 1000; continue; end
    
                        aSpikesNo = sum(uTnoRa >= bndsBaseline(1)   &  uTnoRa <= bndsBaseline(2));
                        bSpikesNo = sum(uTnoRb >= bndsBaseline(1)   &  uTnoRb <= bndsBaseline(2));
    
                        coSpikeCount(:, ibnd) = [aSpikesCo, bSpikesCo];
                        coSpikeTimesAB{1, ibnd} = uTcoRa(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                        coSpikeTimesAB{2, ibnd} = uTcoRb(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));
                        coSpikePhaseAB{ibnd}(1,:) = rippleBandPhase(chA,bnds(ibnd,1) :bnds(ibnd,2) );
                        coSpikePhaseAB{ibnd}(2,:) = rippleBandPhase(chB,bnds(ibnd,1) :bnds(ibnd,2) );

                        noSpikeCount(:, ibnd) = [aSpikesNo bSpikesNo];
                        bndsBaselineAll(ibnd, :) = bndsBaseline;
    
                    end
    
                   
                    controlMask = bounds2mask(bndsBaselineAll, length(coRipMask));
                    uTcontrola = uTa(controlMask(round(uTa)));
                    uTcontrolb = uTb(controlMask(round(uTb)));
        
                    s = uTcontrola - uTcontrolb';
                    N = histcounts(s, -100.5:100.5);
    
               
                    coRipSpike{uA,uB} = coSpikeCount;
                    coRipDur{uA,uB} = dur;
                    noRipDur{uA,uB} = sum(noRipMask);
                    noRipSpike{uA,uB} = noSpikeCount;
                    coFireProbNoR{uA,uB} = N;
                    
                    

%                     coRipMask = bounds2mask(bnds, length(coRipMask), coRipPad);
%                     uTcoRa = uTa(coRipMask(round(uTa)));
%                     uTcoRb = uTb(coRipMask(round(uTb)));
                    s = uTcoRa - uTcoRb';
                    N = histcounts(s, -100.5:100.5);
                    ccgCoR{uA,uB} = N;

                    coSpikeTimes{uA, uB} = coSpikeTimesAB;
                    coSpikePhase{uA, uB} = coSpikePhaseAB;

                    s = uTnoRa - uTnoRb';
                    N = histcounts(s, -100.5:100.5);
                    ccgNoR{uA,uB} = N;
    
                    interactionType{uA,uB} = sprintf('%s-%s', typeA, typeB);
                end
            end
        end
        clc
        fprintf('... %s co firing unit %i / %i ... \n', subject, uAi, length(units))
         cUnit=cUnit+1;
    end
    
    nCoSpikeCoAll = [];
    nCoSpikeNoAll = [];
    
    distanceAll = [];
    c = 1;
    
    nCoSpikeNo = [];
    nCoSpikeCo = [];
    %     dur = 0;
    for uAi = 1:length(units)
        for uBi = 1:length(units)
    %         uA = units(uAi); uB = units(uBi);
            uA = uAi; uB = uBi;
            chA = units{uA,1};
            chB = units{uB,1};
            d = pdist([micrObj.chan_coords(chA,:); micrObj.chan_coords(chB,:)]);

            
            coRip = coRipSpike{uA,uB};
            noRip = noRipSpike{uA,uB};
            dur = sum(coRipDur{uA,uB}) / 1e3;

            if ~isempty(coRip) 
                sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                nCoSpikeCo = size(sp,2); %/dur; % / (dur);
                
                sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                nCoSpikeNo = size(sp,2); %/dur; % / (dur);
%                     dur = dur + (sum(coRipDur{uA,uB}) / 1e3);
    
                nCoSpikeCoAll = [nCoSpikeCoAll, nCoSpikeCo];
                nCoSpikeNoAll = [nCoSpikeNoAll, nCoSpikeNo];
                distanceAll = [distanceAll, d];
            end
            
                
    
        end
    end
    uCh = units(:,end);
    filename = sprintf('%s_coRipple_coFire_noCoRipPad.mat', subject);
    save(fullfile(exportDir, filename), 'nCoSpikeCoAll', 'nCoSpikeNoAll', 'distanceAll', ...
                                        'coRipSpike', 'coRipDur', 'noRipSpike', 'noRipDur', 'coFireProbNoR', ...
                                        'ccgNoR', 'ccgCoR','uCh', ...
                                        'coSpikeTimes', 'coSpikePhase', '-v7.3')
    

    
end

%%
pdfExport = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
tractLenStruct = load(fullfile(exportDir, '..', 'averageConnectivity_tractLengths.mat'));
ctxParc = {'ACC', 'SMA', 'OFC' ,'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};
nCoSpikeCoAllSubj = []; nCoSpikeNoAllSubj = [];
distanceAllSubj =[]; 
ccgAllcoRpeak = [];
ccgAllcoRtrough = [];
ccgAllnoR = [];
compStringAllSubj = []; phLag = [];       
phLagsAll = nan(1, 1e6); cLags = 1;
coRipPad = 0;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    nCo = []; nNo = []; durationAll = [];
    filename = sprintf('%s_coRipple_coFire_noCoRipPad.mat', subject);
    load(fullfile(exportDir, filename))
    tractDistance = []; compString = [];
     for uA = 1:size(coRipDur, 1)
        for uB = 1:size(coRipDur, 2)
            coRip = coRipSpike{uA,uB};
            noRip = noRipSpike{uA,uB};
            dur = (sum(coRipDur{uA,uB}) + 2*coRipPad*length(coRipDur{uA,uB})) / 1e3;
            durNo = noRipDur{uA,uB} / 1e3;

            if (~isempty(coRip) || ~isempty(noRip)) && ~strcmp(uCh{uA},uCh{uB}) %&& strcmp(interactionType{uA, uB}, interact) % || strcmp(interactionType{uA, uB}, 'int-pyr'))
                sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                nCoSpikeCo = size(sp,2); 
                
                sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                nCoSpikeNo = size(sp,2);
                
                nCo = [nCo, nCoSpikeCo];
                nNo = [nNo, nCoSpikeNo];
                durationAll = [durationAll, dur];

                parcelA = uCh{uA}(2:end-1);
                parcelB = uCh{uB}(2:end-1);
                hemA = uCh{uA}(1);
                hemB = uCh{uB}(1);

                broadmanA = broadman{strcmp(ctxParc, parcelA)};
                HCPa = sprintf('%s_%s', hemA, broadmanA);
                broadmanB = broadman{strcmp(ctxParc, parcelB)};
                HCPb = sprintf('%s_%s', hemB, broadmanB);


                d = tractLenStruct.tractLengths(strcmp(tractLenStruct.parcelIDs, HCPa), ...
                                                strcmp(tractLenStruct.parcelIDs, HCPb));
                if strcmp(hemA, hemB); hem = 'ipsi';
                else; hem = 'contra'; end

                parcelString = [{parcelA}, {parcelB}];
                parcelString = sort(parcelString);

                compString = [compString {sprintf('%s %s %s', parcelString{1},parcelString{2}, hem)}];
                tractDistance = [tractDistance d];

                
%                 phLag =


                if d > 0
                    phMat = coSpikePhase{uA, uB};
                    emptyPh = cellfun(@(X) isempty(X), phMat);
                    phMat(emptyPh) = {zeros(2,50)};
                    phLags = cellfun(@(X) circ_mean([X(1,:) - X(2,:)]'), phMat);
                    phLags(emptyPh) = nan;
                    phLagsAll(cLags:cLags+length(phLags)-1) = phLags;
                    cLags = cLags + length(phLags);

                    zeroLag = phLags >= - pi/2 & phLags <= pi/2;
                    nzeroLag = phLags <= - pi/2 | phLags >= pi/2;
                    
                    spMat = coSpikeTimes{uA, uB};
                    N = ccgCoR{uA, uB};
                    uTa = cell2mat(spMat(1,zeroLag)); uTb = cell2mat(spMat(2,zeroLag));
%                     if ~isempty(uTa) && ~isempty(uTb)
%                         s = uTa - uTb';
%                         N = histcounts(s, -100.5:100.5);
%                     else
%                         N = zeros(1,201);
%                     end

%                     dur = sum(coRipDur{uA,uB}(zeroLag)) / 1e3;
                    N = N/dur;
                    ccgAllcoRpeak = [ccgAllcoRpeak; N];

                    uTa = cell2mat(spMat(1,nzeroLag)); uTb = cell2mat(spMat(2,nzeroLag));
                    if ~isempty(uTa) && ~isempty(uTb)
                        s = uTa - uTb';
                        N = histcounts(s, -100.5:100.5);
                    else
                        N = zeros(1,201);
                    end

                    dur = sum(coRipDur{uA,uB}(nzeroLag)) / 1e3;
                    N = N/dur;
                    ccgAllcoRtrough = [ccgAllcoRtrough; N];
    
                    N = ccgNoR{uA, uB};
                    N = N/durNo;
                    ccgAllnoR = [ccgAllnoR; N];
                end



            end
        end
     end
     
%     zAll = zscore([nCoSpikeCoAll(distanceAll > 0) nCoSpikeNoAll(distanceAll > 0)]);
%     nCoSpikeCoAll = zAll(1:sum(distanceAll > 0));
%     nCoSpikeNoAll = zAll(sum(distanceAll > 0)+1:end);

    zAll = zscore([nCoSpikeCoAll(tractDistance > 0) nCoSpikeNoAll(tractDistance > 0)]);
    nCoSpikeCoAll = zAll(1:sum(tractDistance > 0));
    nCoSpikeNoAll = zAll(sum(tractDistance > 0)+1:end);
    compStringAllSubj = [compStringAllSubj compString(tractDistance > 0)];
    
    nCoSpikeCoAllSubj = [nCoSpikeCoAllSubj nCoSpikeCoAll];
    nCoSpikeNoAllSubj = [nCoSpikeNoAllSubj nCoSpikeNoAll];
    distanceAllSubj = [distanceAllSubj tractDistance(tractDistance > 0)];
    disp(subject)

end
disp('computing CCGs complete.')
 %%

close all


figure('Position',  [846 779 754 542]); 
subplot(2,2,1)
[pl1, pa1] = boundedline(-100:100, mean(ccgAllcoRpeak, 'omitnan'), std(ccgAllcoRpeak, 'omitnan')/sqrt(size(ccgAllcoRpeak,1))); hold on; 
pl1.Color = clr(8,:);
pa1.FaceColor = clr(8,:);
pl1.LineWidth = 2;
pa1.FaceAlpha = 0.3;

[pl2, pa2] = boundedline(-100:100, mean(ccgAllnoR, 'omitnan'), std(ccgAllnoR, 'omitnan')/sqrt(size(ccgAllnoR,1))); hold on; 
% [pl2, pa2] = boundedline(-100:100, mean(ccgAllcoRtrough, 'omitnan'), std(ccgAllcoRtrough, 'omitnan')/sqrt(size(ccgAllcoRtrough,1))); hold on; 
pl2.Color = clr(1,:);
pa2.FaceColor = clr(1,:);
pl2.LineWidth = 2;
pa2.FaceAlpha = 0.3;
grid on;
xlim([-100 100])

ylabel('prop.')
xlabel('time [ms]')
legend([pl1, pl2], 'co-ripple', 'no-ripple', 'location', 'best')


subplot(2,2,2)
D = smoothdata(ccgAllcoRpeak, 2, 'gaussian', 1);

[pl1, pa1] = boundedline(-100:100, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
pl1.Color = clr(8,:);
pa1.FaceColor = clr(8,:);
pl1.LineWidth = 2;
pa1.FaceAlpha = 0.3;

% D = smoothdata(ccgAllcoRtrough, 2, 'gaussian', 5);
D = smoothdata(ccgAllnoR, 2, 'gaussian', 5);
[pl2, pa2] = boundedline(-100:100, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
pl2.Color = clr(1,:);
pa2.FaceColor = clr(1,:);
pl2.LineWidth = 2;
pa2.FaceAlpha = 0.3;
xlim([-50 50])
ylabel('prop.')
xlabel('time [ms]')
grid on;

subplot(2,2,3)
keep = sum(ccgAllcoRtrough, 2) >= 0;
D1 = smoothdata(ccgAllcoRpeak(keep,:), 2, 'gaussian', 10);
D2 = smoothdata(ccgAllcoRtrough(keep,:), 2, 'gaussian', 10);
D = mean(D1, 'omitnan')./mean(D2,'omitnan');
D = smoothdata(ccgAllcoRpeak./ccgAllnoR, 2, 'gaussian', 10);
% D = smoothdata(ccgAllcoR, 2, 'gaussian', 10);
D(~isfinite(D)) = nan;
[pl1, pa1] = boundedline(-100:100, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
% pl1 = plot(-100:100, D); hold on; 
pl1.Color = clr(8,:);
pl1.LineWidth = 2;
% pa1.FaceColor = clr(8,:);
% pa1.FaceAlpha = 0.3;
vline(-99:11:99)
% ylim([2.5 4.5])
xlim([-100 100])
ylabel('coR / noR')
xlabel('time [ms]')

subplot(2,2,4)

[pl1, pa1] = boundedline(-100:100, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
% pl1 = plot(-100:100, D); hold on; 
pl1.Color = clr(8,:);
pl1.LineWidth = 2;
% pa1.FaceColor = clr(8,:);
% pa1.FaceAlpha = 0.3;
% ylim([3 4.5])

vline(-99:11:99)

xlim([-50 50])
ylabel('coR / noR')
xlabel('time [ms]')
fig = gcf;
fig.Color = 'w';

%%
figure('Position', [1000 966 263 355]);

nCoSpikeCoDist = [];
nCoSpikeNoDist = [];
nCoSpikeCoDistSEM = [];
nCoSpikeNoDistSEM = [];
bins = quantile(distanceAllSubj, 0:1/6:1); %:0.01:0.9;


for iB = 1:(length(bins)-1)
    ii = distanceAllSubj >= bins(iB) & distanceAllSubj < bins(iB+1);
    
    nCoSpikeCoDist = [nCoSpikeCoDist, mean(nCoSpikeCoAllSubj(ii))];
    nCoSpikeNoDist = [nCoSpikeNoDist, mean(nCoSpikeNoAllSubj(ii))];

    nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(nCoSpikeCoAllSubj(ii), 1)];
    nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(nCoSpikeNoAllSubj(ii),1)];


end

xVal = movmean(bins, 2);
xVal(1) = [];



clr = brewermap(10,'Paired');

p = polyfit(xVal, nCoSpikeCoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);
b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
b1.MarkerFaceColor = clr(8,:);
b1.Color = clr(8,:);
b1.LineWidth = 1.0;




p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);

b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
b2.MarkerFaceColor = clr(1,:);
b2.Color = clr(1,:);
b2.LineWidth = 1.0;
b2.MarkerSize = 5;


ylabel('CoFire Event Rate [z-score]')
xlabel('fiber tract distance between units [mm]')

fig = gcf;
fig.Color = 'w';
box off

%         xlim([0 4000])

ax = gca;

% savepdf(gcf, fullfile(pdfExport, 'coFireDistance_inclSameBundle.pdf'))
% savepdf(gcf, fullfile(pdfExport, 'coFireDistance.pdf'))
% savepdf(gcf, fullfile(pdfExport, 'unitCCGs_inclSameBundle.pdf'))
% savepdf(gcf, fullfile(pdfExport, 'unitCCGs.pdf'))

compStrings = unique(compStringAllSubj);

figure('Position', [990 518 1822 791]);
for iC = 1:length(compStrings)
    ii = strcmp(compStringAllSubj, compStrings{iC});
    xVal = mean(distanceAllSubj(ii));
    
    if contains(compStrings{iC}, 'ipsi'); subplot(1,2,1);
    else; subplot(1,2,2); end
    xlim([0 280])
    ylim([-0.6 1.2])

    vl = vline(xVal);
    vl.Color = 'k';
    vl.LineWidth = 0.5;

    b2 = errorbar(xVal, mean(nCoSpikeCoAllSubj(ii)), computeSEM(nCoSpikeCoAllSubj(ii),1), 'o--'); hold on;
    b2.MarkerFaceColor = clr(8,:);
    b2.Color = clr(8,:);
    b2.LineWidth = 1.0;
    b2.MarkerSize = 5;

    b2 = errorbar(xVal, mean(nCoSpikeNoAllSubj(ii)), computeSEM(nCoSpikeNoAllSubj(ii),1), 'o--'); hold on;
    b2.MarkerFaceColor = clr(1,:);
    b2.Color = clr(1,:);
    b2.LineWidth = 1.0;
    b2.MarkerSize = 5;

    

    text(xVal+5, mean(nCoSpikeCoAllSubj(ii)), compStrings{iC}); hold on;
    
    box off
    ylabel('co-firing rate [z-score]')
    xlabel('fiber tract distance [mm]')

    
end

fig = gcf;
fig.Color = 'w';
savepdf(gcf, fullfile(pdfExport, 'coFireDistance_parcel.pdf'))

figure('Position', [213 530 1055 348]);

subplot(1,3,1)
iiGroup = cellfun(@(X) contains(X, 'ipsi'), compStringAllSubj);
plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
plotGroupDist = distanceAllSubj(iiGroup);
nCoSpikeCoDist = [];
nCoSpikeNoDist = [];
nCoSpikeCoDistSEM = [];
nCoSpikeNoDistSEM = [];
bins = quantile(plotGroupDist, 0:1/6:1); %:0.01:0.9;


for iB = 1:(length(bins)-1)
    ii = plotGroupDist >= bins(iB) & plotGroupDist < bins(iB+1);
    
    nCoSpikeCoDist = [nCoSpikeCoDist, mean(plotGroupCo(ii))];
    nCoSpikeNoDist = [nCoSpikeNoDist, mean(plotGroupNo(ii))];

    nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(plotGroupCo(ii), 1)];
    nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(plotGroupNo(ii),1)];


end

xVal = movmean(bins, 2);
xVal(1) = [];
clr = brewermap(10,'Paired');

p = polyfit(xVal, nCoSpikeCoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);
b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
b1.MarkerFaceColor = clr(8,:);
b1.Color = clr(8,:);
b1.LineWidth = 1.0;




p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);

b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
b2.MarkerFaceColor = clr(1,:);
b2.Color = clr(1,:);
b2.LineWidth = 1.0;
b2.MarkerSize = 5;

ylim([-0.4 0.8])
xlim([0 300])
box off
ylabel('CoFire Event Rate [z-score]')
xlabel('fiber tract distance between units [mm]')
title('ipsi')

subplot(1,3,2)
iiGroup = cellfun(@(X) contains(X, 'contra'), compStringAllSubj);
plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
plotGroupDist = distanceAllSubj(iiGroup);
nCoSpikeCoDist = [];
nCoSpikeNoDist = [];
nCoSpikeCoDistSEM = [];
nCoSpikeNoDistSEM = [];
bins = quantile(plotGroupDist, 0:1/6:1); %:0.01:0.9;


for iB = 1:(length(bins)-1)
    ii = plotGroupDist >= bins(iB) & plotGroupDist < bins(iB+1);
    
    nCoSpikeCoDist = [nCoSpikeCoDist, mean(plotGroupCo(ii))];
    nCoSpikeNoDist = [nCoSpikeNoDist, mean(plotGroupNo(ii))];

    nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(plotGroupCo(ii), 1)];
    nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(plotGroupNo(ii),1)];


end

xVal = movmean(bins, 2);
xVal(1) = [];



clr = brewermap(10,'Paired');

p = polyfit(xVal, nCoSpikeCoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);
b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
b1.MarkerFaceColor = clr(8,:);
b1.Color = clr(8,:);
b1.LineWidth = 1.0;




p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);

b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
b2.MarkerFaceColor = clr(1,:);
b2.Color = clr(1,:);
b2.LineWidth = 1.0;
b2.MarkerSize = 5;

xlim([0 300])
ylim([-0.4 0.8])

ylabel('CoFire Event Rate [z-score]')
xlabel('fiber tract distance between units [mm]')
title('contra')

box off

subplot(1,3,3)
iiGroup = cellfun(@(X) ( contains(X, 'ACC') | ...
                       contains(X, 'SMA') | ...
                       contains(X, 'OFC') ) & ...
                      ~contains(X, 'HIP') & ...
                      ~contains(X, 'AMY')  , compStringAllSubj);
plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
plotGroupDist = distanceAllSubj(iiGroup);
nCoSpikeCoDist = [];
nCoSpikeNoDist = [];
nCoSpikeCoDistSEM = [];
nCoSpikeNoDistSEM = [];
bins = quantile(plotGroupDist, 0:1/6:1); %:0.01:0.9;


for iB = 1:(length(bins)-1)
    ii = plotGroupDist >= bins(iB) & plotGroupDist < bins(iB+1);
    
    nCoSpikeCoDist = [nCoSpikeCoDist, mean(plotGroupCo(ii))];
    nCoSpikeNoDist = [nCoSpikeNoDist, mean(plotGroupNo(ii))];

    nCoSpikeCoDistSEM = [nCoSpikeCoDistSEM, computeSEM(plotGroupCo(ii), 1)];
    nCoSpikeNoDistSEM = [nCoSpikeNoDistSEM, computeSEM(plotGroupNo(ii),1)];


end

xVal = movmean(bins, 2);
xVal(1) = [];



clr = brewermap(10,'Paired');

p = polyfit(xVal, nCoSpikeCoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);
b1 = errorbar(xVal, nCoSpikeCoDist, nCoSpikeCoDistSEM, 'o-'); hold on;
b1.MarkerFaceColor = clr(8,:);
b1.Color = clr(8,:);
b1.LineWidth = 1.0;




p = polyfit(xVal, nCoSpikeNoDist,1); %linear fit
x1 = linspace(min(xVal),max(xVal));
y1 = polyval(p,x1);

b2 = errorbar(xVal, nCoSpikeNoDist, nCoSpikeNoDistSEM, 'o--'); hold on;
b2.MarkerFaceColor = clr(1,:);
b2.Color = clr(1,:);
b2.LineWidth = 1.0;
b2.MarkerSize = 5;

xlim([0 300])
ylim([-0.4 0.8])

ylabel('CoFire Event Rate [z-score]')
xlabel('fiber tract distance between units [mm]')
title('cortico-cortico')


fig = gcf;
fig.Color = 'w';
box off

% savepdf(gcf, fullfile(pdfExport, 'coFireDistance_groups.pdf'))


figure('Position', [555 465 231 413]);
xCo = []; xCoSEM = [];
xNo = []; xNoSEM = [];
xCo(1) = mean(nCoSpikeCoAllSubj); xCoSEM(1) = std(nCoSpikeCoAllSubj)/sqrt(length(nCoSpikeCoAllSubj));
xNo(1) = mean(nCoSpikeNoAllSubj); xNoSEM(1) = std(nCoSpikeNoAllSubj)/sqrt(length(nCoSpikeNoAllSubj));

iiGroup = cellfun(@(X) ( contains(X, 'ACC') | ...
                       contains(X, 'SMA') | ...
                       contains(X, 'OFC') ) & ...
                      ~contains(X, 'HIP') & ...
                      ~contains(X, 'AMY')  , compStringAllSubj);

plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
xCo(4) = mean(plotGroupCo); xCoSEM(4) = std(plotGroupCo)/sqrt(length(plotGroupCo));
xNo(4) = mean(plotGroupNo); xNoSEM(4) = std(plotGroupNo)/sqrt(length(plotGroupNo));

iiGroup = cellfun(@(X) contains(X, 'ipsi'), compStringAllSubj);
plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
xCo(2) = mean(plotGroupCo); xCoSEM(2) = std(plotGroupCo)/sqrt(length(plotGroupCo));
xNo(2) = mean(plotGroupNo); xNoSEM(2) = std(plotGroupNo)/sqrt(length(plotGroupNo));

iiGroup = cellfun(@(X) contains(X, 'contra'), compStringAllSubj);
plotGroupCo = nCoSpikeCoAllSubj(iiGroup);
plotGroupNo = nCoSpikeNoAllSubj(iiGroup);
xCo(3) = mean(plotGroupCo); xCoSEM(3) = std(plotGroupCo)/sqrt(length(plotGroupCo));
xNo(3) = mean(plotGroupNo); xNoSEM(3) = std(plotGroupNo)/sqrt(length(plotGroupNo));






b1 = errorbar(xCo, xCoSEM, 'o'); hold on;
b1.MarkerFaceColor = clr(8,:);
b1.Color = clr(8,:);
b1.LineWidth = 1.0;


b2 = errorbar(xNo, xNoSEM, 'o'); hold on;
b2.MarkerFaceColor = clr(1,:);
b2.Color = clr(1,:);
b2.LineWidth = 1.0;
b2.MarkerSize = 5;

ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'All Sites', 'Ipsi', 'Contra', 'Cortico-cortical'};
ylim([-0.35 0.4])
xlim([0 5])
box off;

ylabel('CoFire Event Rate [z-score]')
fig = gcf;
fig.Color = 'w';
% savepdf(gcf, fullfile(pdfExport, 'coFireDistance_groupsConcat.pdf'))










