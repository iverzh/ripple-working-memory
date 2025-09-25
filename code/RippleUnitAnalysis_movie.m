
close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%
fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
flst = dir(fullfile(fldr, '*LFP*.mat'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;

% % 
subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P48CS_R1', 'P48CS_R2','P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};


study = 'Sternberg';
study = 'bmovie';
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
% dataDirectory =     sprintf('/space/seh10/6/halgdev/projects/iverzh/data/%s/preprocess/data_1kHz', study);
dataDirectory =     sprintf('/space/seh10/6/halgdev/projects/iverzh/data/%s/preprocess/', study);
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end

state = 'wake';
contactType = 'micro';

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfiles = dir(fullfile(dataDirectory, sprintf('*%s*', 'micro.mat')));
LFPfiles = {LFPfiles.name}';
taskFiles = dir(fullfile(dataDirectory, sprintf('*task*')));
taskFiles = {taskFiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../out', '*micro*'));
micrFiles = {micrFiles.name}';
timeShift = 0; %shift unit times to match LFP;
fs = 1e3;


channelCurationFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/assets';
badChan = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Bad Channels.csv'));
badChan = table2cell(badChan(3:end,3:end-1));   
bundleNotes = readtable(fullfile(channelCurationFolder, 'Ueli Movie Datasat Curation - Overall Bundle Notes.csv'));
bundleLab = table2cell(bundleNotes(2,3:end-1));
bundleNotes = table2cell(bundleNotes(3:end,3:end-1));   

recordingState = 'wake';
location = 'NC';
c = 0;
%%
% figure('Position', [1219 1160 1552 413]); 
% figure;
ccgPeaks = true;
saveLFPplots = false;
fldr =     sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/artifactSearch/%s', study);
pctOverlapAll = [];
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    units = LoadSpikeTimes(subject,'RutishauserLab', study);
    c = c + length(units);
%     for iU = 1:size(units,1)
%         subplot(1,2,1)
%         plot(cell2mat(units(iU,4))');
%         subplot(1,2,2)
%         times = units{iU,2}*1e3;
%         acg = times - times';
%         acg(acg == 0) = [];
%         histogram(acg(:), -(coRipPad+0.5):(coRipPad+0.5))
% 
%         waitforbuttonpress; clf;
%     end
    if ccgPeaks
        ccg = nan(size(units,1));
        pctOverlap = nan(size(units,1));
        Nab = cell(size(units,1));
        NautoA = cell(size(units,1));
        NautoB = cell(size(units,1));
        
        for a = 1:length(ccg)
            for b = a+1:length(ccg)
                if strcmp(study, 'bmovie')
                    locA = units{a,end}(1:4);
                    locB = units{b,end}(1:4);
                    
                elseif  strcmp(study, 'Sternberg')
                    locA = units{a,end};
                    locB = units{b,end};
                    
                end
                if a == b; continue; end
                if ~strcmp(locA, locB); continue; end
                timesA = units{a,2}*1e3;
                timesB = units{b,2}*1e3;
                ccgAB = timesA - timesB';
                
                autoA = timesA - timesA';
                autoA(autoA == 0) = nan;
                autoB = timesB - timesB';
                autoB(autoB == 0) = nan;
                
                [Nab{a,b}, timesCCG] = histcounts(ccgAB(:), -(20+0.25):0.5:(20+0.25));
                timesCCG = movmean(timesCCG,2); timesCCG(1) = [];
                baseline = mean(Nab{a,b}(abs(timesCCG) > 1 & abs(timesCCG) < 20));
                middle = mean(Nab{a,b}(abs(timesCCG) < 1));

                [NautoA{a,b}, ~] = histcounts(autoA(:), -(100+0.5):(100+0.5));
                [NautoB{a,b}, times] = histcounts(autoB(:), -(100+0.5):(100+0.5));
                times = movmean(times,2); times(1) = [];

                
                
                if ( length(timesA) > 200 ) && ( length(timesB) > 200)
                    Zab = zscore(Nab{a,b});
                    ccg(a,b) = max(Zab(abs(timesCCG) < 1));
                    pctOverlap(a,b) = sum(abs(ccgAB(:))<=1)/length(timesB) * 100;

                else
                    ccg(a,b) = nan; 
                end
                
                pctOverlapAll = [pctOverlapAll pctOverlap(a,b)]; 
                
                if ccg(a,b) > 5
                     
                    if strcmp(study, 'Sternberg')
                    [iiOvB, iiOvA] = find(abs(ccgAB) < 1);
                    iiA = 1:length(timesA);
                    iiA(ismember(iiA, iiOvA)) = [];
                    
                    iiB = 1:length(timesB);
                    iiB(ismember(iiB, iiOvB)) = [];

                    ccgAB(abs(ccgAB(:))<1) = nan;
                    [Nplot, ~] = histcounts(ccgAB(:),  -(100+0.5):(100+0.5));
 
                    figure('Position',[1310 118 965 805]);
                    subplot(2,2,1)
                    yyaxis left
                    plot(times, Nplot); 
                    yyaxis right
                    plot(timesCCG, Zab); 
                    title(sprintf('AB CCG  %02f %% <1ms', pctOverlap(a,b)))
                    subplot(2,2,2)
                    plot(times, NautoA{a,b}); 
                    title(sprintf('A auto unit %i (%i APs)',  a, length(timesA)))
                    subplot(2,2,3)
                    plot(times, NautoB{a,b}); 
                    title(sprintf('B auto unit %i (%i APs)',  b, length(timesB)))
                    
                    subplot(2,2,4)
                    wavA = mean(squeeze(units{a,4}(:,iiA,:)));
                    wavAov = mean(squeeze(units{a,4}(:,iiOvA,:)));
                    wavB = mean(squeeze(units{b,4}(:,iiB,:)));
                    wavBov = mean(squeeze(units{b,4}(:,iiOvB,:)));
                    plot(wavA, 'r-'); hold on;
                    plot(wavAov, 'r--'); hold on;
                    plot(wavB, 'b-');
                    plot(wavBov, 'b--'); hold on;

                    legend('A waveform','A peak wav', 'B waveform', 'B peak wav') 
                    
                    elseif  strcmp(study, 'bmovie')

% 
                        ccgAB(abs(ccgAB(:))<1) = nan;
                        [Nplot, ~] = histcounts(ccgAB(:),  -(100+0.5):(100+0.5));

                        figure('Position',[1310 118 965 805]);
                        subplot(2,2,1)
                        yyaxis left
                        plot(times, Nplot); 
                        yyaxis right
                        plot(timesCCG, Zab); 
                        title(sprintf('AB CCG  %02f %% <1ms', pctOverlap(a,b)))
                        subplot(2,2,2)
                        plot(times, NautoA{a,b}); 
                        title(sprintf('A auto unit %i (%i APs)',  a, length(timesA)))
                        subplot(2,2,3)
                        plot(times, NautoB{a,b}); 
                        title(sprintf('B auto unit %i (%i APs)',  b, length(timesB)))
                        %                     
                        subplot(2,2,4)
                        wavA = units{a,4};
                        wavB = units{b,4};
                        plot(wavA, 'r-'); hold on;
                        plot(wavB, 'b-');

                        legend('A waveform','B waveform')
                    
                    
                     end
 
% 
% 
                    sgtitle(subject)
                    savepdf(gcf, fullfile(fldr, sprintf('%s_unit%i_unit%i.pdf', subject, a, b)));
                    close

                end
            end
        end

        figure;
        imagesc(ccg, [0 5]); colorbar;
        title(subject)
        savepdf(gcf, sprintf('%s/%s_ccgPeaks.pdf', fldr, subject));
        
        figure;
        imagesc(pctOverlap, [0 10]); colorbar;
        title(subject)
        savepdf(gcf, sprintf('%s/%s_pctOverlap.pdf', fldr, subject));
        
        save(sprintf('%s/%s_ccgPeaks.mat', fldr, subject), 'ccg', 'pctOverlap', 'Nab', 'NautoA', 'NautoB', '-v7.3')
        
        pause(0.5); close all;
    end

    if saveLFPplots
    
        f = contains(LFPfiles, subject);
        micrObj = load(fullfile(dataDirectory, LFPfiles{f}));
        
         [b,a] = butter(3,[70 100]/(fs/2));
         data = micrObj.lfp_data;
         rippleBand = nan(size(data));
         rippleBandPhase = nan(size(data));
         for ch = 1:size(data,1)
             rippleBand(ch,:) = abs(hilbert(filtfilt(b,a,data(ch,:))));
    %          rippleBandPhase(ch,:) = angle(hilbert(rippleBand(ch,:))); 
             clc% zero phase
             fprintf('\n\n\n ... calculating ripple phase ...\n')
             fprintf('     %s [ %02i / %02i ]\n', subject, ch, size(data, 1))
         end
        
        figure('Position', [1 92 3008 1481])
        for ch = 1:size(micrObj.lfp_data,1)
            subplot(8,10, ch)
            plot(micrObj.lfp_data(ch,:))
    %         plot(rippleband(ch,:))
            
        end
        fldr =     '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/artifactSearch';
        savepdf(gcf, fullfile(fldr, sprintf('%s_microLFP.pdf', subject)))
    end
    
%     tag = [recordingState,'_',location,'_','IISremv_v3'];
%         filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, rippleContactType, tag);
%         load(fullfile(matExportFolder, filename));
%     close all
    
    
    
    
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
    units = LoadSpikeTimes(subject,'RutishauserLab', study);

%     tag = [recordingState,'_',location,'_',rippleC    ontactType];
%     filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    tag = [recordingState,'_',location,'_','IISremv_singleWire_v3'];
    tag = [recordingState,'_',location,'_','1kHz_clip_z25'];
%     filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, rippleContactType, tag);
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    load(fullfile(matExportFolder, filename));
    chan_labels = str2double(rippleStats.chanLabels);
%     f = contains(micrFiles, subject);
%     micrObj = load(fullfile(dataDirectory, '../out', micrFiles{f}));

    f = contains(LFPfiles, subject);
    micrObj = load(fullfile(dataDirectory, LFPfiles{f}));
    if strcmp(study, 'Sternberg'); micrObj.times = 0:1/1e3:length(micrObj.data)/1e3; end
    f = contains(bpFiles, subject);
%     bpObj = load(fullfile(dataDirectory, '../out', bpFiles{f}));
    
    
    
    
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
                chanLabel = units{ui,1};% - 96;
                chanLabel = find(ismember(chan_labels,chanLabel));
                if isempty(chanLabel); chanLabel = 1; end

                if strcmp(contactType, 'macro')
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
    
        if strcmp(contactType, 'micro'); dataChan = micrObj.data(chanLabel, :);
        elseif strcmp(contactType, 'macro'); dataChan = bpObj.data(:,chanLabel)';
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
                Npyr(uc,:) = zscore(Npyr(uc,:));
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

win = 50;
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
spikingWindow = N_Hz(:, 601-win:601+win);
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

ax = gca;
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
filename = sprintf('%s_%s_%s_ripple-unit_noAxis.pdf', subject, state, contactType);
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

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
fs = 1e3;
cd(exportDir); pause(0.5);
nUnit = 0;
nUnitNoLFP = 0;
for v = 3
    cUnit = 0; 
    saveLFP = true;
    padLength = 0; %sec
    for coRipPadMS = 1000 %[1000, 200, 100] %0:50:200
        coRipPad = round(coRipPadMS * fs / 1e3);
        unitCCGstep = fs / 1e3;
        histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
        for subj = length(subj_list_full):-1:1
            subject = subj_list_full{subj};
            units = LoadSpikeTimes(subject,'RutishauserLab', study);
        
    %         tag = [recordingState,'_',location,'_',rippleContactType];
    %         filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
%             tag = sprintf([recordingState,'_',location,'_','IISremv_singleWire_v%i'], v);
%             filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, contactType, tag);
            
%             tag = sprintf([recordingState,'_',location,'_','IISremv_singleWire_v%i'], v);
%             filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, contactType, tag);

            tag = sprintf([recordingState,'_',location,'_','1kHz_template_z25'], v);
            filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
            load(fullfile(matExportFolder, filename));
            chan_labels = rippleStats.chanLabels;
    
            if strcmp(study, 'Sternberg')
                chan_labels = cellfun(@(X) str2double(X), chan_labels);
            end
        % %     f = contains(micrFiles, subject);
        % %     micrObj = load(fullfile(dataDirectory, '../out', micrFiles{f}));
        % 
            f = contains(LFPfiles, subject);
            micrObj = load(fullfile(dataDirectory, LFPfiles{f}));
            if strcmp(study, 'Sternberg'); 
                micrObj.times = 0:1/1e3:length(micrObj.data)/1e3;
                micrObj.lfp_data = micrObj.data;
            end
            timesCCG = micrObj.times * rippleStats.fs;
        %     
        %     f = contains(bpFiles, subject);
        %     bpObj = load(fullfile(dataDirectory, '../out', bpFiles{f}));
            
            rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
            for chRipp = 1:size(rippMask,1) 
                if ~isempty(rippleStats.window{chRipp}) 
                    iS = round(timesCCG(rippleStats.window{chRipp}(:,1)));
                    iE = round(timesCCG(rippleStats.window{chRipp}(:,2)));
    
    %                 iS = rippleStats.window{chRipp}(:,1);
    %                 iE = rippleStats.window{chRipp}(:,2);
        
                
                    for ii = 1:length(iE)            
                        if iS(ii) <= 0 || iE(ii) <= 0; continue; end
        
                        rippMask(chRipp,iS(ii):iE(ii)) = 1;
                    end
                end
                
            end
        
            [b,a] = butter(3,[70 100]/(rippleStats.fs/2));
            data = micrObj.lfp_data;
    
            data = data(:, (padLength*fs+1):end);
            rippleBand = nan(size(data));
            rippleBandPhase = nan(size(data));
            for ch = 1:size(data,1)
                for nf = 60:60:fs/2 %240
                    Wo = nf/(fs/2);
                    BW = Wo/35;
                    [bN,aN] = iirnotch(Wo, BW);
                    data(ch,:) = filtfilt(bN,aN,data(ch,:));
                end
    %             rippleBand(ch,:) = filtfilt(b,a,data(ch,(padLength*fs):end));
                rippleBand(ch,:) = filtfilt(b,a,data(ch,:));
                rippleBandPhase(ch,:) = angle(hilbert(rippleBand(ch,:))); 
                clc% zero phase
                fprintf('\n\n\n ... calculating ripple phase ...\n')
                fprintf('     %s [ %02i / %02i ]\n', subject, ch, size(data, 1))
            end
        
        
            coRipSpike = cell(length(units),length(units));
            coRipDur = cell(length(units),length(units));
            controlDur = cell(length(units),length(units));
            noRipDur = cell(length(units),length(units));
            noRipAspike = cell(length(units),length(units));
            noRipSpike = cell(length(units),length(units));
            interactionType = cell(length(units),length(units));
            coFireProbNoR = cell(length(units),length(units));
            ccgNoR = cell(length(units),length(units),2);
            ccgCoR = cell(length(units),length(units));
            ccgNoRcontrol = cell(length(units),length(units));
            coSpikeTimes = cell(length(units),length(units));
            coSpikePhase = cell(length(units),length(units));
            coSpikeLFP = cell(length(units),length(units));
            LFPall = cell(size(rippMask,1),size(rippMask,1), 2);
            rbAll = cell(size(rippMask,1),size(rippMask,1), 2);
            for uAi = 1:size(units, 1)
                if strcmp(study, 'Sternberg')
                    nUnit = nUnit + 1;

                    chA = find(chan_labels == units{uAi,1});
                    if isempty(chA); nUnitNoLFP = nUnitNoLFP + 1; end
                end
                
                
                for uBi = 1:size(units, 1)
                    uA = uAi; uB = uBi;
                    if strcmp(study, 'Sternberg')
                        chA = find(chan_labels == units{uA,1}); chB = find(chan_labels == units{uB,1});
                    else
                        chA = find(strcmp(chan_labels, units{uA,end})); chB = find(strcmp(chan_labels, units{uB,end}));
                    end
    
                    typeA = units{uA,3}; typeB = units{uB,3};
                    if isempty(chA) || isempty(chB); continue; end
                    
                    if strcmp(study, 'bmovie')
                        chLocA = regexprep(rippleStats.chanLabels{chA}, '\d', '');
                        chNumA = regexp(rippleStats.chanLabels{chA}, '\d+', 'match');
                        chNumA = str2double(chNumA);
                        bundleChA = find(contains(bundleLab, chLocA));

                        chLocB = regexprep(rippleStats.chanLabels{chB}, '\d', '');
                        chNumB = regexp(rippleStats.chanLabels{chB}, '\d+', 'match');
                        chNumB = str2double(chNumB);
                        bundleChB = find(contains(bundleLab, chLocB));

                        bchA = badChan(subj, bundleChA);
                        bchA = str2double(split(bchA,','));
                        bchB = badChan(subj, bundleChB);
                        bchB = str2double(split(bchB,','));
                    end
                    

    
                    if chA ~= chB && any(strcmp(typeA,{'pyr', 'int'})) && any(strcmp(typeB,{'pyr', 'int'}))
%                         if (contains(bundleNotes(subj, bundleChA), 'Bad') || any(ismember(bchA, chNumA))) || ...
%                            (contains(bundleNotes(subj, bundleChB), 'Bad') || any(ismember(bchB, chNumB)))
% %                             rippMask(chRipp,:) = nan(1,length(rippMask));
%                             continue
%                         end
                        uTa = (units{uA,2}) * rippleStats.fs; uTa = uTa(uTa >= 0.5);
                        uTb = (units{uB,2}) * rippleStats.fs; uTb = uTb(uTb >= 0.5);
    
    %                     uTa = round(uTa); uTb = round(uTb);
            
                        coRipMask = rippMask(chA, :) & rippMask(chB, :);
                        noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :);
            
                        
            
                        uTcoRa = uTa(coRipMask(round(uTa)));
                        uTcoRb = uTb(coRipMask(round(uTb)));
            
                        uTnoRa = uTa(noRipMask(round(uTa)));
                        uTnoRb = uTb(noRipMask(round(uTb)));
            
                        controlMask = zeros(size(coRipMask));
            
                        if ~isempty(uTcoRa) && ~isempty(uTcoRb)
                            
                            
                            bnds = mask2bounds(coRipMask);
                            dur = [bnds(:,2) - bnds(:,1)] / rippleStats.fs * 1e3;
                            bnds(dur < 25, :) = [];
                            dur(dur < 25) = [];
        
                            if size(bnds,1) < 2; continue; end
                            
                            coSpikeCount = nan(2,length(bnds));
                            coSpikeTimesAB = cell(7,length(bnds));
                            coSpikePhaseAB = cell(1,length(bnds));
                            coSpikeLFPAB = cell(1,length(uTcoRa)); cnt = 1;
                            noSpikeCount = nan(2,length(bnds));
                            LFPA = nan(length(bnds), 201);
                            LFPB = nan(length(bnds), 201);
                            rbA = nan(length(bnds), 201);
                            rbB = nan(length(bnds), 201);
                            bndsBaselineAll = nan(size(bnds)); 
                            tme = 0;
                            for ibnd  = 1:length(bnds)
                                aSpikesCo = sum(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                                bSpikesCo = sum(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));
            
                                shift = round(1e3 * rippleStats.fs / 1e3);
                                bndsBaseline = bnds(ibnd,:) - shift;
                                if any(bndsBaseline <= 0); continue; end
            
                                loopCount = 0;
                                while ~all(noRipMask(bndsBaseline(1):bndsBaseline(2))) && loopCount < 1000  
                                    bndsBaseline = bndsBaseline - 100;
                                    loopCount = loopCount + 1;
                                    if any(bndsBaseline < 1)
                                        loopCount = 1000;
                                        break
                                    end
                                end
            
%                                 if loopCount > 1000; continue; end

%                                 
            
                                aSpikesNo = sum(uTnoRa >= bndsBaseline(1)   &  uTnoRa <= bndsBaseline(2));
                                bSpikesNo = sum(uTnoRb >= bndsBaseline(1)   &  uTnoRb <= bndsBaseline(2));
                                
                                uTcoRaBnd = uTa(uTa >= bnds(ibnd,1)  &  uTa <= bnds(ibnd,2));
                                coSpikeCount(:, ibnd)     = [aSpikesCo, bSpikesCo];
                                coSpikeTimesAB{1, ibnd}   = uTa(uTa >= (bnds(ibnd,1)-coRipPad)  &  uTa <= (bnds(ibnd,2)+coRipPad));
                                coSpikeTimesAB{2, ibnd}   = uTcoRaBnd;
                                coSpikeTimesAB{3, ibnd}   = uTb(uTb >= (bnds(ibnd,1)-coRipPad)  &  uTb <= (bnds(ibnd,2)+coRipPad));
                                
                                uTcontrolRaBnd = uTa(uTa >= bndsBaseline(1)  &  uTa <=  bndsBaseline(2) );
                                coSpikeTimesAB{4, ibnd}   = uTcontrolRaBnd;
                                coSpikeTimesAB{5, ibnd}   = uTb(uTb >= ( bndsBaseline(1) -coRipPad)  &  uTb <= (bndsBaseline(2)+coRipPad));
                                
                                coSpikeTimesAB{6, ibnd}   = rippleBandPhase(chA,round(coSpikeTimesAB{2, ibnd}));
                                coSpikeTimesAB{7, ibnd}   = rippleBandPhase(chB,round(coSpikeTimesAB{3, ibnd}));
    
                                coSpikePhaseAB{ibnd}(1,:) = rippleBandPhase(chA,bnds(ibnd,1) :bnds(ibnd,2) );
                                coSpikePhaseAB{ibnd}(2,:) = rippleBandPhase(chB,bnds(ibnd,1) :bnds(ibnd,2) );
        
                                noSpikeCount(:, ibnd) = [aSpikesNo bSpikesNo];
                                bndsBaselineAll(ibnd, :) = bndsBaseline;
    
                                if saveLFP && (aSpikesCo > 0) && (bSpikesCo > 0)
                                    for iSp = 1:length(uTcoRaBnd)
                                        if (uTcoRaBnd(iSp) < 201) || (uTcoRaBnd(iSp) + 200 > length(data)); continue;
                                        end
                                        spk = round(uTcoRaBnd(iSp));
                                        coSpikeLFPAB{cnt}(1,:) = [ibnd data(chA,spk-200 :spk+200 )];
                                        coSpikeLFPAB{cnt}(2,:) = [ibnd data(chB,spk-200 :spk+200 )];
    
                                        coSpikeLFPAB{cnt}(3,:) = [ibnd rippleBand(chA,spk-200 :spk+200 )];
                                        coSpikeLFPAB{cnt}(4,:) = [ibnd rippleBand(chB,spk-200 :spk+200 )];
    
                                        cnt = cnt + 1;
        
                                    end
                                end
                                center = round(mean(bnds(ibnd, :)));
                                LFPA(ibnd,:) = data(chA,center-100 :center+100 );
                                LFPB(ibnd,:) = data(chB,center-100 :center+100 );
                                rbA(ibnd,:)  = rippleBand(chA,center-100 :center+100 );
                                rbB(ibnd,:)  = rippleBand(chB,center-100 :center+100 );
                                
            
                            end
                            
                            LFPall{chA,chB,1} = LFPA;
                            LFPall{chA,chB,2} = LFPB;
                            rbAll{chA,chB,1} = rbA;
                            rbAll{chA,chB,2} = rbB;
                            
                            
            
                            controlMask = bounds2mask(bndsBaselineAll, length(coRipMask));
                            uTcontrola = uTa(controlMask(round(uTa)));
                            controlMask = bounds2mask(bndsBaselineAll, length(coRipMask), coRipPad);
                            uTcontrolb = uTb(controlMask(round(uTb)));
                
                            s = uTcontrola - uTcontrolb';
                            Ncontrol = histcounts(s, histEdges);
            
                       
                            coRipSpike{uA,uB} = coSpikeCount;
                            coRipDur{uA,uB} = dur;
                            controlDur{uA,uB} = bndsBaselineAll(:,2) - bndsBaselineAll(:,1);
                            noRipDur{uA,uB} = sum(noRipMask);
                            noRipAspike{uA,uB} = length(uTnoRa);
                            noRipSpike{uA,uB} = noSpikeCount;
    %                         coFireProbNoR{uA,uB} = N;
    
                            coRipMask = bounds2mask(bnds, length(coRipMask));
                            uTcoRa = uTa(coRipMask(round(uTa)));
                            coRipMask = bounds2mask(bnds, length(coRipMask), coRipPad);
                            uTcoRb = uTb(coRipMask(round(uTb)));
                            s = uTcoRa - uTcoRb';
                            Nab = histcounts(s, histEdges);
                            ccgCoR{uA,uB} = Nab;
        
                            coSpikeTimes{uA, uB} = coSpikeTimesAB;
                            coSpikePhase{uA, uB} = coSpikePhaseAB;
                            coSpikeLFP{uA, uB}   = coSpikeLFPAB;
        
                            s = uTnoRa - uTnoRb';
                            Nab = histcounts(s, histEdges);
                            ccgNoR{uA,uB,1} = uTnoRa;
                            ccgNoR{uA,uB,2} = uTnoRb;
                            ccgNoRcontrol{uA,uB} = Ncontrol;
                            
                            
                        end
                    end
                end
                clc
                fprintf('... %s co firing unit %i / %i ... \n co ripple padding %03i\n', subject, uAi, length(units), coRipPad)
                 cUnit=cUnit+1;
            end
%             
            nCoSpikeCoAll = [];
            nCoSpikeNoAll = [];
            
            distanceAll = [];
            c = 1;
            
            nCoSpikeNo = [];
            nCoSpikeCo = [];
            %     dur = 0;
            for uAi = 1:size(units,1)
                for uBi = 1:size(units,1)
            %         uA = units(uAi); uB = units(uBi);
                    uA = uAi; uB = uBi;
    %                 chA = find(chan_labels == units{uA,1}); chB = find(chan_labels == units{uB,1});
                    chA = find(strcmp(chan_labels, units{uA,end})); chB = find(strcmp(chan_labels, units{uB,end}));
    
                    if isempty(chA) || isempty(chB); continue; end
                    d = pdist([micrObj.chan_coords(chA,:); micrObj.chan_coords(chB,:)]);
        
                    
                    coRip = coRipSpike{uA,uB};
                    noRip = noRipSpike{uA,uB};
                    dur = sum(coRipDur{uA,uB}) / rippleStats.fs;
        
                    if ~isempty(coRip) 
                        sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                        nCoSpikeCo = size(sp,2); %/dur; % / (dur);
                        
                        sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                        nCoSpikeNo = size(sp,2); %/dur; % / (dur);
        %                     dur = dur + (sum(coRipDur{uA,uB}) / rippleStats.fs);
            
                        nCoSpikeCoAll = [nCoSpikeCoAll, nCoSpikeCo];
                        nCoSpikeNoAll = [nCoSpikeNoAll, nCoSpikeNo];
                        distanceAll = [distanceAll, d];
                    end
                    
                        
            
                end
            end
            uCh = units(:,1);
            uLoc = units(:,end); %cellfun(@(X,Y) [X '_' num2str(Y)], units(:,end), units(:,1), 'UniformOutput', false);
            filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_1kHz_template.mat', subject, coRipPad);
            save(fullfile(exportDir, filename), 'nCoSpikeCoAll', 'nCoSpikeNoAll', 'distanceAll', ...
                                                'coRipSpike', 'coRipDur', 'controlDur', 'noRipSpike', ...
                                                'noRipDur', 'noRipAspike', 'coFireProbNoR', ...
                                                'ccgNoR', 'ccgCoR','ccgNoRcontrol', 'uCh', 'uLoc',...
                                                'coSpikeTimes', 'coSpikePhase', 'coSpikeLFP', 'LFPall', 'rbAll', '-v7.3')
            
        
            
        end
    end
end
%%

exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles';
pdfExport = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));
fs = 1e3; %1e3, 400;
ctxParc = {'ACC', 'SMA', 'OFC' ,'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};
for v = 6
    for ang = 2 %2:5 %2:6
        tic
        
        
        
        
               
        phLagsAll = nan(1, 1e6); cLags = 1;
        cAP =1;
        coRipPad = 1 * fs;
        durTrough = 0; durPeak = 0;
        unitCCGstep = fs / 1e3;
        histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
        coFd = 20;
        minDur = 50;
        
        
        lfpA =  nan(1e6, length(-100:100)); lfpB =  nan(1e6, length(-100:100)); 
        rbA =  nan(1e6, length(-100:100)); rbB = nan(1e6, length(-100:100));
        compStringAllSubj = []; phLag =  nan(1,1e6);
        nCoSpikeCoAllSubj = []; nCoSpikeNoAllSubj =  [];
        distanceAllSubj =[]; 
        ccgAllcoR = nan(1e6, length(histEdges)-1);
        ccgAllcoRpeak =  nan(1e6, length(histEdges)-1);
        ccgAllcoRtrough =  nan(1e6, length(histEdges)-1);
        ccgAllnoRAll = nan(1e6, 2001);
        ccgAllnoRepoch =  nan(1e6, length(histEdges)-1);
        durAllcoR =  nan(1,1e6);
        durAllcoRpeak =  nan(1,1e6);
        durAllcoRtrough =  nan(1,1e6);
        durAllnoRAll =  nan(1,1e6);
        durAllnoRepoch =  nan(1,1e6);
        nAspikeAllcoR =  nan(1,1e6);
        nAspikeAllcoRpeak =  nan(1,1e6);
        nAspikeAllcoRtrough =  nan(1,1e6);
        nAspikeAllnoRAll =  nan(1,1e6);
        nAspikeAllnoRepoch =  nan(1,1e6);
        compStringAll_LFP = cell(1,1e6);
        compStringAll_phase = cell(1,1e6);
        aSpPhase = nan(1,1e6);
        bSpPhase = nan(1,1e6);
        Nphase = zeros(1, length(histEdges)-1);
        c = []; countLFP = 0;countPhase = 0;countCCG = 1;
        for subj = 1:length(subj_list_full) %[1:17 20:length(subj_list_full)] % [1:15, 18:length(subj_list_full)]
            subject = subj_list_full{subj};
            nCo = []; nNo = []; durationAll = [];
%             if ~contains(subject, 'CS'); continue; end
            filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_IISremv_SU_v%i.mat', subject,fs, v);
%             filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_1kHz_template.mat', subject,fs);
            load(fullfile(exportDir, filename)); %, 'nCoSpikeCoAll', 'nCoSpikeNoAll', 'distanceAll', ...
    %                                             'coRipSpike', 'coRipDur', 'controlDur', 'noRipSpike', 'noRipDur', 'coFireProbNoR', ...
    %                                             'ccgNoR', 'ccgCoR','ccgNoRcontrol', 'uCh', ...
    %                                             'coSpikeTimes','coSpikePhase'); %, 'coSpikeLFP',)
            f = contains(LFPfiles, subject);
            micrObj = load(fullfile(dataDirectory, LFPfiles{f}), 'chan_labels');
            
            tag = [recordingState,'_',location,'_','1kHz_template_z25'];

            filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
            load(fullfile(matExportFolder, filename));
            chan_labels = str2double(rippleStats.chanLabels);
            
            if contains(subject, 'sub'); study = 'Sternberg';
            else; study = 'bmovie'; end
    
            
            if strcmp(study, 'Sternberg')
                splitPattern = '(_right|_left)';
                locations = regexp(uLoc', splitPattern, 'split')';
                locations = cellfun(@(X) X{:}, locations,  'UniformOutput', false);
                locations(cellfun(@(X) isempty(X), locations)) = [];
                locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
                locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
                locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
                locations = strrep(locations, 'amygdala', 'AMY');
                locations = strrep(locations, 'hippocampus', 'HIP');


                hem = regexp(uLoc, splitPattern, 'match')';
                hem = cellfun(@(X) X{:}, hem',  'UniformOutput', false);

                hemi = contains(hem, 'left');
                locations(hemi) = cellfun(@(X) ['L' X], locations(hemi),  'UniformOutput', false);
                hemi = contains(hem, 'right');
                locations(hemi) = cellfun(@(X) ['R' X], locations(hemi),  'UniformOutput', false);
                
                task = load(fullfile(dataDirectory,'../task/', sprintf('%s_task.mat', subject)));
                load1acc = sum(task.response_accuracy(task.loads == 1))/ sum(task.loads == 1);
                load3acc = sum(task.response_accuracy(task.loads == 3))/ sum(task.loads == 3);
%                 fprintf('load 1 performance %.2f    load 3 performance %.2f.\n', load1acc, load3acc);
                
%                 if load1acc < 0.9 || load3acc < 0.8; fprintf('skipping %s due to poor task performance.\n', subject); continue; end
            end

            tractDistance =  []; compString = [];
            phLagsSubj = []; 

            cSubj = 0;
             for uA = 1:size(coRipDur, 1)
                for uB = 1:size(coRipDur, 2)



                    coRip = coRipSpike{uA,uB};
                    noRip = noRipSpike{uA,uB};
                    dur = (sum(coRipDur{uA,uB})); %+ 2*0*length(coRipDur{uA,uB})) /fs;
                    durNo = noRipDur{uA,uB} / fs;

%                     if (~isempty(coRip) || ~isempty(noRip)) && uCh{uA} == uCh{uB} %&& strcmp(interactionType{uA, uB}, interact) % || strcmp(interactionType{uA, uB}, 'int-pyr'))
                    if dur > 0 && uCh{uA} ~= uCh{uB} %&& strcmp(interactionType{uA, uB}, interact) % || strcmp(interactionType{uA, uB}, 'int-pyr'))
                        sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                        nCoSpikeCo = size(sp,2); 

                        sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                        nCoSpikeNo = size(sp,2);

                        nCo = [nCo, nCoSpikeCo];
                        nNo = [nNo, nCoSpikeNo];
                        durationAll = [durationAll, dur];

                        if strcmp(study, 'bmovie')
                            parcelA = uLoc{uA}(2:end-1);
                            parcelB = uLoc{uB}(2:end-1);
                            hemA = uLoc{uA}(1);
                            hemB = uLoc{uB}(1);

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
                        else
%                             d = -1;

                            hemA = locations{uA}(1);
                            hemB = locations{uB}(1);
                            if strcmp(hemA, hemB); hem = 'ipsi';
                            else; hem = 'contra'; end
                            
                            if strcmp(locations{uA},  locations{uB}); d = 0; 
                            else; d = 1; end

                            compString = [compString {sprintf('%s %s %s', locations{uA}(2:end),locations{uB}(2:end), hem)}];
%                             compString = [compString; {locations{uA}, locations{uB}}];
                            tractDistance = [tractDistance d];
                            
                        end


        %                 phLag =


                        if d >= 0 %&& ~contains(compString(end), {'AMY', 'HIP'}) 
                            phMat = coSpikePhase{uA, uB};
                            lfpMat = coSpikeLFP{uA, uB};
                            lfpMat(cellfun(@(X) isempty(X), lfpMat)) = [];
                            if ~isempty(lfpMat)
                                iCoRip = cell2mat(cellfun(@(X) X(1,1), lfpMat, 'UniformOutput', false) );

                                lfpMatA = cell2mat(cellfun(@(X) X(1,2:end)', lfpMat, 'UniformOutput', false) );
                                lfpMatB = cell2mat(cellfun(@(X) X(2,2:end)', lfpMat, 'UniformOutput', false) );
                                rbMatA = cell2mat(cellfun(@(X) X(3,2:end)', lfpMat, 'UniformOutput', false) );
                                rbMatB = cell2mat(cellfun(@(X) X(4,2:end)', lfpMat, 'UniformOutput', false) );
                                cSubj = cSubj + size(lfpMatA,2);

                            else
                                iCoRip = [];
                                lfpMatA = [];
                                lfpMatB = [];
                                rbMatA = [];
                                rbMatB = [];
                            end
                            emptyPh = cellfun(@(X) isempty(X), phMat);
                            phMat(emptyPh) = {zeros(2,50)};
                            phLags = cellfun(@(X) circ_mean([X(1,:) - X(2,:)]'), phMat);
                            phLags(emptyPh) = nan;


                            zeroLag = (phLags >= - pi/ang & phLags <= pi/ang) & coRipDur{uA,uB}' >= minDur;
                            nzeroLag = (phLags < - (ang-1)*pi/ang | phLags > (ang-1)*pi/ang)  & coRipDur{uA,uB}' >= minDur;
                            allSpike = coRipDur{uA,uB}' >= minDur;
    %                         nzeroLag = (phLags < - pi/ang | phLags > pi/ang)  & coRipDur{uA,uB}' >= 50;
        %                     nzeroLag = (phLags > pi/6 | phLags < -pi/6)  & coRipDur{uA,uB}' >= 50;
        %                     nzeroLag = phLags <= - pi/2 & phLags >= -pi  & coRipDur{uA,uB}' >= 50;
                            allLag = phLags((zeroLag | nzeroLag) & (coRip(1,:) > 0 & coRip(2,:) > 0));
%                             allLag = [phLags(coRip(1,:) > 0 & coRip(2,:) > 0)];
                            phLagsSubj = [phLagsSubj allLag];
                            phLagsAll(cLags:cLags+length(allLag)-1) = allLag;
                            cLags = cLags + length(allLag);

                            spMat = coSpikeTimes{uA, uB};
        %                     ;
                            uTa = cell2mat(spMat(2,zeroLag)); uTb = cell2mat(spMat(3,zeroLag));
                            if ~isempty(uTa) && ~isempty(uTb)
                                s = uTa - uTb';
                                cAP = cAP + length(uTa);
                                Nz = histcounts(s, histEdges);
                            else
                                Nz = zeros(1,length(histEdges)-1);
                            end
                            nAspikeAllcoRpeak(countCCG) = length(uTa);


        %                     durPeak = [durPeak + sum(coRipDur{uA,uB}(zeroLag)) / fs];
                            durAll = sum(coRipDur{uA,uB}) / fs;
                            durPeak = sum(coRipDur{uA,uB}(zeroLag)) / fs;
    %                         Nz = Nz/durPeak;
                            Nnorm = Nz; %/sum(Nz); durPeak; %sum(Nz);
%                             Nnorm = Nz/durPeak; %sum(Nz);
                            ccgAllcoRpeak(countCCG, :) = Nnorm;
                            durAllcoRpeak(countCCG) = durPeak;

                            uTa = cell2mat(spMat(2,nzeroLag)); uTb = cell2mat(spMat(3,nzeroLag));
                            if ~isempty(uTa) && ~isempty(uTb)
                                s = uTa - uTb';
                                Nnz = histcounts(s, histEdges);
                            else
                                Nnz = zeros(1,length(histEdges)-1);
                            end
                            nAspikeAllcoRtrough(countCCG) = length(uTa);
                            nAspikeAllcoR(countCCG) =  nAspikeAllcoRtrough(countCCG) + nAspikeAllcoRpeak(countCCG);


                            uTa = cell2mat(spMat(4,:)); uTb = cell2mat(spMat(5,:));
                            if ~isempty(uTa) && ~isempty(uTb)
                                s = uTa - uTb';
                                Ncontrol = histcounts(s, histEdges);
                            else
                                Ncontrol = zeros(1,length(histEdges)-1);
                            end
                            nAspikeAllnoRepoch(countCCG) = length(uTa);



        %                     durTrough = [durTrough + sum(coRipDur{uA,uB}(nzeroLag)) / fs];
                            durTrough = sum(coRipDur{uA,uB}(nzeroLag)) / fs;
                            Nnorm = Nnz; %/sum(Nnz); %durTrough; %sum(Nnz);
%                             Nnorm = Nnz/durTrough; %durTrough; %sum(Nnz);
                            ccgAllcoRtrough(countCCG, :) = Nnorm;
                            durAllcoRtrough(countCCG) = durTrough;

                            Nall = Nz + Nnz; %ccgCoR{uA, uB};
        %                     N = N/dur;
        %                     Nnorm = Nall/sum(Nall); %(durTrough+durPeak); %sum(Nall);
%                             Nnorm = Nall/(durTrough+durPeak); %sum(Nall);
        %                     Nnorm = Nall/dur; %sum(Nall);
                            Nnorm = Nall; %/sum(Nall);
                            ccgAllcoR(countCCG, :) = Nnorm;
                            durAllcoR(countCCG) = (durTrough+durPeak);

                            uTnoRa = ccgNoR{uA,uB,1};
                            uTnoRb = ccgNoR{uA,uB,2};
                            s = uTnoRa - uTnoRb';
                            Nab = histcounts(s, histEdges);
                            

                            ccgAllnoRAll(countCCG, :) = Nab;
                            durAllnoRAll(countCCG) = durNo;
                            nAspikeAllnoRAll(countCCG) = length(uTnoRa);
                            
                            ccgAllnoRepoch(countCCG, :) = Ncontrol;
                            durAllnoRepoch(countCCG) = durAll;

%                             iiZ = ismember(iCoRip, find(zeroLag));
%                             iiZ = true(1, length(iCoRip));

%                             lfpA(countLFP+1:countLFP + sum(iiZ),:) = lfpMatA(:,iiZ)';
%                             lfpB(countLFP+1:countLFP + sum(iiZ),:) = lfpMatB(:,iiZ)';
%                             rbA(countLFP+1:countLFP + sum(iiZ),:)  = rbMatA(:,iiZ)';
%                             rbB(countLFP+1:countLFP + sum(iiZ),:)  = rbMatB(:,iiZ)';
                            
                            chA = find(micrObj.chan_labels == uCh{uA}); chB = find(micrObj.chan_labels == uCh{uB});
%                             if ~isempty(LFPall{chA,chB,1})
%                                 iiZ = size(LFPall{chA,chB,1}, 1);
%                                 lfpA(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,1};
%                                 lfpB(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,2};
%                                 rbA(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,1};
%                                 rbB(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,2};
%                                 
%                                 compStringAll_LFP(countLFP+1:countLFP + iiZ) = repmat(compString(end), [1 iiZ]);
% 
% 
%                                 countLFP = countLFP + iiZ;
% 
%                             end
                            
                            if d > 0
                                uTa = cell2mat(spMat(2,allSpike)); uTb = cell2mat(spMat(3,allSpike));
                                uPha = cell2mat(spMat(6,allSpike)); uPhb = cell2mat(spMat(7,allSpike));
                                s = (uTa - uTb');
%                                 [b, a] = find(abs(s) <= coFd & abs(s) >= 5);
                                [b, a] = find(s <= coFd & s >= 5);
                                if ~isempty(a)
                                    aii = find(isnan(aSpPhase), 1, 'first');
                                    aSpPhase(aii:aii+(length(a)-1)) = uPha(a);
                                    bii = find(isnan(bSpPhase), 1, 'first');
                                    bSpPhase(bii:bii+(length(b)-1)) = uPhb(b);
                                    
                                    Nphase = Nphase + histcounts(s(b,a), histEdges);

                                end
                                
                                
                            end
                           
                            
                            
                            compStringAll_phase(countPhase+1:countPhase + length(allLag)) = repmat(compString(end), [1 length(allLag)]);

%                             countLFP = countLFP + sum(iiZ);
                            countPhase = countPhase + length(allLag);
                            countCCG = countCCG + 1;

    %                         figure; 
    %                         subplot(2,1,1)
    %                         plot(rbMatA(:,iiZ));
    %                         subplot(2,1,2)
    %                         plot(rbMatB(:,iiZ));

    %                         waitforbuttonpress; close;

                        end



                    end
                end
             end
             
%              uCh = cell2mat(uCh);
%              for iA = 1:length(uCh)
%                  for iB = chA+1:length(uCh)
%                      chA = uCh(iA); chB = uCh(iB); 
%                     if ~isempty(LFPall{chA,chB,1})
%                         iiZ = size(LFPall{chA,chB,1}, 1);
%                         lfpA(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,1};
%                         lfpB(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,2};
%                         rbA(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,1};
%                         rbB(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,2};
% 
%                         countLFP = countLFP + iiZ;
% 
%                     end
%                  end
%              end
             
             

    %          figure;
    % %     plot(-100:100, mean(lfpA([1:1.5e4],:)), 'r'); hold on;
    % %     plot(-100:100, mean(lfpB([1:1.5e4],:)), 'b'); hold on;
    % %     plot(-100:100, mean(lfpA(2.5e4:end,:)), 'r'); hold on;
        % %     plot(-100:100, mean(lfpB(2.5e4:end,:)), 'b'); hold on;
        %         plot(-100:100, mean(lfpA), 'r'); hold on;
        %         plot(-100:100, mean(lfpB), 'b'); hold on;
        %         vline(-99:11:99)
        %         xlim([-50 50])
        %         title(subject)
        % 
        %         waitforbuttonpress; close;

                 c = [c cSubj];

            %     zAll = zscore([nCoSpikeCoAll(distanceAll > 0) nCoSpikeNoAll(distanceAll > 0)]);
            %     nCoSpikeCoAll = zAll(1:sum(distanceAll > 0));
            %     nCoSpikeNoAll = zAll(sum(distanceAll > 0)+1:end);

%                 zAll = zscore([nCoSpikeCoAll(tractDistance >= 0) nCoSpikeNoAll(tractDistance >= 0)]);
%                 nCoSpikeCoAll = zAll(1:sum(tractDistance >= 0));
%                 nCoSpikeNoAll = zAll(sum(tractDistance >= 0)+1:end);
                compStringAllSubj = [compStringAllSubj compString(tractDistance >= 0)];

                nCoSpikeCoAllSubj = [nCoSpikeCoAllSubj nCoSpikeCoAll];
                nCoSpikeNoAllSubj = [nCoSpikeNoAllSubj nCoSpikeNoAll];
                distanceAllSubj = [distanceAllSubj tractDistance(tractDistance >= 0)];
                fprintf('%s   angle %i coFire ver %i\n', subject, ang, v)
                try
            %         figure; polarhistogram(phLagsSubj);%waitforbuttonpress; 
            %         savepdf(gcf, fullfile(pdfExport, sprintf('%s_coRipPL.pdf', subject)));
            %         close all;
                end
                
                
        end
        disp('computing CCGs complete.')
        toc
        
        filename = sprintf('coRipple_coFire_ang%i_v%i_fireProp.mat', ang, v);
%         save(fullfile(exportDir, '..', filename), '-v7.3')
        
    end
end

figure;
ph1 = polarhistogram(aSpPhase, -pi:pi/4:pi, 'FaceColor', 'b'); hold on; 
ph2 = polarhistogram(bSpPhase, -pi:pi/4:pi, 'FaceColor', [255/255 197/255 1/255]); hold on;
cm = circ_mean(aSpPhase(~isnan(aSpPhase))');
pp = polarplot([cm cm], [0 max(ph1.BinCounts)]); hold on;
pp.LineWidth = 2;
pp.Color = ph1.FaceColor;
cm = circ_mean(bSpPhase(~isnan(bSpPhase))');
pp = polarplot([cm cm], [0 max(ph2.BinCounts)]); hold on;
pp.Color = ph2.FaceColor;
pp.LineWidth = 2;

[pval, med, P] = circ_cmtest(aSpPhase(~isnan(aSpPhase)), bSpPhase(~isnan(bSpPhase)));

title(sprintf('All  p = %.2f', pval));


figure;
dat = aSpPhase - bSpPhase;
ph1 = polarhistogram(dat, -pi:pi/4:pi, 'FaceColor', [0.8 0.8 0.8]); hold on; 
cm = circ_mean(dat(~isnan(dat))');
pp = polarplot([cm cm], [0 max(ph1.BinCounts)]); hold on;
pp.LineWidth = 2;
pp.Color = ph1.FaceColor;

[pval, z] = circ_rtest(dat(~isnan(dat))');


title(sprintf('All  p = %.2f', pval))

fig = gcf;
fig.Color = 'w';

savepdf(gcf, fullfile(pdfExport, sprintf('PhaseLFP_all.pdf')))
    %%
close all
clr = brewermap(10,'Paired');
smWin = 9; %*fs/1e3;
smoothMethod = 'gaussian';
tag = 'bmovie_perSpike';
% tag = 'Sternberg_perSpike';
% tag = 'Sternberg';
% tag = 'bothStudies';
compStringAll_LFP(cellfun(@(X) isempty(X), compStringAll_LFP))     = [];
compStringAll_phase(cellfun(@(X) isempty(X), compStringAll_phase)) = [];

xtimes = -fs:fs;
xl = [-75 75];
pairLocations = regexp(compStringAllSubj, ' ', 'split')';
pairLocationsLFP = regexp(compStringAll_LFP, ' ', 'split')';
pairLocationsPhase = regexp(compStringAll_phase, ' ', 'split')';
ii = 3;
figure('Position',  [846 419 1158 902]);
for iRa = 1:length(ctxParc) 
    for iRb = iRa:length(ctxParc)
        
        
        if iRa ~= iRb
            regionBundle = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocations);
            regionBundleLFP = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocationsPhase);
        else

            regionBundle = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'contra'), pairLocations);
            regionBundleLFP = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'contra'), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'contra'), pairLocationsPhase);
        end
        crossHem = cellfun(@(X) ~strcmp(X(3), 'ipsi'), pairLocations);
        keep = regionBundle;
        keepLFP = regionBundleLFP;
        keepPhase = regionBundlePhase;
        
        %     keep = true(1, size(ccgAllcoR, 1));  
        % keep = keep & sum(ccgAllnoRAll, 2) > 0;
        % keep = sum(ccgAllcoRtrough, 2) > 0;
        %     keep = sum(ccgAllcoR, 2) > 0;
        xT = xtimes >= xl(1) & xtimes <= xl(2);

        keep = keep;% & sum(ccgAllcoR(1:length(keep),xT), 2) > 0;
        fprintf('%s -- %s %i \n',  ctxParc{iRa}, ctxParc{iRb}, sum(keep))
        % times 

        D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
%         D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';
        D4 = ccgAllnoRAll(keep, :) ./ durAllnoRAll(keep)';

%         D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)';
%         D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)';
%         D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)';
%         D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)';
%         D4 = ccgAllnoRAll(keep, :) ./ nAspikeAllnoRAll(keep)';
        
%         D1 = ccgAllcoRpeak(keep, :);
%         D2 = ccgAllcoRtrough(keep, :);
%         D3 = ccgAllcoR(keep, :);
%         D4 = ccgAllnoRepoch(keep, :);
        
        index = (iRb - 1) * 5 + iRa;
        subplot(5,5, index)
        titles = {'zero lag', 'phase offset', 'combined'};
% 
% 
        eval(sprintf('D = D%i;', ii))
        D = smoothdata(D, 2, smoothMethod, smWin);
        Dnull = smoothdata(D4, 2, smoothMethod, smWin);
%             D = D./Dnull;
%             D(~isfinite(D)) = nan;

        [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        pl1.Color = clr(8,:);
        pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;

        [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(Dnull, 'omitnan'), std(Dnull, 'omitnan')/sqrt(size(Dnull,1))); hold on; 
        pl2.Color = clr(1,:);
        pa2.FaceColor = clr(1,:);
        pl2.LineWidth = 2;
        pa2.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3);

%             xlim([-50 50]*fs/1e3)
%             yl = [min(mean(D(:,xT), 'omitnan'))-0.01 max(mean(D(:,xT), 'omitnan'))+0.01];
        xlim((xl*fs)/1e3)
%             ylim(yl)
        ylabel('firing rate')
        xlabel('time [ms]')
%             % grid on;
        fig = gcf;
        fig.Color = 'w';
% 
% %             title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))

        title(sprintf('%i', sum(ccgAllcoR(keep, xT), 'all')))
            
        
        
%         polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        
        
        if iRb == length(ctxParc); xlabel(ctxParc{iRa}); end
        if iRa == 1; ylabel(ctxParc{iRb}); end

       
    end
end
savepdf(gcf, fullfile(pdfExport, sprintf('coRipCCG_%s_angle_%i_GRID_outside_%s.pdf', titles{ii}, ang, tag)))

figure('Position',  [846 419 1158 902]);
for iRa = 1:length(ctxParc)
    for iRb = iRa %:length(ctxParc)
        
        
        if iRa ~= iRb
            regionBundle = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocations);
            regionBundleLFP = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) any(contains(X(1), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        any(contains(X(2), { ctxParc{iRb}, ctxParc{iRa}})) &  ...
                                        ~strcmp(X{1},X{2}), pairLocationsPhase);
        else
            regionBundle = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocations);
            regionBundleLFP = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocationsPhase);

        end
        crossHem = cellfun(@(X) ~strcmp(X(3), 'ipsi'), pairLocations);
        keep = regionBundle;
        keepLFP = regionBundleLFP;
        keepPhase = regionBundlePhase;

        xT = xtimes >= xl(1) & xtimes <= xl(2);

        keep = keep;% & sum(ccgAllcoR(1:length(keep),xT), 2) > 0;
        fprintf('%s -- %s %i \n',  ctxParc{iRa}, ctxParc{iRb}, sum(keep))

        D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
%         D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';
        D4 = ccgAllnoRAll(keep, :) ./ durAllnoRAll(keep)';

        D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)';
        D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)';
%         D4 = ccgAllnoRAll(keep, :) ./ nAspikeAllnoRAll(keep)';

                
%         D1 = ccgAllcoRpeak(keep, :);
%         D2 = ccgAllcoRtrough(keep, :);
%         D3 = ccgAllcoR(keep, :);
%         D4 = ccgAllnoRepoch(keep, :);
        
        index = (iRb - 1) * 5 + iRa;
        subplot(5,5, index)
        titles = {'zero lag', 'phase offset', 'combined'};

        eval(sprintf('D = D%i;', ii))
        D = smoothdata(D, 2, smoothMethod, smWin);
        Dnull = smoothdata(D4, 2, smoothMethod, smWin);
%             D = D./Dnull;
%             D(~isfinite(D)) = nan;

        [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        pl1.Color = clr(8,:);
        pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;

        [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(Dnull, 'omitnan'), std(Dnull, 'omitnan')/sqrt(size(Dnull,1))); hold on; 
        pl2.Color = clr(1,:);
        pa2.FaceColor = clr(1,:);
        pl2.LineWidth = 2;
        pa2.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3);

%             xlim([-50 50]*fs/1e3)
%             yl = [min(mean(D(:,xT), 'omitnan'))-0.01 max(mean(D(:,xT), 'omitnan'))+0.01];
        xlim((xl*fs)/1e3)
%             ylim(yl)
        ylabel('firing rate')
        xlabel('time [ms]')
%             % grid on;
        fig = gcf;
        fig.Color = 'w';
% 
% %             title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))

        title(sprintf('%i', sum(ccgAllcoR(keep, xT), 'all')))
            
        
                
        
        if iRb == length(ctxParc); xlabel(ctxParc{iRa}); end
        if iRa == 1; ylabel(ctxParc{iRb}); end

       
    end
end
savepdf(gcf, fullfile(pdfExport, sprintf('coRipCCG_%s_angle_%i_GRID_within_%s.pdf', titles{ii}, ang, tag)))

figure('Position',  [1739 625 822 717]);
figure('Position', [560 188 1028 690]);

for iRa = 1 %:length(ctxParc)
    for iRb = iRa %:length(ctxParc)
        sameBundle = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocations);
        sameBundleLFP = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsLFP);
        sameBundlePhase = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsPhase);
        
        
       
        crossHem = cellfun(@(X) ~strcmp(X(3), 'ipsi'), pairLocations);

        keep = ~sameBundle;
        keepLFP = ~sameBundleLFP;
        keepPhase = ~sameBundlePhase;

        xT = xtimes >= xl(1) & xtimes <= xl(2);

        keep = keep & sum(ccgAllcoR(1:length(keep),xT), 2) > 0;
%         keep = keep & nAspikeAllcoR(1:length(keep))' > 1;
        fprintf('%s -- %s %i \n',  ctxParc{iRa}, ctxParc{iRb}, sum(keep))

        D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
%         D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';
        D4 = ccgAllnoRAll(keep, :) ./ durAllnoRAll(keep)';
        figure(3)

        
        D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)' * 1e3;
        D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)' * 1e3;
        D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)' * 1e3;
        D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)' * 1e3;
        D4 = ccgAllnoRAll(keep, :) ./ nAspikeAllnoRAll(keep)' * 1e3;

%         D1 = ccgAllcoRpeak(keep, :);
%         D2 = ccgAllcoRtrough(keep, :);
%         D3 = ccgAllcoR(keep, :);
%         D4 = ccgAllnoRepoch(keep, :);
        
        index = (iRb - 1) * 5 + iRa;
%         subplot(5,5, index)
        titles = {'zero lag', 'phase offset', 'combined'};

% 
% 
        eval(sprintf('D = D%i;', ii))
        D = smoothdata(D, 2, smoothMethod, smWin);
        Dnull = smoothdata(D4, 2, smoothMethod, smWin);
%             D = D./Dnull;
%             D(~isfinite(D)) = nan;
        subplot(2,2,1)
        [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        pl1.Color = clr(8,:);
        pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;

        [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(Dnull, 'omitnan'), std(Dnull, 'omitnan')/sqrt(size(Dnull,1))); hold on; 
        pl2.Color = clr(1,:);
        pa2.FaceColor = clr(1,:);
        pl2.LineWidth = 2;
        pa2.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3);
        xlim((xl*fs)/1e3)
        ylabel('across bundle')
        title(sprintf('total spikes: %i', sum(ccgAllcoR(keep, xT), 'all')))


        subplot(2,2,2)
        plot(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan')./mean(Dnull, 'omitnan'))

%             xlim([-50 50]*fs/1e3)
%             yl = [min(mean(D(:,xT), 'omitnan'))-0.01 max(mean(D(:,xT), 'omitnan'))+0.01];
        vline([-99:11:99]*fs/1e3);
        xlim((xl*fs)/1e3)
%             ylim(yl)
        ylabel('firing rate')
        xlabel('time [ms]')
%             % grid on;
        fig = gcf;
        fig.Color = 'w';
        ylabel('coR/noR')

% 
% %             title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))

            
            
        
        
%         polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        
        
        
         figure(4);
        %     plot(-100:100, mean(lfpA([1:1.5e4],:)), 'r'); hold on;
        %     plot(-100:100, mean(lfpB([1:1.5e4],:)), 'b'); hold on;
        %     plot(-100:100, mean(lfpA(2.5e4:end,:)), 'r'); hold on;
        %     plot(-100:100, mean(lfpB(2.5e4:end,:)), 'b'); hold on;
        subplot(2,2,1)
        lfpApl = lfpA(keepLFP,:) * 1e6; rbApl = rbA(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-100:100, mean(lfpApl, 'omitnan'),std(lfpApl, 'omitnan')/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-5 10])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location A')


        subplot(2,2,2)
        [pl1, pa1] = boundedline(-100:100, mean(rbApl, 'omitnan'),std(rbApl, 'omitnan')/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')

        subplot(2,2,3)
        lfpBpl = lfpB(keepLFP,:) * 1e6; rbBpl = rbB(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-100:100, mean(lfpBpl, 'omitnan'),std(lfpBpl, 'omitnan')/sqrt(size(lfpBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-3 20])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location B')
        subplot(2,2,4)
        [pl1, pa1] = boundedline(-100:100, mean(rbBpl, 'omitnan'),std(rbBpl, 'omitnan')/sqrt(size(rbBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')
        
        figure; 
        polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        title('phase lags during co-spike co-ripples')

       
    end
end
figure(3)
savepdf(gcf, fullfile(pdfExport, sprintf('coRipCCG_%s_angle_%i_all_outside_%s.pdf', titles{ii}, ang, tag)))


% figure('Position',  [2010 898 404 440]);
figure('Position', [560 188 1028 690]);

for iRa = 1 %:length(ctxParc)
    for iRb = iRa %:length(ctxParc)
        sameBundle = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocations);
        sameBundleLFP = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsLFP);
        sameBundlePhase = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsPhase);
        
        
       
        crossHem = cellfun(@(X) ~strcmp(X(3), 'ipsi'), pairLocations);

        keep = sameBundle;
        keepLFP = sameBundleLFP;
        keepPhase = sameBundlePhase;

        xT = xtimes >= xl(1) & xtimes <= xl(2);

        keep = keep & sum(ccgAllcoR(1:length(keep),xT), 2) > 0;
%         keep = keep & nAspikeAllcoR(1:length(keep))' > 1;

        fprintf('%s -- %s %i \n',  ctxParc{iRa}, ctxParc{iRb}, sum(keep))
        
        figure(3)

        D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
%         D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';
        D4 = ccgAllnoRAll(keep, :) ./ durAllnoRAll(keep)';

        D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)' * 1e3;
        D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)' * 1e3;
        D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)' * 1e3;
        D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)' * 1e3;
        D4 = ccgAllnoRAll(keep, :) ./ nAspikeAllnoRAll(keep)' * 1e3;

%         D1 = ccgAllcoRpeak(keep, :);
%         D2 = ccgAllcoRtrough(keep, :);
%         D3 = ccgAllcoR(keep, :);
%         D4 = ccgAllnoRepoch(keep, :);
        index = (iRb - 1) * 5 + iRa;
%         subplot(5,5, index)
        titles = {'zero lag', 'phase offset', 'combined'};

% 
% 
        eval(sprintf('D = D%i;', ii))
        D = smoothdata(D, 2, smoothMethod, smWin);
        Dnull = smoothdata(D4, 2, smoothMethod, smWin);
%             D = D./Dnull;
%             D(isinf(D)) = nan;
        subplot(2,2,3)
        [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        pl1.Color = clr(8,:);
        pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;




        [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(Dnull, 'omitnan'), std(Dnull, 'omitnan')/sqrt(size(Dnull,1))); hold on; 
        pl2.Color = clr(1,:);
        pa2.FaceColor = clr(1,:);
        pl2.LineWidth = 2;
        pa2.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3);
        xlim((xl*fs)/1e3)  
        ylabel('within bundle')
        title(sprintf('total spikes: %i', sum(ccgAllcoR(keep, xT), 'all')))



        subplot(2,2,4)
        plot(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan')./mean(Dnull, 'omitnan'))
        vline([-99:11:99]*fs/1e3);
        xlim((xl*fs)/1e3)%             ylim(yl)
        ylabel('firing rate')
        xlabel('time [ms]')
%             % grid on;
        fig = gcf;
        fig.Color = 'w';
        ylabel('noR/coR')
% 
% %             title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))

          
        
        
%         polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        
        
%         ylabel('within bundle')
        
         figure(6);
        subplot(2,2,1)
        lfpApl = lfpA(keepLFP,:) * 1e6; rbApl = rbA(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-100:100, mean(lfpApl, 'omitnan'),std(lfpApl, 'omitnan')/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-5 10])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location A')


        subplot(2,2,2)
        [pl1, pa1] = boundedline(-100:100, mean(rbApl, 'omitnan'),std(rbApl, 'omitnan')/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')

        subplot(2,2,3)
        lfpBpl = lfpB(keepLFP,:) * 1e6; rbBpl = rbB(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-100:100, mean(lfpBpl, 'omitnan'),std(lfpBpl, 'omitnan')/sqrt(size(lfpBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-3 20])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location B')
        subplot(2,2,4)
        [pl1, pa1] = boundedline(-100:100, mean(rbBpl, 'omitnan'),std(rbBpl, 'omitnan')/sqrt(size(rbBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')
        
        figure; 
        polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        title('phase lags during co-spike co-ripples')



       
    end
end
figure(3)


savepdf(gcf, fullfile(pdfExport, sprintf('coRipCCG_%s_angle_%i_all_within_%s_Aspike.pdf', titles{ii}, ang, tag)))
%%
for iRa = 1 %:length(ctxParc)
    for iRb = iRa %:length(ctxParc)
        sameBundle = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocations);
        sameBundleLFP = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsLFP);
        sameBundlePhase = cellfun(@(X) strcmp(X(1), X(2)) & strcmp(X(3), 'ipsi'), pairLocationsPhase);
        
        if iRa ~= iRb
            regionBundle = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}), pairLocations);
            regionBundleLFP = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}), pairLocationsPhase);
        else
            regionBundle = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocations);
            regionBundleLFP = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocationsLFP);
            regionBundlePhase = cellfun(@(X) strcmp(X(1), ctxParc{iRa}) &  strcmp(X(2), ctxParc{iRb}) & strcmp(X(3), 'ipsi'), pairLocationsPhase);
        end
        crossHem = cellfun(@(X) ~strcmp(X(3), 'ipsi'), pairLocations);
%         keep = regionBundle;
%         keepLFP = regionBundleLFP;
%         keepPhase = regionBundlePhase;
        keep = sameBundle;
        keepLFP = sameBundleLFP;
        keepPhase = sameBundlePhase;
        %     keep = true(1, size(ccgAllcoR, 1));  
        % keep = keep & sum(ccgAllnoRAll, 2) > 0;
        % keep = sum(ccgAllcoRtrough, 2) > 0;
        %     keep = sum(ccgAllcoR, 2) > 0;
        keep = keep; % & sum(ccgAllcoR(1:length(keep),:), 2) < 5;
        fprintf('%s -- %s %i \n',  ctxParc{iRa}, ctxParc{iRb}, sum(keep))
        % times 
%         D1 = ccgAllcoRpeak(keep, :) ./ sum(ccgAllcoRpeak(keep, :), 2);
%         D2 = ccgAllcoRtrough(keep, :) ./ sum(ccgAllcoRtrough(keep, :), 2);
%         D3 = ccgAllcoR(keep, :) ./ sum(ccgAllcoR(keep, :), 2);
        % D4 = ccgAllnoRAll(keep, :) ./ sum(ccgAllnoRAll(keep, :), 2);
        D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
        D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
        D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
        D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';

        titles = {'zero lag', 'phase offset', 'combined'};
        for ii = 1:3
            figure('Position',  [846 779 219 542]);
            subplot(2,1,1)
            eval(sprintf('D = D%i;', ii))
            D = smoothdata(D, 2, smoothMethod,smWin);
            [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
            pl1.Color = clr(8,:);
            pa1.FaceColor = clr(8,:);
            pl1.LineWidth = 2;
            pa1.FaceAlpha = 0.3;

        %     D = smoothdata(D4(:,800:1200), 2, smoothMethod,smWin);
            D = smoothdata(D4, 2, smoothMethod,smWin);
            [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
            % [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(ccgAllcoRtrough, 'omitnan'), std(ccgAllcoRtrough, 'omitnan')/sqrt(size(ccgAllcoRtrough,1))); hold on; 
            pl2.Color = clr(1,:);
            pa2.FaceColor = clr(1,:);
            pl2.LineWidth = 2;
            pa2.FaceAlpha = 0.3;
            vline([-99:11:99]*fs/1e3);
            xlim([-300 300]*fs/1e3)

            ylabel('firing rate')
            xlabel('time [ms]')
            legend([pl1, pl2], 'co-ripple', 'no-ripple', 'location', 'best')


            subplot(2,1,2)
            eval(sprintf('D = D%i;', ii))
            D = smoothdata(D, 2, smoothMethod, smWin);
            [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
            pl1.Color = clr(8,:);
            pa1.FaceColor = clr(8,:);
            pl1.LineWidth = 2;
            pa1.FaceAlpha = 0.3;

        %     D = smoothdata(D4(:,800:1200), 2, smoothMethod, smWin);
            D = smoothdata(D4, 2, smoothMethod, smWin);
            [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
            pl2.Color = clr(1,:);
            pa2.FaceColor = clr(1,:);
            pl2.LineWidth = 2;
            pa2.FaceAlpha = 0.3;
            vline([-99:11:99]*fs/1e3);

            xlim([-50 50]*fs/1e3)
            ylabel('firing rate')
            xlabel('time [ms]')
            % grid on;
            fig = gcf;
            fig.Color = 'w';

            title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))
            savepdf(gcf, fullfile(pdfExport, sprintf('coRipCCG_%s_angle_%i_%s-%s_%s.pdf', titles{ii}, ang,  ctxParc{iRa}, ctxParc{iRb}, tag)))
        end
        
        figure; polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);


        %     figure('Position',  [846.0000 948 754.0000 373]); 
        %     
        %     subplot(2,3,1)
        %     Dpeak = mean(D1, 'omitnan')./mean(D4, 'omitnan');
        %     Dpeak(~isfinite(Dpeak)) = nan;
        %     Dpeak = smoothdata(Dpeak, 2, 'gaussian', smWin);
        %     % D = smoothdata(ccgAllcoR, 2, 'gaussian', smWin);
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     pl1 = plot(-coRipPad:unitCCGstep:coRipPad, Dpeak); hold on; 
        %     clrP = clr(8,:)*1.5; clrP(clrP >1) = 1;
        %     pl1.Color = clrP;
        %     pl1.LineWidth = 2;
        %     % pa1.FaceColor = clr(8,:);
        %     % pa1.FaceAlpha = 0.3;
        %     vline([-99:11:99])
        %     % ylim([2.5 4.5])
        %     xlim([-100 100])
        %     
        %     subplot(2,3,2)
        %     Dall = mean(D3, 'omitnan')./mean(D4, 'omitnan');
        %     % Dall = D3./D4; %mean(D3, 'omitnan')./mean(D4, 'omitnan');
        %     Dall(~isfinite(Dall)) = nan;
        %     Dall = smoothdata(Dall, 2, 'gaussian', smWin);
        %     % D = smoothdata(ccgAllcoR, 2, 'gaussian', smWin);
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     % pl2 = plot(-coRipPad:unitCCGstep:coRipPad, mean(Dall, 'omitnan')); hold on; 
        %     pl2 = plot(-coRipPad:unitCCGstep:coRipPad, Dall); hold on; 
        %     pl2.Color = clr(8,:);
        %     pl2.LineWidth = 2;
        %     % pa1.FaceColor = clr(8,:);
        %     % pa1.FaceAlpha = 0.3;
        %     vline([-99:11:99])
        %     % ylim([2.5 4.5])
        %     xlim([-100 100])
        %     
        %     subplot(2,3,3)
        %     Dtr = mean(D2, 'omitnan')./mean(D4, 'omitnan');
        %     Dtr(~isfinite(Dtr)) = nan;
        %     Dtr = smoothdata(Dtr, 2, 'gaussian', smWin);
        %     % D = smoothdata(ccgAllcoR, 2, 'gaussian', smWin);
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     pl3 = plot(-coRipPad:unitCCGstep:coRipPad, Dtr); hold on; 
        %     pl3.Color = clr(8,:)*0.7;
        %     pl3.LineWidth = 2;
        %     vline([-99:11:99])
        %     % ylim([2.5 4.5])
        %     xlim([-100 100])
        %     ylabel('coR / noR')
        %     xlabel('time [ms]')
        %     
        %     % legend([pl1, pl3, pl2], 'zero-lag', 'pi-lag', 'all', 'location', 'best')
        %     
        %     subplot(2,3,4)
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     pl1 = plot(-coRipPad:unitCCGstep:coRipPad, Dpeak); hold on; 
        %     pl1.Color = clrP;
        %     pl1.LineWidth = 2;
        %     vline([-99:11:99])
        %     xlim([-50 50] )
        %     % pa1.FaceColor = clr(8,:);
        %     % pa1.FaceAlpha = 0.3;
        %     % ylim([3 4.5])
        %     % Dall = smoothdata(Dall, 2, 'gaussian', smWin);
        %     % D = smoothdata(ccgAllcoR, 2, 'gaussian', smWin);
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     % pl2 = plot(-coRipPad:unitCCGstep:coRipPad, mean(Dall, 'omitnan')); hold on; 
        %     subplot(2,3,5)
        %     
        %     pl2 = plot(-coRipPad:unitCCGstep:coRipPad, Dall); hold on; 
        %     pl2.Color = clr(8,:);
        %     pl2.LineWidth = 2;
        %     vline([-99:11:99])
        %     xlim([-100 100] )
        %     
        %     subplot(2,3,6)
        %     
        %     % Dtr = smoothdata(Dtr, 2, 'gaussian', smWin);
        %     % D = smoothdata(ccgAllcoR, 2, 'gaussian', smWin);
        %     % [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
        %     pl1 = plot(-coRipPad:unitCCGstep:coRipPad, Dtr); hold on; 
        %     pl1.Color = clr(8,:)*0.7;
        %     pl1.LineWidth = 2;
        %     
        %     vline([-99:11:99])
        %     xlim([-100 100] )
        %     ylabel('coR / noR')
        %     xlabel('time [ms]')
        %     fig = gcf;
        %     fig.Color = 'w';
        %     savepdf(gcf, fullfile(pdfExport,  sprintf('coRipCCG_ratio_angle_%i.pdf',  ang)))


        figure; polarhistogram(phLagsAll(keepPhase), -pi:pi/20:pi);
        title(sprintf('%s <-> %s', ctxParc{iRa}, ctxParc{iRb}))

        savepdf(gcf, fullfile(pdfExport,  sprintf('phaselags_angle_%i_%s-%s_%s.pdf',  ang, ctxParc{iRa}, ctxParc{iRb}, tag)))


        figure('Position', [560 188 1028 690]);
        %     plot(-100:100, mean(lfpA([1:1.5e4],:)), 'r'); hold on;
        %     plot(-100:100, mean(lfpB([1:1.5e4],:)), 'b'); hold on;
        %     plot(-100:100, mean(lfpA(2.5e4:end,:)), 'r'); hold on;
        %     plot(-100:100, mean(lfpB(2.5e4:end,:)), 'b'); hold on;
        subplot(2,2,1)
        lfpApl = lfpA(keepLFP,:) * 1e6; rbApl = rbA(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-200:200, mean(lfpApl),std(lfpApl)/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-5 10])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location A')


        subplot(2,2,2)
        [pl1, pa1] = boundedline(-200:200, mean(rbApl),std(rbApl)/sqrt(sum(keepLFP)), 'r'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')

        subplot(2,2,3)
        lfpBpl = lfpB(keepLFP,:) * 1e6; rbBpl = rbB(keepLFP,:) * 1e6;
        [pl1, pa1] = boundedline(-200:200, mean(lfpBpl),std(lfpBpl)/sqrt(size(lfpBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        %     ylim([-3 20])
        ylabel('LFP (microV)')
        xlabel('time from spike in A')
        legend(pl1, 'location B')
        subplot(2,2,4)
        [pl1, pa1] = boundedline(-200:200, mean(rbBpl),std(rbBpl)/sqrt(size(rbBpl,1)), 'b'); hold on;
        %     pa1.FaceColor = clr(8,:);
        pl1.LineWidth = 2;
        pa1.FaceAlpha = 0.3;    
        vline([-99:11:99]*fs/1e3)
        xlim([-200 200]*fs/1e3)
        ylabel('rippleband (microV)')
        xlabel('time from spike in A')
        sgtitle(sprintf('%s <-> %s', ctxParc{iRa}, ctxParc{iRb}))

        savepdf(gcf, fullfile(pdfExport,  sprintf('coRip_LFP_angle_%i_%s-%s_%s.pdf',  ang, ctxParc{iRa}, ctxParc{iRb}, tag)))
        
%         close all
    end
end

    


%%
ii = 10;
close all;
% 
% figure; 
% subplot(4,1,1); plot(D3(ii,:)); hold on;  
% subplot(4,1,2); plot(D1(ii,:)); hold on;  
% subplot(4,1,3); plot(D2(ii,:));
% subplot(4,1,4); plot(D4(ii,:));
% 
% figure;
% iiM = max(ccgAllcoRpeak,[], 2);
plD = ccgAllcoRpeak(keep,:);

figure; 
plot(-200:200, smoothdata(sum(plD(randperm(size(plD,1), 1000),:)), 'gaussian', 7));
vline([-99:11:99])
xlim([-100 100])
% plD = ccgAllcoRpeak(keep,:);
% for ii = 1:sum(keep)
%     plot(-200: 200, plD(ii,:)); hold on;
%     vline([-99:11:99]);
%     xlim([-100 100])
%     waitforbuttonpress; clf; 
% end


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
%     ylim([-0.6 1.2])

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

% ylim([-0.4 0.8])
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
% ylim([-0.4 0.8])

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
% ylim([-0.4 0.8])

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










