function Movie_corip_analysis(subjIdcs,runBoot,zeroPh,ctOpt,HGburst,stimRespIdx,bilat, taskContrast)

if nargin<2
    runBoot=1
end
if nargin<3
    zeroPh=0
end
if nargin<4
    ctOpt=3;
end
if nargin<5
    HGburst=0
end
if nargin<6
    stimRespIdx=[1];
end
if nargin<7
    bilat=0;
end
if nargin < 8
    taskContrast = 'recognition';
%     taskContrast = 'performance';
end
if HGburst
    HGapnd='_HGburst';
else
    HGapnd='';
end

if zeroPh
    zeroPhApnd='_zeroPh';
else
    zeroPhApnd='';
end

if ctOpt==1
    ctApnd='' %count co-ripple centers
elseif ctOpt==2
    ctApnd='_onset' %count co-ripple starts
elseif ctOpt==3
    ctApnd='_dens' %count co-ripple durations (density)
end

if bilat
    rhApnd='_bilat'
else
    rhApnd='';
end

subjects =  {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};
subjIdcs = 1:numel(subjects);
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
dataDirectory =     '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/preprocess';

LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '*task*'));
taskfiles = {taskfiles.name}';

outDir='/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/code/ConnectionPlots/out';
if ~exist(outDir,'dir');unix(sprintf('mkdir -p %s',outDir));end

close all
for s = 1:numel(subjects)
    subj = subjects{s};
    clc
    fprintf('====================================================\n')
    fprintf('          SUBJECT %s                          \n',subj)
    
    
    tic;fprintf('loading data...\n')
    load('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/MG49/data_1kHz/MG49_1kHz_unitsRemoved_WN.mat', 'data');
    Fs = 1000;
    
    modifier = '';
    tag = ['wake_NC_',modifier];
    filename = sprintf('%s_LFP_macro_ripple_stats_%s.mat', subj, tag);
    load(fullfile(matExportFolder, filename))
    
    modifier = '';
    tag = ['wake_NC_',modifier];
    filename = sprintf('%s_LFP_micro_ripple_stats_%s.mat', subj, tag);
    microObj = load(fullfile(matExportFolder, filename));

    f = contains(LFPfilesMacr, subj);
    macrObj = load(fullfile(dataDirectory, LFPfilesMacr{f}), 'times');
    LFPtimeMacr = macrObj.times;

    f = contains(LFPfilesMicr, subj);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}), 'times');
    LFPtimeMicr = micrObj.times;

    rippleStats.recordingType = [repmat({'macro'}, [1 length(rippleStats.locs)]) repmat({'micro'}, [1 length(microObj.rippleStats.locs)])];
    rippleStats.locs = [cellfun(@(X) round(LFPtimeMacr(X) * Fs), rippleStats.locs, 'UniformOutput', false), ...
                        cellfun(@(X) round(LFPtimeMicr(X) * Fs), microObj.rippleStats.locs, 'UniformOutput', false)];
    rippleStats.window(cellfun(@(X) isempty(X), rippleStats.window)) = {[1 1]};
    microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
    rippleStats.window = [cellfun(@(X) round([LFPtimeMacr(X(:,1)); LFPtimeMacr(X(:,2))] * 1e3)', rippleStats.window, 'UniformOutput', false), ...
                          cellfun(@(X) round([LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))] * 1e3)', microObj.rippleStats.window, 'UniformOutput', false)];

    rippleStats.chanLabels = [rippleStats.chanLabels; microObj.rippleStats.chanLabels];
    rippleStats.macroTimes = LFPtimeMacr;
    rippleStats.microTimes = LFPtimeMicr;
    
    rippleStats.recordingLength = max([rippleStats.recordingLength, microObj.rippleStats.recordingLength]);



%     [labParc,~,~,~,~,rejChan,coords]=retrieveChanParcels_MG(subj);
%     toc
%     
%     tic;fprintf('channel bookkeeping...')
%     
%     rmCh =[rejChan; true(size(labParc,1)-numel(rejChan),1)]|strcmp(labParc(:,1),'N/A');
%     if ~bilat
%         rmCh=rmCh|strcmp(labParc(:,3),'rh');
%     end
%     labParc(rmCh,:)=[];
%     data(:,rmCh)=[];
%     nCh=sum(~rmCh);
%     parcMap=readcell('/home/jgarret/TaskAnalysis/destrieuxParcMap_split.xls');
%     [reg,~,iReg]=unique(parcMap(:,1));
    
    
    toc
    f = contains(taskfiles, subj);
    trials = load(fullfile(dataDirectory, taskfiles{f}));

    %% analysis
    respApnd={'','_resp'};
    % allTrl=loadTrialInfo_v2(subj,0);
    trlOut = trials.start_time(2:end); resp = round(trials.response_time(2:end) * 1e3);
    switch taskContrast
        case 'performance'
            codes = trials.response_correct(2:end);
        case 'recognition'
            codes = contains(trials.stimulus_file(2:end), 'new');
    end
    trl = [];
    trl(:,1) = round(trlOut * Fs);
    trl(:,2) = codes;
    trl(:,3) = trials.response_confidence(2:end);
    
    for r=stimRespIdx
        tic;fprintf('trial bookkeeping...')
        
        if r==1
            fprintf('stimulus locked...')
            novel = {0, 1, [0 1], [0 1], [0 1]};
            confidence = {[1:3], [1:3], 3, [1:2], [1:3]};
            switch taskContrast
                case 'performance'
                    condNames={'correct','incorrect','base200','base100'};

                case 'recognition'
                    condNames={'familiar','novel','confident','unsure','all','base200','base100'};
            end
            stimTimes=cellfun(@(x,y) trl(any(trl(:,2)==x,2) & any(trl(:,3)==y,2),1), novel, confidence, 'uni',0)';
            allStimTimes=trl(any(trl(:,2)==[novel{[1 2]}],2),1);
        elseif r==2
            fprintf('response locked...')
            codes={1};
            condNames={'TP'};
            stimTimes=cellfun(@(x) resp(any(resp(:,2)==x,2),1),codes,'uni',0);
            allStimTimes=resp(any(resp(:,2)==[codes{[2 4]}],2),1);
        end
        
        
        
        
        stimTimes{end+1}=allStimTimes;
        stimTimes{end+1}=allStimTimes;
        rip_bounds_mask = false(length(rippleStats.locs), rippleStats.recordingLength);
        gCh=1:length(rippleStats.locs); %find(~rmCh)';
        nCh = length(gCh);
        tic;fprintf('populating ripples...\n')
        for ch = gCh
            %     rip_mask(ch,rippleStats.locs{ch}) = true;
            for rip=1:size(rippleStats.window{ch},1)
                if  rippleStats.window{ch}(rip,1) > 0
                    rip_bounds_mask(ch==gCh, rippleStats.window{ch}(rip,1):rippleStats.window{ch}(rip,2)) = true;
                end
                %         rip_peak_mask(rippleStats.locs{ch}(rip),ch==gCh) = true;
            end
        end
%         rip_bounds_mask = rip_bounds_mask(:,mask);
%         corip_bounds_mask = rip_bounds_mask & permute(rip_bounds_mask,[1 3 2]);
        nS = length(rip_bounds_mask);
        coripBnds=cell(nCh,nCh);
        corip_cent_mask = false(nS,nCh,nCh);
        toc
        
        if zeroPh
            [b,a]=butter(3,[70 110]/(Fs/2),'bandpass');
            RBang=single(angle(hilbert(filtfilt(b,a,double(data)))));
            coripAng = RBang - permute(RBang,[1 3 2]);
        end
        
        tic;fprintf('finding co-ripples\n')
        for c1=1:nCh
            for c2=1:c1
                corip_bounds_mask = rip_bounds_mask(c1,:) & rip_bounds_mask(c2,:);
                coripBnds{c1,c2}=mask2bounds(corip_bounds_mask);
                coripBnds{c1,c2}(diff(coripBnds{c1,c2},1,2)<25,:)=[];
                if zeroPh
                    for rip=size(coripBnds{c1,c2}):-1:1
                        crTmp=coripBnds{c1,c2}(rip,:);
                        tmpAng = double(coripAng(crTmp(1):crTmp(2),c1,c2));
                        eveAng = angle(mean(exp(1i*tmpAng)));
                        %                 if eveAng>pi/4 || eveAng<-pi/4
                        if eveAng>pi/3 || eveAng<-pi/3
                            coripBnds{c1,c2}(rip,:)=[];
                        end
                    end
                end
                if ctOpt==1
                    corip_cent_mask(:,c1,c2) = bounds2mask(floor(mean(coripBnds{c1,c2},2)),nS);
                elseif ctOpt==2
                    corip_cent_mask(:,c1,c2) = bounds2mask(coripBnds{c1,c2}(:,1),nS);
                elseif ctOpt==3
                    corip_cent_mask(:,c1,c2) = bounds2mask(coripBnds{c1,c2},nS);
                end
            end
            if mod(c1, 10) == 1
%                 clc
%                 fprintf('computing coRipple Rates for session %s\n', subj)
%                 fprintf('channel %i / %i\n', c1, nCh)
                
                prog = ceil(c1/nCh*10);
                Bar = repmat('*', [1 prog]);
                Remain = repmat('-', [1 10-prog]);
                progBar = [Bar, Remain];
                fprintf('%s\n', progBar)
            end
        end
        toc
        
        tic;fprintf('allocating (co-)ripples into trials x bins...')
        coripCond=cell(numel(stimTimes),1);
        if r==1
            binVec=-1000:100:3000;
        elseif r==2
            binVec=-1000:100:1000;
        end
        nBin=numel(binVec)-1;
        sampVec=cell2mat(arrayfun(@(x) binVec(x)+1:binVec(x+1),(1:nBin)','uni',0));
        corip_cent_tmp=permute(corip_cent_mask,[1 4 5 2 3]);
        for t=1:numel(stimTimes)
            fprintf('%d/%d...',t,numel(stimTimes))
            condStim=permute(stimTimes{t},[2 3 1]);
            nStimCond=numel(condStim);
            if nStimCond < 2
                continue
            end
            
            if t==length(stimTimes)-1
%                 condSamp = repmat([[-599:-500];[-499:-400]], [15 1]) + condStim;
                condSamp = repmat([[-199:-100];[-99:0]], [20 1]) + condStim;
            elseif t==length(stimTimes)
%                 condSamp = repmat(-499:-400, [30 1]) + condStim;
                condSamp = repmat(-99:0, [40 1]) + condStim;
            else
                condSamp=sampVec+condStim;
            end
            
            coripCond{t}=false(nCh,nCh,nBin,nStimCond);
            coripCondTmp=reshape(corip_cent_tmp(condSamp(:),:,:,:,:),[size(condSamp) nCh nCh]);
            coripCond{t}=permute(sum(coripCondTmp,2),[4 5 1 3 2]);
            if any(sum(coripCond{t},4)>2^16-1,'all')
                error('co-ripple metric exceeds uint16 val')
            end
        end
        toc
        
        tic;fprintf('saving data...\n')
        save(sprintf('%s/coripChanDat%s_%s%s%s%s%s.mat',outDir,respApnd{r},subj,rhApnd,HGapnd,zeroPhApnd,ctApnd),...
            'coripCond','gCh','condNames','subj','condNames','-v7.3','-nocompression')
        toc
        
    end
end

clear coripCond corip_cent_tmp coripCondTmp corip_cent_mask corip_bounds_mask rip_bounds_mask data
if runBoot
    Movie_genCoripDist(subjIdcs,1,any(stimRespIdx==1),any(stimRespIdx==2),zeroPh,ctOpt,HGburst,bilat)
end

end