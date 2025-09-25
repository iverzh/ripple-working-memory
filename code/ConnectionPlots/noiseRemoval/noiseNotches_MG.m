function noiseNotches_MG(subjIdcs,applyFilters)
% apply notches to noise bands in power spectrum for FW data

%subjIdcs=[68 86 130 133 190 233];

%subjIdcs= [311 297 233 226 190 133 130 111 101 -89 86 68];
%subjIdcs=[311 297 -240 226 111 101 -89];
% subjIdcs=[442 321 313 259 255 170 90 87 71 70 67];
    nS=numel(subjIdcs);
    
    if nargin<2
        applyFilters=0;
    end
%nHrz=4;
nHrz=ceil(nS/(4/(1+double(applyFilters))));
nVrt=ceil(nS*(1+double(applyFilters))/nHrz);
gap=[0.04 0.02];
marg_h=[0.05 0.03];
marg_w=[0.04 0.02];

fRoot='/space/seh8/5/halgdev/projects/jgarret/altnalysis_native_fs/MG/';
% fRoot='/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/';
outDir='/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/';
rsFs=1000;

subjects=genSubjList(subjIdcs,'MG');

notchParams=readcell('/home/jgarret/TaskAnalysis/lang/WN/noiseRemoval/subjNotches_WN.xls');
subjCols=find(cellfun(@(x) ~ismissing(x), notchParams(1,2:end)))+1;

h=figure('Position',get(0,'Screensize'));
for s=1:nS
    subj=subjects{s};
    nSpPre = floor((s-1)/nHrz)*(1+double(applyFilters))*nHrz + mod(s-1,nHrz)+1;
    nSpPost = nSpPre + nHrz;
    inFil=sprintf('%s/%s/tasks/%s_WN_noNotch.mat',fRoot,subj,subj);
    outFil=sprintf('%s/%s/tasks/%s_WN.mat',fRoot,subj,subj);
    outFilRs=sprintf('%s/%s/tasks/%s_WN.mat',outDir,subj,subj);
    
    if ~exist(inFil,'file')
        warning('no file labeled noNotch found, copying base file')
        unix(sprintf('cp %s %s', outFil,inFil));
    end
    if ~exist(outFilRs,'file')
        warning('no resampled analysis file detected, copying base file before appending clean data')
        unix(sprintf('cp %s %s', outFilRs,outFil));
    end
    
    load(inFil,'data','Fs','label','bpOp','bpLabel','trigData')
    if subjIdcs(s)==49
        largeArt=mask2bounds(any(abs(data)>1000,2));
        fprintf('%d large artifacts detected, removing total of %.2f seconds\n',size(largeArt,1),sum(diff(largeArt,1,2)+1250)/1000)
        data(bounds2mask(largeArt,size(data,1),250,1000),:)=0;
    elseif subjIdcs(s)==-51
        data=data-median(data,1);
        [trl,~,rej]=loadTrials(subj,'WN');
%         trl(rej,:)=[];
        trlMsk=bounds2mask(trl(:,1)+[-500 10000],size(data,1));
%         data(~trlMsk,:)=0;
%         data(trlMsk,:)=[];
    else
        data=data-median(data,1);
        [trialDat,~,rej]=loadTrials(subj,'WN');
        trialDat(rej,:)=[];
%         rmDat = ~bounds2mask(trialDat(:,1)+[0 2000],length(data));
%         keepDat = bounds2mask(trialDat(:,1)/2+[-1000 -500],length(data));
%         keepDat = bounds2mask(trialDat(:,1)/2+[-400 0],length(data));
%         data(rmDat,:)=0;
%         data=data(keepDat,:);
%         data(~k
%             data=data(round(size(data,1)/2)+[0:5000],:);
%         largeArt=mask2bounds(mean(abs(data),2)>1000);
%         fprintf('%d large artifacts detected, removing total of %.2f seconds\n',size(largeArt,1),sum(diff(largeArt,1,2)+1250)/1000)
%         data(bounds2mask(largeArt,size(data,1),250,1000),:)=0;
    end
    
    
    Fs=round(Fs);
    subtightplot(nVrt,nHrz,nSpPre,gap,marg_h,marg_w)
    pwelch(data(:,1:5:end),[],[],[],Fs);
    title([subj ' pre'])
    if mod(nSpPre,nHrz)~=1
        ylabel('')
    end
    if applyFilters || nSpPre <= floor((nS-1)/nHrz)*(1+double(applyFilters))*nHrz + mod(nS-1,nHrz)+1
        xlabel('')
    end
    
    %% parse subject params
    subjColIdx=find(cell2mat(notchParams(1,subjCols))==subjIdcs(s));
    
    if ~isempty(subjColIdx)&&notchParams{1,subjCols(subjColIdx)}~=subjIdcs(s)
        error('wrong column')
    end
    
    if isempty(subjColIdx)||~applyFilters;continue;end
    
    if subjColIdx==numel(subjCols)
        lastCol=size(notchParams,2);
    else
        lastCol=subjCols(subjColIdx+1)-1;
    end
    subjParams=notchParams(2:end,subjCols(subjColIdx):lastCol);
    notchBW=cell2mat(subjParams(1,:));
    nBW=lastCol-subjCols(subjColIdx)+1;
    WoSet=cell(nBW,1);
    for n=1:nBW
        WoTmp=subjParams(2:end,n);
        WoSet{n}=cell2mat(WoTmp(cellfun(@(x) ~ismissing(x), WoTmp)))';
    end
    
%     if subjIdcs(s)==133
%         [amp,freq]=pwelch(double(data),[],[],[],Fs);
%         ampMean=mean(amp,2);
%         logAmp=log10(ampMean);
%         fRes=numel(freq)/range(freq);
%         [b,a]=butter(3,0.5/(fRes/2),'high');
%         ampSpikes=filtfilt(b,a,logAmp);
%         noiseFreqs=freq(round(mean(mask2bounds(ampSpikes>0.5),2)));
%         notchBW(end+1)=0.03;
%         nBW=nBW+1;
%         WoSet(end+1)={noiseFreqs};
%     end
    %% apply filters
    data=double(data);
%     if strcmp(subj,'MG51')
    if subjIdcs(s)==-51
        filtFreq=0.4;
        filtOrd=4;
    else
        filtFreq=0.9;
        filtOrd=3;
    end
        tic;fprintf('antialiasing filter...\n')
        [b,a]=butter(filtOrd,filtFreq,'low');
        data=filtfilt(b,a,data);toc
%     end
    fprintf('notch filters...\n');tic
    Nq = Fs/2;
    for n=1:nBW
        for w=1:numel(WoSet{n})
            Wo = WoSet{n}(w)/Nq;
            if Wo>0.97
                continue
            end
            bw = notchBW(n)/Nq;
            
            [b,a] = iirnotch(Wo,bw);
            data = filtfilt(b,a,data);
        end
    end;toc
    

    
    
    %% plot and save
    subtightplot(nVrt,nHrz,nSpPost,gap,marg_h,marg_w)
    pwelch(data(:,1:5:end),[],[],[],Fs)
    title([subj ' post'])
    if mod(nSpPost,nHrz)~=1
        ylabel('')
    end
    if nSpPost <= floor((nS-1)/nHrz)*(1+double(applyFilters))*nHrz + mod(nS-1,nHrz)+1
        xlabel('')
    end
    disp('saving')
    dataDoub=data;
    data=single(dataDoub);
    save(outFil,'data','label','bpOp','bpLabel','Fs','trigData')
    
    if Fs~=rsFs
        [P,Q]=rat(rsFs/Fs);
        data=resample(dataDoub,P,Q);
        trigData=resample(double(trigData),P,Q,0);
        data=single(data);
        Fs=rsFs;
    end
%     save(outFilRs,'data','Fs','label','bpOp','bpLabel','-append')
    save(outFilRs,'data','Fs','label','bpOp','bpLabel','trigData','-v7.3','-nocompression')
end

savepdf(h,'/home/jgarret/TaskAnalysis/lang/WN/noiseRemoval/prePostNoiseNotch.pdf')

end