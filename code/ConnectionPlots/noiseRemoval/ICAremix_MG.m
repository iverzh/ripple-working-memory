function ICAremix_MG(icaMethod,subj,trlSlct,filtDat)

% icaMethod='fastica';
% % icaMethod='runica';
% % subj='NY133';
% % subj='NY255';
% subj='NY86';
% % subj='NY70';
if nargin<3
    trlSlct=false
end
if trlSlct
    trlSlctApnd='_trlSlct';
else
    trlSlctApnd='';
end
if nargin<4
    filtDat=false
end
if filtDat
    filtApnd='_filtDat';
else
    filtApnd='';
end

% datDir='/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/';
datDir='/space/seh8/5/halgdev/projects/jgarret/altnalysis_native_fs/MG/';

load(sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/WN/ica/%s_ICA_%s%s%s.mat',subj,icaMethod,trlSlctApnd,filtApnd),'comp','labParcs','deadChans')
% lab=labParcs(:,4);
load(sprintf('%s/%s/tasks/%s_WN_preICA.mat',datDir,subj,subj),'label')
lab=label;


wts=comp.topo;

nCh=size(wts,1);
nComp=size(wts,2);

compDat=comp.trial{1}';

% Fs=1000;
Fs=comp.fsample;




% trl=loadTrialInfo_v2(subj,0);
% trl=trl(trl(:,1)>3,1);
%
% load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/v2/%s_FW_ft_cln_v2.mat',subj), 'rejChan');
% chanIdx=find(~rejChan);
%%
%{
rip_mask=any(rip_mask,1);
ripIdx=find(rip_mask);

rip_bounds_mask=any(rip_bounds_mask,1);
rip_bounds=mask2bounds(rip_bounds_mask);
ripIdx=mean(rip_bounds,2);
%%
%}


%% RUN THIS SECTION (e.g. in debug mode) TO PLOT THINGS
if 0
    
    %%
%     [trl,~,rej]=loadTrials('MG51','WN');
%     [b,a]=butter(3,[70 190]/(Fs/2),'bandpass');
%     gTrl=trl(~rej,:);
%     tTrl=gTrl(:,1)/2;
    nPerPg=15;
    rsFctr=5;
    for n=1:ceil(nComp/nPerPg)
        figure(n)
        clf
        subtightplot(1,1,1,[0 0],[0 0],[0.02 0.001])
        set(gcf,'Position',get(0,'Screensize'))
        idxVec=(1:nPerPg)+(nPerPg*(n-1));
        idxVec(idxVec>nComp)=[];
        lamplot(compDat(1:rsFctr:end,idxVec),1,1)
        yticklabels(idxVec)
        %     waitforbuttonpress
        hold on
%             for t=1:numel(tTrl)
%                 line(tTrl(t)*[1 1], ylim)
%             end
        
    end
    %%
    for n=1%1:ceil(nComp/nPerPg)
        figure(n)
               hold on
            for t=1:numel(tTrl)
                line(tTrl(t)*[1 1], ylim)
            end 
    end
    
    %%
    
    [trl,~,rej]=loadTrials('MG51','WN');
    [b,a]=butter(3,[70 190]/(Fs/2),'bandpass');
    gTrl=trl(~rej,:);
    tTrl=gTrl(:,1)/2;
    trlMsk=(tTrl+(-500:1000))';
    tmp=permute(compDat,[1 3 2]);
    tmp2=reshape(tmp(trlMsk,:),[size(trlMsk),124]);
    tmp3=squeeze(mean(tmp2,2));
    tmp4=abs(hilbert(filtfilt(b,a,double(tmp3))));
    tmp5=smoothdata(tmp4,1,'gaussian',30);
    
    nPerPg=15;
    rsFctr=5;
    for n=1:ceil(nComp/nPerPg)
        figure(n)
        clf
        subtightplot(1,1,1,[0 0],[0 0],[0.02 0.001])
        set(gcf,'Position',get(0,'Screensize'))
        idxVec=(1:nPerPg)+(nPerPg*(n-1));
        idxVec(idxVec>nComp)=[];
        lamplot(-(tmp5(:,idxVec)-mean(tmp5(:,idxVec),1)),1,1)
        yticklabels(idxVec)
        %     waitforbuttonpress
        hold on
        %     for r=1:numel(ripIdx)
        %         line(ripIdx(r)*[1 1], ylim)
        %     end
    end
end

%}
%%
%{
nPerPg=15;
rsFctr=5;
for n=1:ceil(nComp/nPerPg)
    figure(n+10)
    clf
    subtightplot(1,1,1,[0 0],[0 0],[0.02 0.001])
    set(gcf,'Position',get(0,'Screensize'))
    idxVec=(1:nPerPg)+(nPerPg*(n-1));
    idxVec(idxVec>nComp)=[];
    lamplot(data(1:rsFctr:end,idxVec),1,1)
    yticklabels(label(idxVec))
%     waitforbuttonpress
    hold on
%     for r=1:numel(ripIdx)
%         line(ripIdx(r)*[1 1], ylim)
%     end
end
end
%}

%%
switch subj
    case 'MG49'
        trgIdx=1;
        switch icaMethod
            case 'runica'
                artComps=[76];
        end
    case 'MG51'
%         switch icaMethod
%             case 'fastica'
%                 artComps=[5 31 44 66 71 99 100 104 111 117]; %obvious artifacts (104 also has lots of line noise)
%                 artComps=[artComps 45 56]; %components with lots of line noise
%                 artComps=[artComps 113]; %components with large HF component
% %                 artComps=[artComps 36 35]; %components with low weight variance
%                 artComps=[artComps 23 15 45 3 62 105]; %components with high 9th percentile weight
%             case 'runica'
%                 % approach: ignore large artifacts in favor of removing
%                 % them later, instead get rid of coherent high frequency
%                 % noise
%                 artComps=[1 2 9 11 17 ... high frequency oscillatory (line?) noise
%                     12 ... high frequency with low variance in component weights
%                     37 ... high frequency with lowish variance in component weights, but may not be important
%                     3 5 6 7 8 12 13 15 ... giant weight [wtVals(:,1),wtVals(:,2)]=sort(median(abs(wts),1),'descend');
%                     ];
%                     

        switch icaMethod
            case 'fastica'
                [trialDat,~,rej]=loadTrials(subj,'WN');
                trialDat(rej,:)=[];
                rmDat = ~bounds2mask(trialDat(:,1)+[0 2000],length(compDat));
%                 % approach: ignore large artifacts in favor of removing
%                 % them later, instead get rid of coherent high frequency
%                 % noise
                artComps=[13 16 ... high frequency oscillatory (line?) noise
                    55 115 ... super high frequency unvarying noise
                    13 70 21 81 29 73 44 55 16 ... median deviation very small (<=0.05)...
                        ... relative to median (median based variance), formula:
                        ... [tmp(:,1),tmp(:,2)]=sort(abs(mad(wts,1,1)./median(wts,1)));
                    16 73 13 ... super high frequency: [b,a]=butter(3,200/(Fs/2),'high'); ...
                        ... [tmp(:,1),tmp(:,2)]=sort(mean(abs(hilbert(filtfilt(b,a,double(compDat)))),1),'descend');
                        ... or [tmp(:,1),tmp(:,2)]=sort(mean(abs(hilbert(filtfilt(b,a,double(compDat(~rmDat,:))))),1),'descend');
                     ... high frequency with lowish variance in component weights, but may not be important
                     13 81 53 16 70 21 73 29 44 ... giant median weight [tmp(:,1),tmp(:,2)]=sort(median(abs(wts),1),'descend');
                     9 ... 14 ... components well represented in G1-53 and not G57-58
                    ];
                artComps=unique(artComps);
            case 'runica'
                artComps=[1 2 4 8 9 11 12 17 ... high frequency oscillatory (line?) noise
                     ... super high frequency unvarying noise
                     7 5 9 3 12 ... median deviation very small (<=0.05)...
                        ... relative to median (median based variance), formula:
                        ... [tmp(:,1),tmp(:,2)]=sort(abs(mad(wts,1,1)./median(wts,1)));
                     12 9 1 ... super high frequency: [b,a]=butter(3,200/(Fs/2),'high'); ...
                        ... [tmp(:,1),tmp(:,2)]=sort(mean(abs(hilbert(filtfilt(b,a,double(compDat)))),1),'descend');
                        ... or [tmp(:,1),tmp(:,2)]=sort(mean(abs(hilbert(filtfilt(b,a,double(compDat(~rmDat,:))))),1),'descend');
                     1 3 5 6 7 9 12 13 14 ... giant median weight [tmp(:,1),tmp(:,2)]=sort(median(abs(wts),1),'descend');
                     ... components well represented in G1-53 and not G57-58
                    ];
                artComps=sort(unique(artComps));
                artComps=[];
                
        end
    case 'MG51v2'
        switch icaMethod
            case 'runica'
                
            case 'fastica'
                
        end
end

if 0
    fprintf('removing line noise... ');tic
    Nq = Fs/2;
    notchBW=[3];% 20 2];
    WoSet={[60:60:Nq], [495],[499.9]};
    for n=1:numel(notchBW)
        for w=1:numel(WoSet{n})
            Wo = WoSet{n}(w)/Nq;
            bw = notchBW(n)/Nq;
    
            [b,a] = iirnotch(Wo,bw);
            compDat = single(filtfilt(b,a,double(compDat)));
        end
    end
    toc
end


wts=comp.topo;
wts(:,artComps)=0;

data=zeros(length(comp.trial{1}),numel(lab));
% data(:,~deadChans) = (wts * comp.trial{1})';
data(:,~deadChans) = (wts * compDat')';

% Fs=1000;
%
% Nq = Fs/2;
% notchBW=[3 20 2];
% WoSet={[60:60:Nq], [495],[499.9]};
% for n=1:numel(notchBW)
%     for w=1:numel(WoSet{n})
%         Wo = WoSet{n}(w)/Nq;
%         bw = notchBW(n)/Nq;
%
%         [b,a] = iirnotch(Wo,bw);
%         data = single(filtfilt(b,a,double(data)));
%     end
% end

%%
%{
wtsNorm=sqrt(wts.^2./mean(wts.^2,2));
figure();
scatter(repmat((1:nCh)',[nComp 1]), wtsNorm(:))
hold on
scatter(1:nCh,max(wtsNorm,[],2),[],'r')
legend('other comps','max weighted comp')
ylabel('component power')
xlabel('channel')

% wtTr=wts';
% figure();
% scatter(repmat((1:nComp)',[nCh 1]), wtTr(:).^2)
% hold on
% scatter(1:nComp,max(wtTr.^2,[],2),[],'r')
% legend('other chans','max chan weight')
% ylabel('component power')
% xlabel('channel')
%}



%%
%{
[b,a]=butter(4,[70 100]/Nq);
rbDat=filtfilt(b,a,double(data));
figure()
set(gcf,'Position',get(0,'Screensize'))
intvl=3;
for t=1:numel(trl)
    dat_tmp=rbDat(trl(t)+(1:1000),chanIdx(1:intvl:end));
    dat_tmp=dat_tmp./std(dat_tmp,0,1);
    lamplot(dat_tmp,1,0)
    xlabel('ms')
    yticklabels(lab(chanIdx(1:intvl:end)))
    waitforbuttonpress
end
%}
%%
%
% sampIdx=[950:980];
%
% [b,a]=butter(4,[70 100]/Nq);
% rbComp_tmp=abs(hilbert(filtfilt(b,a,double(compDat(trl(t)+(1:1000),:)))));
% rbComp_tmp=rbComp_tmp.*median(wts,1);
% [val,idx]=sort(mean(rbComp_tmp(sampIdx,:),1)'.^2,'descend');
%}
%%

if exist('trgIdx','var')
    mobj=matfile(sprintf('%s/%s/tasks/%s_WN_preICA.mat',datDir,subj,subj));
    data(:,trgIdx)=mobj.data(:,trgIdx);
end
if ~exist(sprintf('%s/%s/tasks/%s_WN_noNotch.mat',datDir,subj,subj),'file')
    warning('no noNotch file, copying preICA file before overwriting data')
    unix(sprintf('cp %s %s',sprintf('%s/%s/tasks/%s_WN_preICA.mat',datDir,subj,subj),...
        sprintf('%s/%s/tasks/%s_WN_noNotch.mat',datDir,subj,subj)));
end
fprintf('saving... ');tic
save(sprintf('%s/%s/tasks/%s_WN_noNotch.mat',datDir,subj,subj),'data','-append')
toc;fprintf('noNotch file copy saved, saving to main file...');tic
save(sprintf('%s/%s/tasks/%s_WN.mat',datDir,subj,subj),'data','-append')
toc
% save(sprintf('/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/%s/tasks/%s_WN.mat',subj,subj),'data','-append')

end
