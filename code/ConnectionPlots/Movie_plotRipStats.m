function WN_plotRipStats(subjIdcs,avgRef)
close all
if nargin==0
    subjIdcs = [49 51];
end
if nargin<2
    avgRef=false;
end
if avgRef
    apndAvgRef='_avgRef';
else
    apndAvgRef='';
end
parcType='ST',
fRoot='/space/seh8/6/halgdev/projects/jacob/NC_ripple/matFiles/MG/';
datRoot='/space/seh8/6/halgdev/projects/jacob/auditory_ripples/';
subjects=genSubjList(subjIdcs,'MG');
nSubj=numel(subjects);
apndHG='';
% ripFreq=cell(nSubj+1,1);
% ripFreq{end}=[];
%grp={};
ripFreq=cell(1,numel(subjects)+1); %over rips
ripDur=ripFreq;
ripAmp=ripFreq;
ripPer=ripFreq;
ripDens=cell(1,numel(subjects)+1); %over channels


metaParcMap=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls',parcType));
regs=unique(metaParcMap(:,1));
nReg=numel(regs);
regFreqs=cell(nReg,1);
regDurs=cell(nReg,1);
regAmps=cell(nReg,1);
regDens=cell(nReg,1);

nCh=nan(nSubj,1);
chFreqs=cell(nSubj,1);
chFreqMeans=cell(nSubj,1);
subjChParc=cell(nSubj,1);

for s=1:nSubj
    subj=subjects{s};
    labParc=retrieveChanParcels_MG(subj,parcType);
    

    [trl,~,rej]=loadTrials(subj,'WN');
    gTrl=trl(~rej,1);
    
%     try
        load(sprintf('%s/%s_ripple_stats_WN_0_trlSlct%s%s.mat',fRoot,subj,apndAvgRef,apndHG),'rippleStats')
        mobj=matfile(sprintf('%s/altnalysis_1kHz_single/MG/%s/tasks/%s_WN.mat',datRoot,subj,subj));
        nS=size(mobj,'data',1);
    gTrlMsk=bounds2mask(gTrl+[-500 2000],nS);

        
                load(sprintf('%s/WN/rejection/%s_WN_ft_cln.mat',datRoot,subj),'rejChan')
            chan_labels=labParc(:,4);
            chanMap=find(~strcmp(labParc(:,1),'N/A') & ~[rejChan; true(size(labParc,1)-numel(rejChan),1)])';
      
% load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/chansOfInterest/subjCoi_allLH_allRH_200%s.mat',subj),'subjCoi')
%     catch
%         
%         load(sprintf('%s/%s_ripple_stats_FW_0.mat',fRoot,subj),'rippleStats')
%     end
    
    
%     chan_labels=labParc(:,4);
%     label=rippleStats.chanLabels';
%     label=cellfun(@(x) [x(regexp(x,'[A-Z]')) sprintf('%02d',str2num(x(regexp(x,'[0-9]'))))],label,'uni',0);
%     switch subj
%             case 'NY233'
%                 label=strrep(label,'RM04','RMT04'); %typo
%             case 'NY190'
%                 label{114}='RFT10';
%                 label(61:64)={'GB59';'GB60';'GB63';'GB64'};
%             case 'NY86'
%                 label=strrep(label,'TO','IO');
%     end
%     
%     %load(sprintf('/home/jgarret/TaskAnalysis/lang/FW/chansOfInterest/%s_regChans.mat',subj));
%     

    
%     chanMap=cellfun(@(x) find(strcmp(x,chan_labels)),subjCoi.(subj));
%     chanMap=chanMap(arrayfun(@(x) any(strcmp(labParc{x,1},strRoi)),chanMap));
    
    gChanParc.(subj)=labParc(chanMap,:);
    
%     nCh(s)=numel(subjCoi.(subj));
    nCh(s)=numel(chanMap);
%     if ~all(arrayfun(@(x) strcmp(gChanParc{x,4},subjCoi.(subj){x}), 1:nCh(s)))
%         error('label mismatch')
%     end
    subjChParc{s}=gChanParc.(subj);
    
    
    chFreqs{s}=cell(nCh(s),1);
    chFreqMeans{s}=nan(nCh(s),1);
    
    
    fn=fieldnames(rippleStats);
    for f=1:numel(fn)
        if size(rippleStats.(fn{f}),2)==numel(chan_labels)
            rippleStats.(fn{f})=rippleStats.(fn{f})(chanMap);
        end
    end
    
    chanRegs=labParc(chanMap,1);
    chanHemi=labParc(chanMap,3);
    
    gRip=cell(nCh(s),1);
    for ch=1:nCh(s)
        gRip{ch}=gTrlMsk(rippleStats.locs{ch});
        rippleStats.oscFreq{ch}(~gRip{ch})=[];
        rippleStats.duration{ch}(~gRip{ch})=[];
        rippleStats.rippleAmp{ch}(~gRip{ch})=[];
    end
    ripFreq{s}=[rippleStats.oscFreq{:}];
    for r=1:nReg
        regFreqs{r} = [regFreqs{r} rippleStats.oscFreq{strcmp(chanRegs,regs{r})&strcmp(chanHemi,'lh')}];
        regDurs{r} = [regDurs{r} rippleStats.duration{strcmp(chanRegs,regs{r})&strcmp(chanHemi,'lh')}];
        regAmps{r} = [regAmps{r} rippleStats.rippleAmp{strcmp(chanRegs,regs{r})&strcmp(chanHemi,'lh')}];
        regDens{r} = [regDens{r} cell2mat(rippleStats.density(strcmp(chanRegs,regs{r})&strcmp(chanHemi,'lh')))];
        
    end
    for ch=1:nCh(s)
        chFreqs{s}{ch}=rippleStats.oscFreq{ch};
        chFreqMeans{s}(ch)=mean(chFreqs{s}{ch});
    end
    %ripFreq{end}=[ripFreq{end} ripFreq{s}];
    
    ripDur{s}=[rippleStats.duration{:}];
    %ripDur{end}=[ripDur{end} ripDur{s}];
    
    ripAmp{s}=[rippleStats.rippleAmp{:}];
%     if ((abs(subjIdcs(s))>111 && abs(subjIdcs(s))<442) || subjIdcs(s)==-89) && subjIdcs(s)~=133
%         ripAmp{s}=ripAmp{s}/10;
%     end
    %ripAmp{end}=[ripAmp{end} ripAmp{s}];
    
    
    ripDens{s}=cell2mat(rippleStats.density)/(sum(gTrlMsk)/nS);
%     ripDens{s}=cellfun(@(x) (sum(x)/sum(gTrlMsk))*60*1000,gRip)';
    
    
%     subjFreq=[rippleStats.oscFreq{:}];
%     nRip=numel(subjFreq);
%     ripFreq=[ripFreq subjFreq];
%     
%     ripDur=[ripDur rippleStats.duration{:}];
%     
%     subjAmp=[rippleStats.rippleAmp{:}];
%     if (abs(subjIdcs(s))>111 && abs(subjIdcs(s))<442) || subjIdcs(s)==-89
%         subjAmp=subjAmp/10;
%     end
%     ripAmp=[ripAmp subjAmp];
%     
%     grp=[grp repmat({subj},[1 nRip])];
end



regFreqMeans=cellfun(@(x) mean(x),regFreqs);
regFreqLabs=regs;
save(sprintf('%s/coripDat/regRipFreqs_%s.mat',datRoot,parcType),'regFreqLabs','regFreqMeans','regFreqs')
% chFreqSubs=subjects;
% save('/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/chanRipMeans/chanRipMeans.mat',...
%     'chFreqs','chFreqMeans','chFreqSubs','gChanParc')

ripFreq(end)={[ripFreq{:}]};
ripDur(end)={[ripDur{:}]};
ripAmp(end)={[ripAmp{:}]};
ripDens(end)={[ripDens{:}]};
% ripDens{end}(ripDens{end}>60)=nan;

subselReg=[1 2 5 6 7 8 11 12 13 14 16 17 19 20 21 23 24 25 26 27 28 29];
% subselReg=subselReg([10 6 12 7 3 5 9 18 1 13 19 2 20 11 4 8 16 17 15 14]);

subselReg=subselReg([6 12 7 3 5 9 22 18 21 1 13 19 2 20 11 4 8 17 15 14]);
regFreqs=regFreqs(subselReg);
regDurs=regDurs(subselReg);
regAmps=regAmps(subselReg);
regDens=regDens(subselReg);
regs=regs(subselReg);
nReg=numel(subselReg);
regFreqs(end+1)={[regFreqs{:}]};
regDurs(end+1)={[regDurs{:}]};
regAmps(end+1)={[regAmps{:}]};
regDens(end+1)={[regDens{:}]};

% subjects=[subjects {'all'}];
subjects=[arrayfun(@(x) sprintf('P%d',x),1:nSubj,'uni',0) {'all'}];
nSubj=nSubj+1;
subsel=[nSubj 1:nSubj-1];
%%
gap=[0.04 0.1];
marg_h=[0.09 0.07];
marg_w=[0.1 0.04];
% nPlt=nSubj+1;
% nHrz=5;
nPlt=4;

%h1=figure('Position',  [261 1 1100 1000])
h1=figure('Position',  [261 1 700 1000])
subtightplot(nPlt,1,1,gap,marg_h,marg_w)
%boxplot(ripFreq, grp)
[~,L]=violin(ripFreq);%,'xlabel',subjects)
xticks(1:nSubj)
xticklabels(subjects(subsel))
ylabel('freq. (Hz)')
ylim([60 120])
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(ripFreq{end}),median(ripFreq{end})))
ax=gca;
ax.FontSize=8;
%axis('square')
% set(ax,'linewidth',1.5);
set(L,'FontSize',8)





subtightplot(nPlt,1,2,gap,marg_h,marg_w)
% boxplot(ripDur, grp)
[~,L]=violin(ripDur(subsel));%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nSubj)
xticklabels(subjects(subsel))
ylim([0 600])
ylabel('dur. (ms)')
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(ripDur{end}),median(ripDur{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);

subtightplot(nPlt,1,3,gap,marg_h,marg_w)
% boxplot(ripAmp, grp)
[~,L]=violin(ripAmp(subsel));%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nSubj)
xticklabels(subjects(subsel))
ylabel('amp. (uV)')
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(ripAmp{end}),median(ripAmp{end})))
ax=gca
ax.FontSize=8;
% set(ax,'linewidth',1.5);

% subtightplot(nPlt,1,4,gap,marg_h,marg_w)
% % boxplot(ripAmp, grp)
% violin(ripPer(1:end-1));%,'xlabel',subjects)
% set(gca,'YScale','log')
% xticks(subjOrdr)
% xticklabels(subjects)
% ylabel('inter-ripple period (ms)')

subtightplot(nPlt,1,4,gap,marg_h,marg_w)
% boxplot(ripAmp, grp)
[~,L]=violin(ripDens(subsel));%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nSubj)
xticklabels(subjects(subsel))
ylabel('rips/min')
ylim([0 max(ylim)])
title(sprintf('per channel, mean=%.1f median=%.1f',nanmean(ripDens{end}),nanmedian(ripDens{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);


set(h1,'color','w');% make figure white
% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');
% wysiwyg export
addpath('/home/bqrosen/matlab/')
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved
          rmpath('/home/bqrosen/matlab/')

savepdf(h1,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/ripSummary%s.pdf',apndHG))
saveas(h1,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/ripSummary%s.png',apndHG))

%%
load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/regShrt_%s.mat',parcType))
regs=[regShrt(subselReg) {'all'}];
nReg=nReg+1;
h2=figure('Position',  [261 1 700 1000])
% h2=figure('Position',  [1 1 1000 900])
subtightplot(nPlt,1,1,gap,marg_h,marg_w)
%boxplot(regFreqs, grp)
[~,L]=violin(regFreqs');%,'xlabel',subjects)
xticks(1:nReg)
xticklabels(regs)
ylabel('freq. (Hz)')
ylim([60 120])
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regFreqs{end}),median(regFreqs{end})))
ax=gca;
ax.FontSize=8;
%axis('square')
% set(ax,'linewidth',1.5);
set(L,'FontSize',8)





subtightplot(nPlt,1,2,gap,marg_h,marg_w)
% boxplot(regDurs, grp)
[~,L]=violin(regDurs');%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylim([0 600])
ylabel('dur. (ms)')
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regDurs{end}),median(regDurs{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);

subtightplot(nPlt,1,3,gap,marg_h,marg_w)
% boxplot(regAmps, grp)
[~,L]=violin(regAmps');%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylabel('amp. (uV)')
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regAmps{end}),median(regAmps{end})))
ax=gca
ax.FontSize=8;
% set(ax,'linewidth',1.5);

% subtightplot(nPlt,1,4,gap,marg_h,marg_w)
% % boxplot(regAmps, grp)
% violin(ripPer(1:end-1));%,'xlabel',subjects)
% set(gca,'YScale','log')
% xticks(subjOrdr)
% xticklabels(subjects)
% ylabel('inter-ripple period (ms)')

subtightplot(nPlt,1,4,gap,marg_h,marg_w)
% boxplot(regAmps, grp)
[~,L]=violin(regDens');%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylabel('rips/min')
ylim([0 max(ylim)])
title(sprintf('per channel, mean=%.1f median=%.1f',nanmean(regDens{end}),nanmedian(regDens{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);


set(h2,'color','w');% make figure white
% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');
% wysiwyg export
addpath('/home/bqrosen/matlab/')
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved
          rmpath('/home/bqrosen/matlab/')

savepdf(h2,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/ripRegSummary%s.pdf',apndHG))
saveas(h2,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/ripRegSummary%s.png',apndHG))

load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/regShrt_%s.mat',parcType))
regs=[regShrt(subselReg) {'all'}];
nReg=nReg+1;
h2=figure('Position',  [261 1 700 1000])
% h2=figure('Position',  [1 1 1000 900])
subtightplot(nPlt,1,1,gap,marg_h,marg_w)
%boxplot(regFreqs, grp)
[~,L]=violin(regFreqs');%,'xlabel',subjects)
xticks(1:nReg)
xticklabels(regs)
ylabel('freq. (Hz)')
ylim([60 120])
title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regFreqs{end}),median(regFreqs{end})))
ax=gca;
ax.FontSize=8;
%axis('square')
% set(ax,'linewidth',1.5);
set(L,'FontSize',8)


%%
regs={'Early Visual','Other'};
regDurs={[regDurs{1:2}], [regDurs{3:end-1}]};
regAmps={[regAmps{1:2}], [regAmps{3:end-1}]};
regDens={[regDens{1:2}], [regDens{3:end-1}]};
regFreqs={[regFreqs{1:2}], [regFreqs{3:end-1}]};
%%
saveSheet=cell(3,4);
h2=figure('Position',  [261 1 400 400])
% h2=figure('Position',  [1 1 1000 900])
subtightplot(2,2,1,gap,marg_h,marg_w)
%boxplot(regFreqs, grp)
[~,L]=violin(regFreqs);%,'xlabel',subjects)
xticks(1:nReg)
xticklabels(regs)
ylabel('freq. (Hz)')
ylim([60 120])
% title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regFreqs{end}),median(regFreqs{end})))
ax=gca;
ax.FontSize=8;
%axis('square')
% set(ax,'linewidth',1.5);
set(L,'FontSize',8)
saveSheet{1,1}=nanmean(regFreqs{1});
saveSheet{2,1}=nanmean(regFreqs{2});
[~,saveSheet{3,1}]=ttest2(regFreqs{1},regFreqs{2});

subtightplot(2,2,2,gap,marg_h,marg_w)
% boxplot(regDurs, grp)
[~,L]=violin(regDurs);%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylim([0 600])
ylabel('dur. (ms)')
% title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regDurs{end}),median(regDurs{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);
saveSheet{1,2}=nanmean(regDurs{1});
saveSheet{2,2}=nanmean(regDurs{2});
[~,saveSheet{3,2}]=ttest2(regDurs{1},regDurs{2});

subtightplot(2,2,3,gap,marg_h,marg_w)
% boxplot(regAmps, grp)
[~,L]=violin(regAmps);%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylabel('amp. (uV)')
% title(sprintf('per ripple, mean=%.1f median=%.1f',mean(regAmps{end}),median(regAmps{end})))
ax=gca
ax.FontSize=8;
% set(ax,'linewidth',1.5);
saveSheet{1,3}=nanmean(regAmps{1});
saveSheet{2,3}=nanmean(regAmps{2});
[~,saveSheet{3,3}]=ttest2(regAmps{1},regAmps{2});

% subtightplot(nPlt,1,4,gap,marg_h,marg_w)
% % boxplot(regAmps, grp)
% violin(ripPer(1:end-1));%,'xlabel',subjects)
% set(gca,'YScale','log')
% xticks(subjOrdr)
% xticklabels(subjects)
% ylabel('inter-ripple period (ms)')

subtightplot(2,2,4,gap,marg_h,marg_w)
% boxplot(regAmps, grp)
[~,L]=violin(regDens);%,'xlabel',subjects)
set(L,'Visible','off')
xticks(1:nReg)
xticklabels(regs)
ylabel('rips/min')
ylim([0 max(ylim)])
% title(sprintf('per channel, mean=%.1f median=%.1f',nanmean(regDens{end}),nanmedian(regDens{end})))
ax=gca;
ax.FontSize=8;
% set(ax,'linewidth',1.5);
saveSheet{1,4}=nanmean(regDens{1});
saveSheet{2,4}=nanmean(regDens{2});
[~,saveSheet{3,4}]=ttest2(regDens{1},regDens{2});



set(h2,'color','w');% make figure white
% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');
% wysiwyg export
addpath('/home/bqrosen/matlab/')
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved
          rmpath('/home/bqrosen/matlab/')

% savepdf(h2,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/FW/ripRegSummary_earlyVisVsOther%s.pdf',apndHG))
% saveas(h2,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/FW/ripRegSummary_earlyVisVsOther%s.png',apndHG))

saveGraphic(h2,sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/ripRegSummary_earlyVisVsOther%s.pdf',apndHG),[],600)

saveSheet=[{'','freq','dur','amp','rate'};[{'visual';'other';'pval'} saveSheet]];
writecell(saveSheet,'~/TaskAnalysis/lang/FW/earlyVisCoripTest.xlsx')


    
end

