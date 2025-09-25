function WN_regPLV

% subjIdcs=[-240 233 226 190 133 130 101 -89 86 68];
subjIdcs=[51];
close all
plvThresh=0.2;
avgRef=0;
coripCenterPhase=0;
plotPerReg=1;
plotSingleRegPr=0;
binSz=100;
parcType='ST';
brainImgDir='/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/';
if coripCenterPhase
    append='_coripCenterPhase';
else
    append='';
end
if avgRef
    apndAvgRef='_avgRef';
else
    apndAvgRef='';
end
alpha=0.01;
fRoot='/space/seh8/6/halgdev/projects/jacob/auditory_ripples/coripDat/';
% load(sprintf('%s/subjRipDat/subjRipDat_200%s_allLH.mat',fRoot,apndAvgRef))
% load(sprintf('%s/coripDat/coripDat_split%s_allLH_200.mat',fRoot,apndAvgRef))
outDir=sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/figs/WN/regPLV%s/',append);
if ~exist(outDir,'dir');unix(sprintf('mkdir -p %s',outDir));end
%%
subjects=genSubjList(subjIdcs,'MG');

nS=numel(subjects);

coripPhDist=cell(nS,1);
coripPhAgg=[];
subjCh=cell(nS,1);
% tmp=[];
for s=1:nS
    subj=subjects{s}
    load(sprintf('%s/%s/phRip_%s.mat',fRoot,subj,subj),'gCh','labParc','phRip')
%     [trl,resp]=loadTrialInfo_v2(subj,0);
    [trl,~,rej]=loadTrials(subj,'WN');
    gTrl=trl(~rej,1);


    ripSubj=phRip;
    ripSubj=cos(ripSubj)+1i*sin(ripSubj);
    nRip=sum(~isnan(ripSubj),2);
%         trlMsk=bounds2mask([0 1500]+gTrl,length(nRip)
        trlMsk=bounds2mask([-200 800]+gTrl,length(nRip));
   
%     gResp=resp(resp(:,2)==2&resp(:,4)>0,1);
%     respMsk=bounds2mask([-400 12]+gResp,length(ripSubj));
    
    labParc=retrieveChanParcels_MG(subj,parcType);
    subjCh{s}=labParc(gCh,[1 4]);
    nCh{s}=size(subjCh{s},1);
    
    
        trlMsk=bounds2mask([0 1500]+gTrl,length(nRip));
%         trlMsk=bounds2mask([200 400]+gTrl,length(nRip));

%     trlMsk=respMsk;
    nRip(~trlMsk)=nan;
    
    ripBnds=mask2bounds(nRip>=2);
    if isempty(ripBnds)
        continue
    end
    ripBnds(ripBnds(:,2)-ripBnds(:,1)<25,:)=[];
    coripPhDist{s}=nan(size(ripBnds,1),nCh{s},nCh{s});
    for r=1:size(ripBnds,1)

        ch=find(~all(isnan(ripSubj(ripBnds(r,1):ripBnds(r,2),:)),1));
        ripPhTmp=nan(nCh{s});
        for c1=ch
            for c2=ch
                if c1==c2
                    continue
                end
                mskPr=all(~isnan(ripSubj(ripBnds(r,1):ripBnds(r,2),[c1 c2])),2);
                bndPr=mask2bounds(mskPr);
                [ripLenPr,ripIdx] = max(diff(bndPr,1,2));
                if isempty(ripLenPr)||ripLenPr<25
                    ripPhTmp(c1,c2)=nan;
                else
%                     tmp=numel(diff(bndPr,1,2)>25);
                    prCoripIdx=(bndPr(ripIdx,1):bndPr(ripIdx,2))+ripBnds(r,1)-1;
%                     prCoripIdx=(round(ripLenPr/4+bndPr(ripIdx,1)):round(bndPr(ripIdx,2)-ripLenPr/4))+ripBnds(r,1)-1;
                    phCh1=angle(ripSubj(prCoripIdx,c1));
                    phCh2=angle(ripSubj(prCoripIdx,c2));
                    angTmp=phCh1-phCh2;
                    if coripCenterPhase
                        angTmp=angTmp(round(numel(angTmp)/2));
                    end
                    ripPhTmp(c1,c2)=angle(mean(cos(angTmp)+1i*sin(angTmp)));
                    %fprintf('')
                end
            end
        end
        
        coripPhDist{s}(r,:,:)=ripPhTmp;
    end
    
end


%% remap coripples to channel pairs

prNames=cell(nS,1);
prRegs = cell(nS,1);
prMap=cell(nS,1);
plvPr = cell(nS,1);
nPr=nan(nS,1);
for s=1:nS
    nPr(s)= nCh{s}*(nCh{s}-1)/2;
    prNames{s}=cell(nPr(s),1); %channel pair names
    prRegs{s}=cell(nPr(s),1); %region pair names
    prMap{s}=nan(nPr(s),2);%map from pairs to channels
    lastIdx=0;
    for c=1:nCh{s}
        for t=c+1:nCh{s}
            prNames{s}{c+t+lastIdx-2}=[subjCh{s}{c,2} ' - ' subjCh{s}{t,2}];
            prRegs{s}{c+t+lastIdx-2}=[subjCh{s}{c,1} ' - ' subjCh{s}{t,1}];
            prMap{s}(c+t+lastIdx-2,:)=[c t];
        end
        lastIdx=lastIdx+t-c-2;
    end
end

prRips=struct();
for s=1:nS
    subj=subjects{s};
    prRips.(subj)=cell(nPr(s),1);
    for p=1:nPr(s)
        
        prChans=prMap{s}(p,:);
        if ~isempty(coripPhDist{s})
%             prRips.(subj){p}=coripPhDist{s}(all(coripChans{s}(:,prChans),2),prMap{s}(p,1),prMap{s}(p,2));
            prRips.(subj){p}=coripPhDist{s}(:,prMap{s}(p,1),prMap{s}(p,2));
            prRips.(subj){p}(isnan(prRips.(subj){p})) = [];
        end
    end
end

% for s=1:nS
%     plvPr{s} = nan(nCh{s},nCh{s});
%     for c1=1:nCh{s}
%         for c2=1:nCh{s}
%             if ~isempty(coripPhDist{s})
%                 distTmp = coripPhDist{s}(all(coripChans{s}(:,[c1 c2]),2),c1,c2);
%             end
%             distTmpAll=vertcat(distTmp{:});
%             if numel(distTmpAll)>=20
%                 plvPr{s}(c1,c2) = abs(mean(1i*sin(distTmpAll) + cos(distTmpAll)));
%             end
%         end
%     end
% end

%% find pairwise frequency differences
%{
load('/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/chanRipMeans/chanRipMeans.mat','chFreqs','chFreqMeans','chFreqSubs')
rmSubj=cellfun(@(x) any(strcmp(x,{'NY311','NY297','NY111'}')),chFreqSubs);
chFreqSubs(rmSubj)=[];
chFreqMeans(rmSubj)=[];
chFreqs(rmSubj)=[];
subjRemap=cellfun(@(x) find(strcmp(subjects,x)),chFreqSubs);
chFreqs=chFreqs(subjRemap);
chFreqMeans=chFreqMeans(subjRemap);
chFreqDiff=cell(nS,1);
chFreqDiffPair=cell(nS,1);
for s=1:nS
    chFreqDiff{s}=abs(chFreqMeans{s}-chFreqMeans{s}');
    chFreqDiffPair{s}=chFreqDiff{s}(prMap{s}(:,1),prMap{s}(:,2));
end
%}

%%

prRipsAll=cell(sum(nPr),1);
prPLVall=nan(sum(nPr),1);
prFreqDiffsAll=nan(sum(nPr),1);

cplx=@(x) cos(x)+1i*sin(x);

nRipPr = nan(sum(nPr),1);
prRipMean = nan(sum(nPr),1);
prPLV = nan(sum(nPr),1);
pPrUni = nan(sum(nPr),1);

for s=1:nS
    idxSubj=1:nPr(s);
    idxAll=idxSubj;
    if s>1
        idxAll=idxAll+sum(nPr(1:(s-1)));
    end
    for p=idxSubj
        prRipsAll(idxAll(p))=prRips.(subjects{s})(p);
        prPLVall(idxAll(p))=abs(mean(cplx(prRipsAll{idxAll(p)})));
        %prFreqDiffsAll(idxAll(p))=chFreqDiffPair{s}(p);
    end
end

for p=1:numel(prRipsAll)
    nRipPr(p) = numel(prRipsAll{p});
    if nRipPr(p)
        prRipMean(p) = circ_mean(prRipsAll{p});
        pPrUni(p) = circ_rtest(prRipsAll{p});
    end
end

pPrUni(nRipPr<20) = nan;
pPrUni(~isnan(pPrUni))=mafdr(pPrUni(~isnan(pPrUni)),'BHFDR',true);
sig=pPrUni<alpha;


%{
h0=figure('Position',get(0,'Screensize')/2);
scatter(abs(prFreqDiffsAll(sig&nRipPr>20)),prPLVall(sig&nRipPr>20),2.5,'filled','r');hold on
scatter(abs(prFreqDiffsAll(~sig&nRipPr>20)),prPLVall(~sig&nRipPr>20),2.5,'filled','b')
mdl=fitlm(abs(prFreqDiffsAll(nRipPr>20)),prPLVall(nRipPr>20));
m=mdl.Coefficients.Estimate(2);
b=mdl.Coefficients.Estimate(1);
p=mdl.Coefficients.pValue(2);
rsq=mdl.Rsquared.Adjusted;
title(sprintf('pairs with N_{corip}>20, b=%.5f R^{2}=%.4f p=%.4f',m,rsq,p))
line(xlim,xlim*m + b,'color','g')
xlabel('absolute frequency difference')
ylabel('PLV')
legend(sprintf('PLV p<%.2f',alpha),sprintf('PLV p>%.2f',alpha))
prettifyPlot(h0,gca)
saveGraphic(h0,outDir,'freqDiffVplv')

h1=figure('Position',get(0,'Screensize')/2);
prFreqDiffsTmp=abs(prFreqDiffsAll(nRipPr>20));
sigTmp=sig(nRipPr>20);
boxplot(prFreqDiffsTmp,~sigTmp)
xticklabels({'PLV p<0.05','PLV p>0.05'})
[h,p,~,stats]=ttest2(prFreqDiffsTmp(sigTmp),prFreqDiffsTmp(~sigTmp));
title(sprintf('p=%.10f',p))
ylim([0 6])
ylabel('absolute frequency difference (Hz)')
prettifyPlot(h1,gca)
saveGraphic(h1,outDir,'freqDiffbySigPLV')
%}


sigPrMean = prRipMean(sig);
% sigPrMean = prRipMean(nRipPr>20);
sigPrMeanMean = circ_mean(sigPrMean);
%% plot 
h4=figure('Position',[1 1 1800 1000]);
histEdgesPrnt={'0','%s/6','%s/3','%s/2','2%s/3','5%s/6','%s','-5%s/6','-2%s/3','-%s/2','-%s/3','-%s/6'};

cmap=colormap('jet');
cmap=cmap(129:end,:);
cmap=cell2mat(arrayfun(@(x) ...
    reshape([repmat(cmap(:,x),[1 2])+[zeros(size(cmap,1),1) [diff(cmap(:,x)); diff(cmap(end-1:end,x))]/2]]',[], 1),1:3,'uni',0));
cmap=cmap.*(1 - ((size(cmap,1)-1):-1:0)/(6*size(cmap,1)))';
% cmap=cmap(1:end-3,:);
cmap(cmap<0)=0;
cmap(cmap>1)=1;

% histEdgesPrnt={'0','%s/12','%s/6','%s/4','%s/3','5%s/12','%s/2',...
%     '7%s/12','2%s/3','3%s/4','5%s/6','11%s/12','%s',...
%     '-11%s/12','-5%s/6','-3%s/4','-2%s/3','-7%s/12','-%s/2',...
%     '-5%s/12','-%s/3','-%s/4','-%s/6','-%s/12'};

histEdgesPrnt=cellfun(@(x) sprintf(x,'\pi'),histEdgesPrnt,'uni',0);
% histEdges=linspace(-1, 1, 12)*pi;
histEdges=linspace(-1, 1, 25)*pi;
%     histClr='c';
%     histClr=[0 0.447 0.741];
    sigPrCplx=cos(sigPrMean)+1i*sin(sigPrMean);
    prPLV=abs(mean(sigPrCplx));
    
    plvTmp=prPLV/0.6;
    plvTmp=size(cmap,1)*plvTmp + 1;
    plvTmp(plvTmp>size(cmap,1))=size(cmap,1);
    plvTmp(plvTmp<1)=1;
    plvTmp=round(plvTmp);
    histClr=cmap(plvTmp,:);
    
    ph=polarhistogram(sigPrMean,histEdges,'FaceColor',histClr,'FaceAlpha',1);
    ph.EdgeAlpha = 1;
    ph.LineWidth = 1.5;
    set(gca,'ThetaTickLabel', histEdgesPrnt);
    prettifyPlot(gcf,gca)
    %set(gcf,'Color',[0 0 0])
    set(gca,'FontSize',32,'LineWidth',8)
    
    ax=gca;
    tmp=ax.RTickLabel;
    ax.RTickLabel=repmat({''},size(tmp));
    ax.RTick=[];
    gridSat=0;
    gridClr=[1 1 1]*gridSat;
    ax.RColor=gridClr;
    ax.ThetaColor=gridClr;
    ax.GridColor=gridClr;
    ax.MinorGridColor=gridClr;
    ax.GridAlpha=0.5;
    ax.LineWidth=0.25;
    ax.RAxis.Limits=[0 max(ph.Values)];
    
    title(sprintf('per pair (>20 corip and sig.), mean = %.2f, PLV=%.2f, n=%d',sigPrMeanMean,prPLV,numel(sigPrMean)))
    
    saveboth(h4,sprintf('%s/pairMeanHist.pdf',outDir))

h5=figure('Position',[1 1 1800 1000]);

allCorip=vertcat(prRipsAll{:});
allCoripMean=circ_mean(allCorip);
%     histClr='c';

    allCoripCplx=cos(allCorip)+1i*sin(allCorip);
    allPLV=abs(mean(allCoripCplx));
    
    plvTmp=allPLV/0.6;
    plvTmp=size(cmap,1)*plvTmp + 1;
    plvTmp(plvTmp>size(cmap,1))=size(cmap,1);
    plvTmp(plvTmp<1)=1;
    plvTmp=round(plvTmp);
    histClr=cmap(plvTmp,:);
    
%     histClr=[0 0.447 0.741];
        ph=polarhistogram(allCorip,histEdges,'FaceColor',histClr,'FaceAlpha',1);
    ph.EdgeAlpha = 1;
    ph.LineWidth = 1.5;
    set(gca,'ThetaTickLabel', histEdgesPrnt);
    prettifyPlot(gcf,gca)
    %set(gcf,'Color',[0 0 0])
    set(gca,'FontSize',32,'LineWidth',8)
    
    ax=gca;
    tmp=ax.RTickLabel;
    ax.RTickLabel=repmat({''},size(tmp));
    ax.RTick=[];
    gridSat=0;
    gridClr=[1 1 1]*gridSat;
    ax.RColor=gridClr;
    ax.ThetaColor=gridClr;
    ax.GridColor=gridClr;
    ax.MinorGridColor=gridClr;
    ax.GridAlpha=0.5;
    ax.LineWidth=0.25;
    ax.RAxis.Limits=[0 max(ph.Values)];
    
    
    title(sprintf('all corip, mean = %.2f, PLV=%.2f, n=%d',allCoripMean, allPLV, numel(allCorip)))
    
    saveboth(h5,sprintf('%s/coripHistAll.pdf',outDir))
    
%% define region pairs and map channel pairs to region pairs


regAxis={'MedOccip','Occipital','Ling','InfOccip','Fus','IT','MTG','pSTG','STG','aSTG','Angular','Supramarginal','Parahip','AT','OFC',...
    'IFG','MFG','SFG','SPL','Subcentral','Precentral','Postcentral'};
% regShrt={'MO','Oc','VO','Ling','Fus','IT','MT','ST','Ang','SM','PH','TP','OF','IF','SF','SP','SC','Pre','Pos'};
hemi='lh';
load(sprintf('%s/regionXYpos_%s_tiltedLateral_%s.mat',brainImgDir,hemi,parcType),'regLabs','regLabsShrt','locs','labOffset')
regShrt=regLabsShrt;


% regionSet=[1 2 5:8 11:14 16 17 19:21 23:28];
% regShrt=regShrt(regionSet);
allRegs=regLabs;
nReg=numel(allRegs);
nRegPr= nReg*(nReg-1)/2;
regPairs=cell(nRegPr,1); %region pair names
regPairShrt=cell(nRegPr,1);
regPairMap=nan(nRegPr,2);%map from region pairs to regions
lastIdx=0;


regMap = cell2mat(cellfun(@(x) find(strcmp(x,regAxis)),allRegs,'uni',0));
flipPlv=false(nRegPr,1);
for n=1:nReg
    for t=n+1:nReg
        regPairs{n+t+lastIdx-2}=[allRegs{n} ' - ' allRegs{t}];
        regPairShrt{n+t+lastIdx-2}=[regShrt{n} '-' regShrt{t}];
        regPairMap(n+t+lastIdx-2,:)=[n t];
        if regMap(t)>regMap(n)
            flipPlv(n+t+lastIdx-2) = true;
        end
    end
    lastIdx=lastIdx+t-n-2;
end

%set up indices for merging regions in second column into regions in first column
mergeRegs = {};%{'Occipital','InfOccip'};
mergeRegIdx = nan(size(mergeRegs));
for m=1:numel(mergeRegs)
    mergeRegIdx(m) = find(strcmp(allRegs,mergeRegs{m}));
end

reg=cell(2,1);
regIdx=nan(2,1);
for s=1:nS
    iReg{s}=nan(numel(nPr(s)),1);
    flipPr{s}=false(numel(nPr(s)),1);
    for p=1:numel(prRegs{s})
        dash=regexp(prRegs{s}{p},'-');%dash=dash(2);
        reg{1}=prRegs{s}{p}(1:dash-2);
        reg{2}=prRegs{s}{p}(dash+2:end);
%         reg{1}=prRegs{s}{p};
%         reg{2}=prRegs{s}{p};
%         
        %merge regions
        for r=1:2
            if ~any(strcmp(reg{r},allRegs));continue;end
            
            regIdx(r)=find(strcmp(reg{r},allRegs));
            if ~isempty(mergeRegs)
                match = mergeRegIdx(:,2)==regIdx(r);
                if any(match)
                    regIdx(r)=mergeRegIdx(match,1);
                end
            end
        end
        
        if isempty(regIdx(1))||isempty(regIdx(2));continue;end
        if any(all([regPairMap(:,1)==regIdx(1) regPairMap(:,2)==regIdx(2)],2))
            iReg{s}(p)=find(all([regPairMap(:,1)==regIdx(1) regPairMap(:,2)==regIdx(2)],2));
        elseif any(all([regPairMap(:,1)==regIdx(2) regPairMap(:,2)==regIdx(1)],2))
            iReg{s}(p)=find(all([regPairMap(:,1)==regIdx(2) regPairMap(:,2)==regIdx(1)],2));
            flipPr{s}(p)=true;
        elseif regIdx(1)==regIdx(2) %within region
            iReg{s}(p)=0;
            flipPr{s}(p)=false;
        else
%             error('channel pair not matched to region pair');
        end
    end
end


%%

regRips=struct();
for s=1:nS
    subj=subjects{s};
    regRips.(subj)=cell(nRegPr,1);
    for r=1:nRegPr
        
        regChanPrs=find(iReg{s}==r);
        for p=1:numel(regChanPrs)
            chanIdx = prMap{s}(regChanPrs(p),:);
%             if flipPr{s}(p);chanIdx=flip(chanIdx);end

%             chRips=find(all(coripChans{s}(:,chanIdx),2));
%             distTmp=nan(numel(chRips),1);
%             for c = 1:numel(chRips)
%                 distTmp(c) = coripPhDist{s}(chRips(c),chanIdx(1),chanIdx(2));
%             end
            
            distTmp = coripPhDist{s}(:,chanIdx(1),chanIdx(2));
            
            
%             distTmp = coripPhDist{s}(all(coripChans{s}(:,chanIdx),2),chanIdx(1),chanIdx(2));
%             distTmp = coripPhDist{s}(:,chanIdx(1),chanIdx(2));
            if flipPr{s}(p); distTmp = -distTmp; end
            distTmp(isnan(distTmp))=[];
            
            regRips.(subj){r}=[regRips.(subj){r}; distTmp];
        end
    end
end

regRipsAll=repmat({[]},[nRegPr 1]);
regPrPLV = nan(nRegPr,1);
regPh = nan(nRegPr,1);

for r=1:nRegPr
    regRipsAll{r}=[];
    for s=1:nS
        regRipsAll{r} = [regRipsAll{r}; regRips.(subjects{s}){r}];
    end
%     regRipsAll{r}(isnan(regRipsAll{r}))=[]; %FIGURE OUT WHERE THESE NANS CAME FROM
    % there should be no nans
    if isempty(regRipsAll{r});continue;end
    regPrPLV(r) = abs(mean(cos(regRipsAll{r}) + 1i*sin(regRipsAll{r})));
    regPh(r) = circ_mean(regRipsAll{r});
    
end



%% test significance of nonuniformity and deviation from zero phase
    
pRegPLV=nan(nRegPr,1);
pRegDir=nan(nRegPr,1);
% nZero=nan(nRegPr,1);
% nAll=nan(nRegPr,1);

for r=1:nRegPr
    if isempty(regRipsAll{r})
        continue
    end
    pRegPLV(r)=circ_rtest(regRipsAll{r});
    isZero=regRipsAll{r}<pi/8 & regRipsAll{r}>-pi/8;
    pRegDir(r)=myBinomTest(sum(isZero),numel(isZero),1/8,'one');
%     pRegDir(r)=circ_medtest(regRipsAll{r},0);
    
end
pRegPLV(~isnan(pRegPLV))=mafdr(pRegPLV(~isnan(pRegPLV)),'BHFDR',true);
pRegDir(pRegPLV>0.05)=nan;
pRegDir(~isnan(pRegDir)) = mafdr(pRegDir(~isnan(pRegDir)),'BHFDR',true);

% regPh(pRegDir>0.05|isnan(pRegDir))=0;

%% transform to matrices

nRipRegPr=cellfun(@(x) numel(x),regRipsAll);

regPlvMat = nan(nReg,nReg);
regPhMat = nan(nReg,nReg);
flipPlvMat = zeros(nReg,nReg);
nRipMat=zeros(nReg,nReg);
pRegPLVmat = zeros(nReg,nReg);
pRegDirMat = nan(nReg,nReg);
regNameMat = cell(nReg,nReg); %sanity check
for r=1:nRegPr
    idx=regPairMap(r,:);
    regPlvMat(idx(1),idx(2)) = regPrPLV(r);
    regPhMat(idx(1),idx(2)) = regPh(r);
    flipPlvMat(idx(1),idx(2)) = double(flipPlv(r));
    nRipMat(idx(1),idx(2)) = nRipRegPr(r);
    pRegPLVmat(idx(1),idx(2)) = pRegPLV(r);
    pRegDirMat(idx(1),idx(2)) = pRegDir(r);
    regNameMat(idx(1),idx(2)) = regPairs(r);
end
regPlvMat(isnan(regPlvMat))=0;
regPlvMat = mirrormat(regPlvMat);

flipPlvMat=mirrormat(flipPlvMat);
flipPlvMat=logical(flipPlvMat);

nRipMat=mirrormat(nRipMat);

regPhMat(isnan(regPhMat))=0;
regPhMat=mirrormat(regPhMat);
regPhMat(regPhMat==0)=nan;

regPhMatOrd = regPhMat.*flipPlvMat; %positive means that earlier/more sensory/less motor region LEADS

pRegPLVmat=mirrormat(pRegPLVmat);
pRegPLVmat(1:nReg+1:end)=nan;

nanMsk=~isnan(pRegDirMat);
pRegDirMat(~nanMsk)=0;
nanMsk=mirrormat(nanMsk);
pRegDirMat=mirrormat(pRegDirMat);
pRegDirMat(~nanMsk)=nan;

tmp=regNameMat';
regNameMat(tril(true(nReg)))=tmp(tril(true(nReg)));

close all
%%
ordToReg=[10 6 12 7 3 5 9  22 18 21 1 13 19 2 20 11 4 8 16 17 15 14];
regToOrd=arrayfun(@(x) find(ordToReg==x),1:numel(ordToReg));
regPrMapOrd(:,1)=regToOrd(regPairMap(:,1));
regPrMapOrd(:,2)=regToOrd(regPairMap(:,2));
regPrMapOrd=sort(regPrMapOrd,2);
numTab=zeros(nReg,nReg);
PLVtab=zeros(nReg,nReg);
sigTab=zeros(nReg,nReg);
for r=1:nRegPr
    
    regPrTmp=regPrMapOrd(r,:);
    
%     numTab{regPrTmp(1),regPrTmp(2)}=sprintf('%d, %.2f, %.2f',numel(regRipsAll{r}),pRegPLV(r),pRegDir(r));
    numTab(regPrTmp(1),regPrTmp(2))=numel(regRipsAll{r});
    PLVtab(regPrTmp(1),regPrTmp(2))=regPrPLV(r);
    if pRegPLV(r)>alpha
        sigVal=0;
    elseif pRegDir(r)>=alpha
        sigVal=1;
    else
        sigVal=2;
    end
    sigTab(regPrTmp(1),regPrTmp(2))=sigVal;
    
end
numTab=mirrormat(numTab);
PLVtab=mirrormat(PLVtab);
sigTab=mirrormat(sigTab);
numTabPrnt=[[{''} regShrt(ordToReg)];[regShrt(ordToReg)' num2cell(numTab)]];
PLVtabPrnt=[[{''} regShrt(ordToReg)];[regShrt(ordToReg)' num2cell(PLVtab)]];
sigTabPrnt=[[{''} regShrt(ordToReg)];[regShrt(ordToReg)' num2cell(sigTab)]];
writecell(numTabPrnt,sprintf('%s/regNcoripTab.xls',outDir));
writecell(PLVtabPrnt,sprintf('%s/regPLVtab.xls',outDir));
writecell(sigTabPrnt,sprintf('%s/regSigPhTab.xls',outDir));


%% plot

histEdgesPrnt={'0','%s/6','%s/3','%s/2','2%s/3','5%s/6','%s','-5%s/6','-2%s/3','-%s/2','-%s/3','-%s/6'};
histEdgesPrnt=cellfun(@(x) sprintf(x,'\pi'),histEdgesPrnt,'uni',0);
if plotPerReg
histEdges=linspace(-1, 1, 13)*pi;
% histEdges=linspace(-5/6, 1, 12)*pi;
% 
regShrt{ordToReg(5)}(1)='M'
regShrt{ordToReg(6)}(1)='L'
hTmp=figure('Position',[1 1 400 200]);
cmap=colormap('jet');
cmap=cmap(129:end,:);
cmap=cell2mat(arrayfun(@(x) ...
    reshape([repmat(cmap(:,x),[1 2])+[zeros(size(cmap,1),1) [diff(cmap(:,x)); diff(cmap(end-1:end,x))]/2]]',[], 1),1:3,'uni',0));
cmap=cmap.*(1 - ((size(cmap,1)-1):-1:0)/(6*size(cmap,1)))';
% cmap=cmap(1:end-3,:);
cmap(cmap<0)=0;
cmap(cmap>1)=1;
colormap(cmap);
cbar=colorbar('southoutside');
cbar.Ticks=[0:0.1:0.6];
caxis([0 0.6])
cbar.Label.String='PLV';
cbar.FontSize=12;
cbar.Label.FontSize=14;

drawnow
saveGraphic(hTmp,sprintf('%s/regHistClrBar',outDir),[],600)
close(hTmp)
% h1=figure('Position',[300 300 1000 700]);
h1=figure('Position',[1 39 910 922]);
set(h1,'Units','inches')
set(h1,'PaperPosition',get(h1,'Position'))
drawnow
set(gcf,'Color',[1 1 1])
colormap(cmap)

% ordToReg=[10 6 12 7 21 3 5 9 18 1 13 19 2 20 11 4 8 16 17 15 14];
% regToOrd=arrayfun(@(x) find(ordToReg==x),1:numel(ordToReg));
% regPrMapOrd(:,1)=regToOrd(regPairMap(:,1));
% regPrMapOrd(:,2)=regToOrd(regPairMap(:,2));
% regPrMapOrd=sort(regPrMapOrd,2);


    

gap=[0.0025 0];
marg_h=[0.02 0.02];
marg_w=[0.02 0.02];
for n=1:nReg
    figure(h1)
            subtightplot(nReg+1,nReg+1,(nReg+1)*(n+1),gap,marg_h,marg_w)
        text(0.5, 0.5, regShrt{ordToReg(n)},'FontSize',14,'HorizontalAlignment','center')
        ax=gca;
        set(ax,'XColor',[1 1 1],'YColor',[1 1 1])
%         xticklabels([])
        axis off

            subtightplot(nReg+1,nReg+1,n,gap,marg_h,marg_w)
        text(0.5, 0.5, regShrt{ordToReg(n)},'FontSize',14,'HorizontalAlignment','center')
        ax=gca;
        set(ax,'XColor',[1 1 1],'YColor',[1 1 1])
%         xticklabels([])
        axis off
        
        
%     if n<nReg
%         subtightplot(nReg,nReg,nReg*(n+1),gap,marg_h,marg_w)
%         text(0.5, 0.5, regShrt{ordToReg(n)})
%         ax=gca;
%         set(ax,'XColor',[1 1 1],'YColor',[1 1 1])
% %         xticklabels([])
%         axis off
%     end
%     if n>1
%         subtightplot(nReg,nReg,n-1,gap,marg_h,marg_w)
%         text(0.5, 0.5, regShrt{ordToReg(n)})
%         ax=gca;
%         set(ax,'XColor',[1 1 1],'YColor',[1 1 1])
% %         xticklabels([])
%         axis off
%     end
end
drawnow
for r=1:nRegPr    
    figure(h1)
    for t=1:2
    if isempty(regRipsAll{r});continue;end
    regPrTmp=regPrMapOrd(r,:);
    
    if t==1
        
    subtightplot(nReg+1,nReg+1,regPrTmp(1)*(nReg+1) + regPrTmp(2),gap,marg_h,marg_w);
    elseif t==2
        
    subtightplot(nReg+1,nReg+1,regPrTmp(2)*(nReg+1) + regPrTmp(1),gap,marg_h,marg_w);
    end
%     subplot(nReg,nReg,(regPrTmp(1)-1)*nReg + regPrTmp(2))
    
    regPlvTmp=regPrPLV(r)/0.6;
    regPlvTmp=size(cmap,1)*regPlvTmp + 1;
    regPlvTmp(regPlvTmp>size(cmap,1))=size(cmap,1);
    regPlvTmp(regPlvTmp<1)=1;
    regPlvTmp=round(regPlvTmp);
%     histClr='c';
%     histClr=[0 0.447 0.741];
    histClr=cmap(regPlvTmp,:);
        ph=polarhistogram(regRipsAll{r},histEdges,'FaceColor',histClr,'FaceAlpha',1);%,'EdgeColor',histClr);
    ph.EdgeAlpha = 0.75;
    ph.LineWidth = 0.25;
%     set(gca,'ThetaTickLabel', histEdgesPrnt);
    set(gca,'ThetaTickLabel', repmat({''},numel(histEdgesPrnt),1));
    prettifyPlot(gcf,gca)
    %set(gcf,'Color',[0 0 0])
    set(gca,'FontSize',4);%,'LineWidth',0.5)
    pDirStr=sprintf('%.2f',pRegDir(r));
    if strcmp(pDirStr,'NaN')
        pDirStr='N/A';
    end
    ax=gca;
    tmp=ax.RTickLabel;
    ax.RTickLabel=repmat({''},size(tmp));
    ax.RTick=[];
    gridSat=0;
    gridClr=[1 1 1]*gridSat;
    ax.RColor=gridClr;
    ax.ThetaColor=gridClr;
    ax.GridColor=gridClr;
    ax.MinorGridColor=gridClr;
    ax.GridAlpha=0.5;
    ax.LineWidth=0.25;
    ax.RAxis.Limits=[0 max(ph.Values)];
%     title(sprintf('%.2f, %s',pRegPLV(r),pDirStr),'FontSize',4);%,'Color','w')
%     title(regPairShrt{r})
%     drawnow
%     numTab{regPrTmp(1),regPrTmp(2)}=sprintf('%d, %.2f, %.2f',numel(regRipsAll{r}),pRegPLV(r),pRegDir(r));
    end
end

%   set(h1,'PaperUnits','inches')
%   PP = get(h1,'paperposition');
%   PP(1:2) = 0.5;
%   set(h1,'papersize',PP(3:4)+.5);
%   set(h1,'InvertHardCopy','off') % preserve background color


saveGraphic(h1,sprintf('%s/regionPLV',outDir),[],600)%,[],300)%,[],600)
% numTabPrnt=[[{''} regShrt(ordToReg)];[regShrt(ordToReg)' numTab]];
% writecell(numTabPrnt,sprintf('%s/regPLVtab.xls',outDir));

if plotSingleRegPr
    h12=figure('Position',[1 1 1800 1000]);
    for r=1:nRegPr    
    if isempty(regRipsAll{r});continue;end
%     histClr='c';
    histClr=[0 0.447 0.741];
        ph=polarhistogram(regRipsAll{r},histEdges,'FaceColor',histClr,'FaceAlpha',1);
    ph.EdgeAlpha = 1;
    ph.LineWidth = 3;
%     set(gca,'ThetaTickLabel', histEdgesPrnt);
    set(gca,'ThetaTickLabel', repmat({''},numel(histEdgesPrnt),1));
    prettifyPlot(gcf,gca)
    %set(gcf,'Color',[0 0 0])
    set(gca,'FontSize',32,'LineWidth',8)
    title(sprintf('%s p_{1}=%.2f p_{2}=%.2f',regPairs{r},pRegPLV(r),pRegDir(r)));%,'Color','w')
%     drawnow
    saveGraphic(h12,sprintf('%s/perRegPr/%s',outDir,regPairs{r}))
    end
end

%% number co-R vs significance

h3=figure('Position',[1 1 1000 1000])
scatter(nRipRegPr, log10(pRegPLV))
title('Num Co-R vs.Significance')
xlabel('Num Co-R')
xticks(0:1000:max(nRipRegPr))
ylabel('log10(p)')
savepdf(h3,sprintf('%s/numCoRipVpVal.pdf',outDir))
xlim([0 3000])
ylim([-10 0])
xticks(0:250:3000)
hold on
line(xlim, [1 1]*log10(0.05), 'Color','r','LineStyle','--')
legend('region pairs','p=0.05')
saveboth(h3,sprintf('%s/numCoRipVpVal_trimmed',outDir))


%%
save(sprintf('%s/regPLV_%s.mat',fRoot,parcType),'regPrPLV','regPh','regPairs','regRipsAll')

%% PLV graph
%

load(sprintf('%s/regionXYpos_%s_tiltedLateral_%s.mat',brainImgDir,hemi,parcType),'regLabs','regLabsShrt','locs','labOffset')
figure(1)
% colormap parula
colormap jet
cmap=colormap;
close(1)
gap=[0.03 0.05];
marg_h=[0.03 0.03];
marg_w=[0.02 0.05];

%region node positions
hemi='lh';
brainImg=imread(sprintf('%s/png/allGreyBrain_%s_tiltedLateral.png',brainImgDir,hemi));
imSz=size(brainImg);
imSz(3)=[];
locs=locs + [60 98];
labOffset=labOffset + [-25 -5];
labLocs=locs+labOffset;
locs(:,2)=(imSz(2)-locs(:,2))*0.96;
labLocs(:,2)=(imSz(2)-labLocs(:,2))*0.96;
labLocs(:,1)=(labLocs(:,1)-imSz(1))*1.01 + imSz(1) + 20;
% locs= locs*0.97 + imSz(1)/2;

% labOffset=labOffset*0.97;

% locs(:,1)=locs(:,1)+40;
sclFctr=200;
locs=locs/sclFctr;
labLocs=labLocs/sclFctr;
% locs(:,2)=locs(:,2)-1.11;
% labOffset=labOffset/sclFctr;
% labLocs=locs+labOffset;

scrnScl=max([1080 1920]./imSz(1:2));
axSz=flip(imSz(1:2)*scrnScl/0.9);


rmIdx=false(numel(regLabs),1);
regLabs(rmIdx)=[];
locs(rmIdx,:)=[];

plotRegs=cellfun(@(x) any(strcmp(x,regLabs)),allRegs);
% plotRoi=printRoi(plotRegs);

gap=[0.03 0.05];
marg_h=[0.03 0.03];
marg_w=[0.02 0.05];

edgeFctr=15;
h2=figure('Position',[1 1 1080 1080]);





image((1:imSz(2))/sclFctr, (1:imSz(1))/sclFctr, brainImg)
hold on


% lgdGrph=graph(lgdGrph);



% lgdPlt=plot(lgdGrph,'MarkerSize',eps,'NodeFontSize',eps,'LineWidth',lgdWt,'EdgeLabel',maxValStr);
% xloc=lgdPlt.XData;
% lgdPlt.XData=lgdPlt.YData;
% lgdPlt.YData=xloc;

plotPlv=regPlvMat(plotRegs,plotRegs);
plotPh=regPhMatOrd(plotRegs,plotRegs);
rmEdge=plotPlv<plvThresh|pRegPLVmat>alpha|isnan(pRegPLVmat);
plotPlv(rmEdge)=0;
maxVal=max(plotPlv,[],'all');
plotPh(rmEdge)=0;
% plotPh(plotPh==0)=[];
plotPh=plotPh(tril(true(numel(regLabs)),-1));
edgeRef=plotPlv(tril(true(numel(regLabs)),-1))==0; % unplotted edges, to remove phase (color) data for
plotPh(edgeRef)=[];

%legend values
lgdMaxPos=[imSz(2)/sclFctr imSz(1)/sclFctr]/10;
lgdLocs = [0 0; lgdMaxPos(1) 0;...
    0 lgdMaxPos(2)/2; lgdMaxPos(1) lgdMaxPos(2)/2;...
    0 lgdMaxPos(2); lgdMaxPos(1) lgdMaxPos(2)];
lgdLocs=lgdLocs+lgdMaxPos/4;
lgdGrph=zeros(6,6);
lgdGrph(1,2)=1;
lgdGrph(3,4)=1;
lgdGrph(5,6)=1;
% lgdWt=linspace(0.1,maxVal,3);
lgdWt=[0.2 0.5 0.8];
lgdGrph(lgdGrph~=0)=lgdWt;
lgdGrph=mirrormat(lgdGrph);
lgdStr=arrayfun(@(x) num2str(x,'%.2f'),lgdWt,'uni',0);

edgeStr = [repmat({''},sum(~~plotPlv,'all')/2,1); lgdStr'];
plotPlv=diagConcat(plotPlv,lgdGrph);
plotPh=[plotPh; zeros(sum(~~lgdGrph,'all')/2,1)];
    printRoi=regShrt';%(regionSet)';
% printRoi={'AT','Ang','Fus','IF','IT','IO','Ling','MF','MT','OF','Occ','PH','Pos','Pre','SF','SP','ST','Sub','SM'}';
plotLabs=[printRoi(~rmIdx); arrayfun(@(x) repmat(' ',1,x), 1:size(lgdGrph,1),'uni',0)'];
locs=[locs;lgdLocs];


% frequency coloring of nodes

load(sprintf('/home/jgarret/TaskAnalysis/lang/FW/regRipFreqs_%s.mat',parcType),'regFreqMeans','regFreqLabs')
regionSet=cellfun(@(x) find(strcmp(x,regFreqLabs)),regLabs);
plotFreqMeans=regFreqMeans(regionSet,1);
plotNodeClrs = (plotFreqMeans-90)*(127/5.5) + 128;

colormap jet
cmap=colormap;
nodeClrs=cmap(round(plotNodeClrs),:);


G=graph(plotPlv,plotLabs,'omitselfloops');

Gplot=plot(G, 'NodeColor', [nodeClrs; repmat([0.1 0.3 1], [6 1])],...
    'EdgeColor',repmat([0.1 0.3 1], numel(G.Edges.Weight),1),'EdgeAlpha',0.7,'NodeFontSize',15,'MarkerSize',20,...
    'NodeLabelColor',[1 1 1],'LineWidth',G.Edges.Weight * edgeFctr, 'EdgeLabel', edgeStr,'EdgeLabelColor',[1 1 1]);
% clrNorm=max(abs(plotPh));
clrNorm=pi;
% plotClr = floor(abs(plotPh)*(255/clrNorm))+1;
plotClr = floor(plotPh*(128/clrNorm))+129;
Gplot.EdgeColor = cmap(plotClr,:);


set(gca,'GridColor',[0 0 0])

Gplot.XData=locs(:,1);
Gplot.YData=locs(:,2);


%set line width normalized by maximum
%Gplot.LineWidth = edgeFctr*G.Edges.Weight/max(cell2mat(cellfun(@(x) max(x.Edges.Weight),G(:,1),'uni',0)),[],'all');

                   figClr=0.8*ones(1,3);
            set(gca,'color',figClr)
            set(gcf,'color',figClr)
            set(gca,'XColor',figClr,'YColor',figClr,'TickDir','out')
            title(sprintf('Co-ripple phase-locking strength and lag, PLV %.3f-%.3f',plvThresh,max(max(plotPlv))))
            cbar=colorbar;
%             colormap parula
            ylim([0 imSz(1)/sclFctr])
            xlim([0 imSz(2)/sclFctr])
%             caxis([0 1]*clrNorm)
%             caxis([-1 1]*clrNorm)
            caxis([-1 1]*(1/2)*(1000/90))
%             cbar.Label.String='mean phase lag (radians)';
            cbar.Label.String='mean phase lag (ms)';
            cbar.Color=[1 1 1];
%             cbar.Ticks=linspace(-pi,pi,13);
            cbar.Ticks=linspace(-5,5,11);
% histEdgesPrnt={'-%s','-5%s/6','-2%s/3','-%s/2','-%s/3','-%s/6','0','%s/6','%s/3','%s/2','2%s/3','5%s/6','%s'};
% histEdgesPrnt=cellfun(@(x) sprintf(x,'\pi'),histEdgesPrnt,'uni',0);
%             cbar.TickLabels=histEdgesPrnt;
%             
            
                nodeLabs=Gplot.NodeLabel;
                replotLabs=find(cellfun(@(x) any(strcmp(nodeLabs,x)),regLabsShrt));
                Gplot.NodeLabel(1:end-6)=arrayfun(@(x) repmat(' ',x,1), 1:numel(Gplot.NodeLabel)-6,'uni',0);
                for r=1:numel(replotLabs)
                    rIdx=replotLabs(r);
                    text(labLocs(rIdx,1),labLocs(rIdx,2),regLabsShrt(rIdx),'color',[1 1 1],'FontSize',40,'FontName','Arial');
                end
            
set(h2,'InvertHardCopy','off','Color',[0 0 0])
figure(h2)
ax=gca;
set(ax,'XColor',[0 0 0],'YColor',[0 0 0])
ax.Title.Color=[1 1 1];
saveGraphic(h2,sprintf('%s/plvBrainGraph_thresh_%.2f',outDir,plvThresh));

h7=figure();
colormap jet
cbar2=colorbar;
caxis([84.5 95.5])
cbar2.Ticks=linspace(min(caxis),max(caxis),11);
set(h7,'Color',[0 0 0])
cbar2.Color=[1 1 1];

saveGraphic(h7,sprintf('%s/plvBrainGraph_freqScale',outDir));
%}

%%
%{
hLgd=figure('Position',[1 1 1080 1080]);
lgdGrph=zeros(6,6);
lgdGrph(1,2)=1;
lgdGrph(3,4)=1;
lgdGrph(5,6)=1;
incrmnt=logspace(-1, 0, 3);
%         lgdGrph=lgdGrph*sclFctr;
lgdGrph=mirrormat(lgdGrph);
lgdGrph=graph(lgdGrph);%,repmat({''},[6 1]));
%maxValEdge=max(max(G.Edges(~arrayfun(@(y) strcmp(G.Edges(y,:).EndNodes{1}, G.Edges(y,:).EndNodes{2}),1:height(G.Edges))',:).Weight),[],'all');
%maxValAll=max(max(G.Edges.Weight),[],'all');
% lgdWt=incrmnt*edgeFctr(iSet)*(maxValEdge/maxValAll);
maxVal=max(G.Edges.Weight);
lgdWt = incrmnt*edgeFctr*maxVal;
% maxValStr=arrayfun(@(x) num2str(x,'%.2E'),maxValEdge*incrmnt,'uni',0);
maxValStr=arrayfun(@(x) num2str(x,'%.2E'),maxVal*incrmnt,'uni',0);
lgdPlt=plot(lgdGrph,'MarkerSize',eps,'NodeFontSize',eps,'LineWidth',lgdWt,'EdgeLabel',maxValStr);
xloc=lgdPlt.XData;
lgdPlt.XData=lgdPlt.YData;
lgdPlt.YData=xloc;
% savepdf(hLgd,sprintf('%s/%sv%s_widLgnd.pdf',figOutDir,respLab{:}));
savepdf(hLgd,'/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/FW/phaseFigs_v3/plvGraphWidLgnd.pdf')

% hLgd2=figure('Position',[1 1 1080 1080]);
% lgdGrph2=zeros(3,3);
% lgdGrph2(1,1)=1;
% lgdGrph2(2,2)=1;
% lgdGrph2(3,3)=1;
% lgdGrph2=graph(lgdGrph2,arrayfun(@(x) num2str(x,'%.2E'),incrmnt*nodeSzNorm,'uni',0),'omitselfloops');
% plot(lgdGrph2,'MarkerSize',30*incrmnt)
% savepdf(hLgd2,sprintf('%s/%sv%s_szLgnd.pdf',figOutDir,respLab{:}));


%}
end

end
