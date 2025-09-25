function tStart = Movie_lumpSubjDat(subjIdcs,parcType,fuseBins,HGburst,zeroPh,stimRespIdx,fdr,fuseBinOffset)
%derp2
% addpath(genpath('.'))
if nargin==0||isempty(subjIdcs)
    subjIdcs=1
end
if nargin<2||isempty(parcType)
%     parcType='fusSplit7';
%     parcType='ST';
%     parcType='UA';
%     parcType = 'combine'
    parcType = 'RutishauserLab'
end
if nargin<3
    fuseBins=1
end
if nargin<8
    fuseBinOffset=0
end
if nargin<4
    HGburst=0;
end
% if nargin<5
    ctOpt=3;
% end
if nargin<5
    zeroPh=0;
end
if nargin<6
    stimRespIdx=[1];
end
if nargin<7
    fdr=0;
end

if HGburst
    HGapnd='_HGburst';
else
    HGapnd='';
end

if zeroPh==1
    zeroPhApnd='_zeroPh';
elseif zeroPh==2
    zeroPhApnd='_zeroPhProp';
elseif zeroPh==0
    zeroPhApnd='';
end

if ctOpt==1
    ctApnd=''; %count co-ripple centers
elseif ctOpt==2
    ctApnd='_onset' %count co-ripple starts
elseif ctOpt==3
    ctApnd='_dens' %count co-ripple durations (density)
end

if fdr
    fdrApnd='';
else
    fdrApnd='_noFDR'
end
%%

try
    parpool(16)
catch
    delete(gcp('nocreate'))
    parpool(16)
end

tStart=datetime
fRoot='/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/code/ConnectionPlots/out';
switch parcType
    case {'UA', 'combine', 'RutishauserLab'}
        allParc=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls','ST'));
    otherwise
        allParc=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls',parcType));
end
allParc=allParc(:,1);
allParc=unique(allParc);
Lparc=allParc;
allParc=cellfun(@(x) cellfun(@(y) [x y],allParc,'uni',0),{'lh-';'rh-'},'uni',0);
allParc=[allParc{1};allParc{2}];


roiIdx=[1:2 5:8 11:14 16 17 19:21 23:27];
switch parcType
    case'ST'
        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:29];
    case'fusSplit2'
        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:28];
    case 'fusSplit3'
        roiIdx=[1 2 3 5 6 7 10 11 12 14 15 17 18 20 21 22 23 24 25];
    case 'fusSplit4'
        roiIdx=[1 2 5 6 9 10 11 12 13 15 16 18 19 20 22 23 24 25 26];
    case'fusSplit6'
        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:28];
    case'fusSplit7'
        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:28];
    case'fusSplit8'
        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:28];
    case 'UA'
        roiIdx = [];
        [clusterNum] = gridCluster(10, 5);

        home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
        chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
        chLocs = dlmread(chmap_file);

%         for r = 1:96
%             roiIdx(r) = clusterNum(chMap == r);
%         end

        roiIdx = 1:max(clusterNum(:));
    case 'combine'
        [clusterNum] = gridCluster(10, 5);

        roiIdx=[1:2 5:8 11:14 16 17 19:21 23:29];
        nECoGroi = length(roiIdx);
        strRoi=allParc(roiIdx);
%         UAparc = max(roiIdx)+1:max(roiIdx)+max(clusterNum(:));
% 
        roiIdx(end+1) = max(roiIdx)+1;

%         strRoi=[strRoi' arrayfun(@(X) {num2str(X)}, UAparc)];
        strRoi=[strRoi' 'UA'];
        Lparc(end+1) = {'UA'};
    case 'RutishauserLab'
        strRoi = {'deepOFC', 'deepACC', 'deepSMA', 'deepAMY', 'deepHIP', 'deepSPE', ...
                  'supOFC', 'supACC', 'supSMA', 'supAMY', 'supHIP', 'supSPE'};
        roiIdx = 1:length(strRoi);



    
end
if strcmp(parcType, 'UA')
    strRoi=arrayfun(@(X) {num2str(X)}, roiIdx);
elseif ~any(strcmp(parcType, {'combine', 'RutishauserLab'}))
    strRoi=allParc(roiIdx);
end
nRoi=numel(roiIdx);
%% collate stimulus-evoked co-rippling
%
%stim

for sr=stimRespIdx
if sr==1
%     subjIdcs=flip([68 86 -89 101 130 133 190 226 233 -240]);
%     subjIdcs=[49 51];
%     subjIdcs=[49 51];
%     subjIdcs=flip([68 86 -89 101 130 133 190 226 -240]);
%     stimCond={'familiar','novel','all','base200','base100'};
    stimCond={'familiar','novel','confident','unsure','all','base200','base100'};

    respApnd='';
    if ~fuseBins
%         testBins=3:10;
%         testBins=2:5;
        testBins=6:30;
        binSz=100;
    else
%         testBins=[2 3 4 5];
        testBins=2:10;
        binSz=400;
        if fuseBinOffset
            testBins(end)=[];
        end
    end
elseif sr==2
%     subjIdcs=[49 51];
    stimCond={'TP'};
    respApnd='_resp';
    if ~fuseBins
        testBins=3:7;
        binSz=100;
    else
        if ~fdr
            testBins=1:5;
            if fuseBinOffset
                testBins(testBins==5)=[];
            end
        else
            testBins=2:4;
        end
        binSz=200;
    end
end
subjects={'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};
nStimCond=numel(stimCond);
nBin=40;



nTestBin=numel(testBins);

nSubj=numel(subjects);

nPairReg=zeros(nRoi,nRoi,nSubj);
nChReg=zeros(nRoi,nSubj);
nTrlSub=zeros(nStimCond,nSubj);

nTrlRip=zeros(nRoi,nStimCond); %total number trials per region and condition for ripple rate
nTrlCorip=zeros(nRoi,nRoi,nStimCond); %total number trials per region pair and condition for co-ripple rate

nSubjRip=zeros(nRoi,1); %number of subjects contributing to each region ripple rate estimation
nSubjCorip=zeros(nRoi,nRoi); %number of subjects contributing to each region pair co-ripple rate estimation

nCorip=nan(nRoi,nRoi,nBin,nStimCond,nSubj);
coripDenom=zeros(nRoi,nRoi,nBin,nStimCond,nSubj);

nRip=nan(nRoi,nBin,nStimCond,nSubj);
ripDenom=zeros(nRoi,nBin,nStimCond,nSubj);

nChAll=zeros(nSubj,1);
nBoot = 1e3;
for s=1:numel(subjects)
    subj=subjects{s}
    load(sprintf('%s/coripChanDat%s_%s%s%s%s.mat',fRoot,respApnd,subj,HGapnd,zeroPhApnd,ctApnd),'coripCond','gCh');
    load(sprintf('%s/coripDist%s_%s%s%s%s.mat',fRoot,respApnd,subj,HGapnd,zeroPhApnd,ctApnd),'coripDistSubj','tests');
    if s==1
%         if sr==1
%             nTest=4;
%         else
            nTest=size(tests,1);
%         end
        coripDist=zeros(nRoi,nRoi,nBin,nBoot,nTest,2,numel(subjects));
        ripDist=zeros(nRoi,nBin,nBoot,nTest,2,numel(subjects));
    end
    if strcmp(parcType, 'UA')
        labParc=retrieveChanParcels_MG(subj,'ST');

        labParc(91:96,:) = labParc(1:6,:);

    elseif strcmp(parcType, 'combine')
        labParc=retrieveChanParcels_MG(subj,'ST');
        labParc(end+1,1) = {'UA'};
        labParc(end,3) = {'lh'};
    elseif strcmp(parcType, 'RutishauserLab')
        filename = sprintf('%s_rippleStats.mat', subj);
        M = load(fullfile('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/ripDetections', filename));
        
        lab = M.rippleStats.chanLabels;
        nL = numel(lab);
        [alph,numr] = deal(cell(nL,1));
        labParc = cell(nL, 2);
        for iL = 1:nL
          nMsk = isstrprop(lab{iL},'digit');
          if ~nMsk(end) % if label does not end in a number
            alph{iL} = lab{iL};
            numr{iL} = NaN;
          else% if label ends in a number
            % consider everything before the number at the end the name aka alph
            % (no matter how many digits that number is)
            idx = find(nMsk,1,'first');
            alph{iL} = lab{iL}(2:idx-1);
            numbp = lab{iL}(idx:end);
            if strcmp(M.rippleStats.recordingType{iL}, 'macro')
                numr{iL} = sscanf(numbp, '%f-%f');
            else
                numr{iL} = str2double(numbp);
            end
          end

          if numr{iL}(1) > 4 && strcmp(M.rippleStats.recordingType{iL}, 'macro')
            depth = 'sup';
          else
            depth = 'deep';
          end
          lParc = sprintf('%s%s', depth,  alph{iL});
          labParc{iL,1} = lParc;
          labParc{iL,2} = lab{iL};
        end
        

    else
        labParc=retrieveChanParcels_MG(subj,parcType);
    end
    labParc=labParc(gCh,:);

    nChAll(s)=numel(gCh);
    
    chanMap=false(numel(gCh),numel(roiIdx));
    if strcmp(parcType, 'UA')
        home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
        chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
        
        chLocs = dlmread(chmap_file);
        
        [clusterNum] = gridCluster(size(chLocs,1, 2), 5);

         for r = 1:size(chanMap,1)
            cluster = clusterNum(chLocs == r);
            chanMap(r,cluster) = true;
        end
    elseif strcmp(parcType, 'combine')  
        for r=1:numel(roiIdx)
            chanMap(:,r)=strcmp(labParc(:,1),Lparc(roiIdx(r)))&strcmp(labParc(:,3),'lh');
        end
           home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
        chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
        
        
    elseif strcmp(parcType, 'RutishauserLab') 
        Lparc = strRoi;
        for r=1:numel(roiIdx)
            chanMap(:,r)=strcmp(labParc(:,1),Lparc(roiIdx(r)));
        end
        
    end
    
    for b=1:nBin
        for ts=1:nTest
            parfor bt=1:nBoot
                for d=1:2
                    coripDistChan=mirrormat(double(coripDistSubj(:,:,b,bt,ts,d)));
                    ripDistChan=coripDistChan(1:size(coripDistChan,1)+1:end); %transfer single ripples
                    ripDist(:,b,bt,ts,d,s)=ripDistChan*chanMap;
                    coripDistChan(1:size(coripDistChan,1)+1:end)=0; %remove single ripples from corip data
                    coripDist(:,:,b,bt,ts,d,s)=(coripDistChan*chanMap)'*chanMap;
                    if ~issymmetric(coripDist(:,:,b,bt,ts,d,s))
                        error('symmetry sanity check fail on simulated test distribution')
                    end
                end
            end
        end
    end
    
    nChReg(:,s)=sum(chanMap,1);
    nPairReg(:,:,s)=nChReg(:,s).*nChReg(:,s)';
    
    
    for c=1:nStimCond
        if isempty(coripCond{c});continue;end
        if any(size(coripCond{c}) == 1); continue; end
        nTrlSub(c,s)=size(coripCond{c},4);
        for b=1:nBin
            coripTrl=nan(nRoi,nRoi,nTrlSub(c,s));
            ripTrl=nan(nRoi,nTrlSub(c,s));
            for t=1:nTrlSub(c,s)
                tmpRip=coripCond{c}(:,:,b,t);
                
                ripTrl(:,t)=tmpRip(1:nChAll(s)+1:end)*chanMap;
                
                tmpRip(1:nChAll(s)+1:end)=0; %remove single ripples
                tmpRip=mirrormat(tmpRip);
                coripTrl(:,:,t) = (tmpRip*chanMap)'*chanMap;
                tmp=coripTrl(:,:,t).';
                tmp2=coripTrl(:,:,t);
                triangDiff=sum(abs(tmp(:)-tmp2(:)));
                if triangDiff>eps*10000
                    error('symmetry sanity check fail on co-ripple counts')
                else
                    tmp(triu(true(size(tmp)),1))=tmp2(triu(true(size(tmp2)),1));
                    coripTrl(:,:,t)=tmp;
                end
            end
            nCorip(:,:,b,c,s)=sum(coripTrl,3);
            nRip(:,b,c,s)=sum(ripTrl,2);
        end
        
        
        coripDenom(:,:,:,c,s)=repmat(nPairReg(:,:,s)*nTrlSub(c,s), [1 1 nBin]); %number of chances for co-ripple to occur
        ripDenom(:,:,c,s)=repmat(nChReg(:,s)*nTrlSub(c,s), [1 nBin]); %number of chances for ripple to occur
        
    end
    nTrlRip = nTrlRip + repmat(logical(nChReg(:,s)),[1 nStimCond]) .* nTrlSub(:,s)';
    nTrlCorip = nTrlCorip + repmat(logical(nPairReg(:,:,s)),[1 1 nStimCond]) .* permute(nTrlSub(:,s),[2 3 1]);
    nSubjRip = nSubjRip + logical(nChReg(:,s));
    nSubjCorip = nSubjCorip + logical(nPairReg(:,:,s));
end
coripDist=sum(coripDist,7);

if fuseBins
    nBinNew=floor((nBin-fuseBinOffset)/2);
    for b=1:nBinNew
        binIdx=[2*(b-1)+1, 2*b]+fuseBinOffset;
        nCorip(:,:,b,:,:)=sum(nCorip(:,:,binIdx,:,:),3);
        coripDenom(:,:,b,:,:)=sum(coripDenom(:,:,binIdx,:,:),3);
        
        nRip(:,b,:,:)=sum(nRip(:,binIdx,:,:),2);
        ripDenom(:,b,:,:)=sum(ripDenom(:,binIdx,:,:),2);
        
        coripDist(:,:,b,:,:,:)=sum(coripDist(:,:,binIdx,:,:,:),3);
        ripDist(:,b,:,:,:)=sum(ripDist(:,binIdx,:,:,:),2);
    end
    nBin=nBinNew;
    nCorip=nCorip(:,:,1:nBin,:,:);
    coripDenom=coripDenom(:,:,1:nBin,:,:);
    nRip=nRip(:,1:nBin,:,:);
    ripDenom=ripDenom(:,1:nBin,:,:);
    coripDist=coripDist(:,:,1:nBin,:,:,:);
    ripDist=ripDist(:,1:nBin,:,:,:);
    %     testBins(testBins==5)=[];
    
end
coripDist=sum(coripDist,7);
ripDist=sum(ripDist,6);

nChRegPool=sum(nChReg,2); %total number of channels per region
nPairRegPool=sum(nPairReg,3); %total number of channel pairs per region pair

%pool and test (co-)rippling rates
nCorip(isnan(nCorip))=0;
% nCorip=mirrormatND(nCorip);
nCoripPool=sum(nCorip,5);
coripDenomPool=sum(coripDenom,5);

coripRate=nCoripPool./coripDenomPool;

nRipPool=nansum(nRip,4);
ripDenomPool=sum(ripDenom,4);
ripRate=nRipPool./ripDenomPool; %pooled co-ripple rate

%% stats
stimTest=tests;
pCoripCompBoot=nan(nRoi,nRoi,nTestBin,nTest,2); %roi x roi x bin x test x direction (pos/neg/two-tailed)
pRipCompBoot=nan(nRoi,nTestBin,nTest,2);
pCoripRegion=nan(nRoi,nTestBin,nTest,2);
for b=1:nTestBin
    for r1=1:nRoi
        for t=1:nTest
%             if any(stimTest(t,1)==[1 6],2) && b>2 && sr==1 && ~fuseBins
%                 continue
%             end
            
%             N1=ripDenomPool(r1,testBins(b),stimTest(t,1));
%             N2=ripDenomPool(r1,testBins(b),stimTest(t,2));
% 
%             currBootSamps=squeeze(double(ripDist(r1,testBins(b),:,t,:)))./[N1 N2]; % bootstrap matrix is upper triangle
%             currBootSamps=currBootSamps(:,1)./currBootSamps(:,2);
%             measRelRate = ripRate(r1,testBins(b),stimTest(t,1)) / ripRate(r1,testBins(b),stimTest(t,2));
% 
%             if isnan(measRelRate)
%                 continue
%             end
%             if any(isnan(currBootSamps))
%                 error('there are bootstrap samples with nan rate corresponding to a non-nan measured rate')
%             end
%             pRipCompBoot(r1,b,t,1) = mean(currBootSamps>=measRelRate);
%             pRipCompBoot(r1,b,t,2) = mean(currBootSamps<=measRelRate);

%             if any(coripDenomPool(r1,testBins(b),stimTest(t,:))==0,4)
%                 pRipCompBoot(r1,b,t,:) = nan;
%             end

            
            
%             N1=sum(coripDenomPool(r1,:,testBins(b),stimTest(t,1)));
%             N2=sum(coripDenomPool(r1,:,testBins(b),stimTest(t,2)));
%             
%             currBootSamps=squeeze(double(sum(coripDist(r1,:,testBins(b),:,t,:),2)))./[N1 N2]; % bootstrap matrix is upper triangle
%             currBootSamps=currBootSamps(:,1)./currBootSamps(:,2);
%             measRelRate = sum(coripRate(r1,:,testBins(b),stimTest(t,1))) / sum(coripRate(r1,:,testBins(b),stimTest(t,2)));
%             
%             if isnan(measRelRate)
%                 continue
%             end
%             if any(isnan(currBootSamps))
%                 error('there are bootstrap samples with nan rate corresponding to a non-nan measured rate')
%             end
%             pCoripRegion(r1,b,t,1) = mean(currBootSamps>=measRelRate);
%             pCoripRegion(r1,b,t,2) = mean(currBootSamps<=measRelRate);
%             

                
            for r2=1:r1
                N1=coripDenomPool(r1,r2,testBins(b),stimTest(t,1));
                N2=coripDenomPool(r1,r2,testBins(b),stimTest(t,2));
                
                currBootSamps=squeeze(double(coripDist(r1,r2,testBins(b),:,t,:)))./[N1 N2]; % bootstrap matrix is upper triangle
                currBootSamps=currBootSamps(:,1)./currBootSamps(:,2);
                measRelRate = coripRate(r1,r2,testBins(b),stimTest(t,1)) / coripRate(r1,r2,testBins(b),stimTest(t,2));
                
                if isnan(measRelRate)
                    continue
                end
                if any(isnan(currBootSamps))
                    error('there are bootstrap samples with nan rate corresponding to a non-nan measured rate')
                end
                pCoripCompBoot(r1,r2,b,t,1) = mean(currBootSamps>=measRelRate);
                pCoripCompBoot(r1,r2,b,t,2) = mean(currBootSamps<=measRelRate);

                if r1==r2
                    pCoripCompBoot(r1,r2,b,t,:)=nan;
                end
                
                if any(coripDenomPool(r1,r2,testBins(b),stimTest(t,:))==0,4)
                    pCoripCompBoot(r1,r2,b,t,:) = nan;
                end
            end
        end
    end
end

if fdr
for t=1:nTest
    for d=1:2
        tmp=pCoripCompBoot(:,:,:,t,d);
        %         tmp(~isnan(tmp))=1-fdr_bky(tmp(~isnan(tmp)),alpha);
        tmp(~isnan(tmp))=mafdr(tmp(~isnan(tmp)),'BHFDR',true);
        pCoripCompBoot(:,:,:,t,d)=tmp;
        
        tmp2=pRipCompBoot(:,:,t,d);
        tmp2=reshape(mafdr(tmp2(:),'BHFDR',true),size(tmp2));
        pRipCompBoot(:,:,t,d)=tmp2;
        
        for b=1:nTestBin
            tmp2=pCoripCompBoot(:,:,b,t,d);
            tmp2(isnan(tmp2))=1;
            tmp2 = 1 - mirrormat(1 - tmp2);
            pCoripCompBoot(:,:,b,t,d)=tmp2;
        end
    end
end
end


% mirror along diagonal
% pCoripCompBoot(isnan(pCoripCompBoot))=1;
% pCoripCompBoot=1 - mirrormatND(1 - pCoripCompBoot);
%2 tail at 0.02
% pCoripCompBoot(:,:,:,:,3)=double(~any(pCoripCompBoot<=alpha,5));
% pCoripCompBoot=pCoripCompBoot(:,:,:,:,testDir);
%%
save(sprintf('%s/lumpDat_%s_%s_%s%s%s_%d_Toffset%d%s%s.mat',fRoot,'All',respApnd,parcType,HGapnd,zeroPhApnd,binSz,fuseBinOffset,ctApnd,fdrApnd),...
    'coripRate','nCorip','coripDenom',...
    'ripRate','nRip','ripDenom',...
    'pCoripCompBoot','pRipCompBoot','coripDist','ripDist',...
    'nChReg','nPairReg','nSubjCorip','tests','stimCond',...
    'strRoi','roiIdx','subjects','testBins',...
    '-v7.3','-nocompression')
end
end