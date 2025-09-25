function WN_plotCoripDat(subjIdcs,fuseBins,parcType,fuseBinOffset,zeroPh,ctOpt,alpha,testDir,resp,runFDR,subBinIdx,plotTests)
rmpath(genpath('/home/jgarret/mat/old/'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/STG/WN'))
%% plot coripDat to show modulation over subjects
%v2 separate plotting from data collation/statistics
if nargin<3||isempty(parcType)
%     parcType='fusSplit2';
%     parcType='split';
%       parcType='ST';
      parcType='UA';
      parcType = 'combine';
%     parcType='fusSplit7'
end
subjects=genSubjList(subjIdcs,'MG');
% subjects={'MG49'};
% if nargin==0||isempty(eveIdxAll)
    eveIdxAll=[1];
% end

% if nargin<2
% %     coEveIdxAll=1:3;
    coEveIdxAll=1;
% end
if nargin<2
    fuseBins=1
end
if ~fuseBins
    fuseBinOffset=0
elseif nargin<4
    fuseBinOffset=0;
end
if nargin<5
    zeroPh=0;
end
if nargin<6
    ctOpt=3;
end
if nargin<7||isempty(alpha)
    alpha=0.05;
end
if nargin<8
    testDir=1;
end
if nargin<9
    resp=0;
end
if nargin<10
    runFDR=1;
end

if ctOpt==1
    ctApnd=''; %count co-ripple centers
elseif ctOpt==2
    ctApnd='_onset' %count co-ripple starts
elseif ctOpt==3
    ctApnd='_dens' %count co-ripple durations (density)
end

if resp
    respApnd='_resp';
else
    respApnd='';
end

if zeroPh==1
    zeroPhApnd='_zeroPh';
elseif zeroPh==2
    zeroPhApnd='_zeroPhProp';
elseif zeroPh==0
    zeroPhApnd='';
    
end

if fuseBins
    binSz=200;
else
    binSz=100;
end

if ~resp
    if fuseBinOffset
%         bins=100:binSz:900;
          bins=-900:binSz:2900;
    else
%         bins=0:binSz:1000;
        bins=-1000:binSz:2000;
    end
else
%     if fuseBinOffset
%         bins=-500:binSz:500;
%     else
%         bins=-600:binSz:4000;
%     end
end
if runFDR
    fdrApnd='_noFDR';
else
    fdrApnd='';
end

% runFDR=0;

%%

eveTypes={'', '_ripHG','_ripNoHG';...
    '_HGburst','_HGrip','_HGnoRip'};
eveNames=reshape([{'co-ripple'}  cellfun(@(x) ['co-' x(2:end)],eveTypes(2:end),'uni',0)],size(eveTypes));
eveNames{2,1}='co-HG-burst';
for eveIdx=eveIdxAll
    for coEveIdx=coEveIdxAll
        
eveApnd=eveTypes{eveIdx,coEveIdx};
eveName=eveNames{eveIdx,coEveIdx}

% clearvars
% datDir='/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/coripDat_v2/';
datDir='/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/STG/WN/out';
load(sprintf('%s/lumpDat_%s_%s_%s%s%s_%d_Toffset%d%s%s.mat',datDir,[subjects{:}],respApnd,parcType,eveApnd,zeroPhApnd,binSz,fuseBinOffset,ctApnd,fdrApnd),...
    'coripRate',...'nCorip','coripDenom',...
    'ripRate',...'nRip','ripDenom',...
    'pCoripCompBoot',...'coripDist',...
    'nChReg',...'nPairReg','nSubjCorip',
    'tests',...
    'strRoi','roiIdx',...
    'subjects','testBins','stimCond'...
    );
if ~resp
%     plotTests= [1 3 4 5 8];
%     plotTests=3:6;
%      plotTests=1:6;
%     plotTests=5:6;
%     plotTests=7;
%     plotTests=[4 6];
%     plotTests=1:size(tests,1);
% plotTests=7;
%     stimCond={'match','mismatch','word','noise','base200','base100'};
    stimCond={'C','I','W','N','B1','B2'};
else
%     plotTests=[1 size(tests,1)];
%     plotTests=1:size(tests,1);
%     plotTests=[1 3];
%     plotTests=3;
end
if ~exist('plotTests','var')
    plotTests=1:size(tests,1);
end
if ctOpt==3
    ripRate=ripRate/100;
    coripRate=coripRate/100;
end
pCoripCompBoot=pCoripCompBoot(:,:,:,:,testDir);
close all
saveFigs=1;
plotPerReg=0;

if testDir==3
    apndTest='_2tail';
elseif testDir==2
    apndTest='_reverseContrast';
else
    apndTest='';
end

if nargin<11||runFDR==0
%     subBinIdx=1:numel(testBins);
    if fuseBins
%         subBinIdx=3:numel(testBins);
          subBinIdx=3:8;
    else
        subBinIdx=6:10;
    end
        
end


    
outDir=sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/UtahArray/STG/WN/out/figs/conXplot_%s_%s_%s%s_%s_p%d_bin%d_offset%d%s%s_%d-%d/',...
    [subjects{:}],respApnd,parcType,zeroPhApnd,eveApnd,alpha*100,binSz,fuseBinOffset,ctApnd,fdrApnd,subBinIdx(1),subBinIdx(end));

hemi='lh';
% apndAvgRef='_avgRef';


if strcmp(parcType, 'UA') 
    allParc=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls','ST'));

    units = LoadSpikeTimes('MG49','CellExplorer', 'WN');
    unitChannels = cellfun(@(X) any(strcmp(X, {'pyr', 'int', 'mult'})), units(:,3));
    unitChannels = unique(cell2mat(units(:,1)));
elseif  strcmp(parcType, 'combine')
        allParc=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls','ST'));

else
    allParc=readcell(sprintf('/home/jgarret/TaskAnalysis/destrieuxParcMap_%s.xls',parcType));
end


allParc=allParc(:,1);
allParc=unique(allParc);
allParc=cellfun(@(x) cellfun(@(y) [x y],allParc,'uni',0),{'lh-';'rh-'},'uni',0);
allParc=[allParc{1};allParc{2}];



nSubj=numel(subjects);

%% set parameters for co-ripple plotting

if ctOpt==3
    densFctr=2;
else
    densFctr=1;
end
% pThresh=pThresh*[1 1 1 1 1 1 1 1 1 1 1 1]; %[0.001 0.001];
% edgeFctr=20;%*[1 1 1 1 1 1 1 1 1 1 1 1 1];%[20 20];
edgeFctr=15*densFctr;
nodeFctr=40;
clrRange=1.5;%[1 1 1 1 1 1 1 1 1 1 1 1 1 1];

if strcmp(parcType, 'UA')
    load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/regShrt_%s.mat','ST'),'regShrt')
    printRoi = strRoi;
%     nChReg = false(1,length(strRoi));
%     nChReg(unitChannels) = true;
elseif  strcmp(parcType, 'combine')
    load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/regShrt_%s.mat','ST'),'regShrt')
    regShrt(end+1) = {'UA'};
    printRoi=regShrt(roiIdx);
    printRoi=cellfun(@(x) [upper(x(1)) x(2:end)],printRoi,'uni',0);
else
    load(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/regShrt_%s.mat',parcType),'regShrt')
    printRoi=regShrt(roiIdx);
    printRoi=cellfun(@(x) [upper(x(1)) x(2:end)],printRoi,'uni',0);
end
    
    

nChRegPool=sum(nChReg,2);
stimTest=tests;
nTest=size(stimTest,1);
testBins=testBins(subBinIdx);
pCoripCompBoot=pCoripCompBoot(:,:,subBinIdx,:);
nTestBin=numel(testBins);

if runFDR
    for t=1:nTest
        tmp=pCoripCompBoot(:,:,:,t);
%                 tmp(~isnan(tmp))=1-fdr_bky(tmp(~isnan(tmp)),alpha);
        STidx=bounds2mask([18 21 22],22);
%         tmp(~STidx,~STidx,:)=nan;
        tmp(~isnan(tmp))=mafdr(tmp(~isnan(tmp)),'BHFDR',true);
        pCoripCompBoot(:,:,:,t)=tmp;
        for b=1:nTestBin
            tmp2=pCoripCompBoot(:,:,b,t);
            tmp2(isnan(tmp2))=1;
            tmp2 = 1 - mirrormat(1 - tmp2);
            pCoripCompBoot(:,:,b,t)=tmp2;
        end
    end
end

%% stim condition contrast
%

%region node positions
brainImgDir='/space/seh8/5/halgdev/projects/jgarret/FW_data/brainView/';
if strcmp(parcType, 'UA') || strcmp(parcType, 'combine')
    load(sprintf('%s/regionXYpos_%s_tiltedLateral_%s.mat',brainImgDir,hemi,'ST'),'regLabs','regLabsShrt','locs','labOffset')

else
    load(sprintf('%s/regionXYpos_%s_tiltedLateral_%s.mat',brainImgDir,hemi,parcType),'regLabs','regLabsShrt','locs','labOffset')

end
regLabsShrt=cellfun(@(x) [upper(x(1)) x(2:end)],regLabsShrt,'uni',0);
brainImg=imread(sprintf('%s/png/allGreyBrain_%s_tiltedLateral.png',brainImgDir,hemi));
imSz=size(brainImg);
imSz(3)=[];
locs=locs + [45 104];
labOffset=labOffset + [-25 -5];
labLocs=locs+labOffset;
locs(:,2)=(imSz(2)-locs(:,2))*0.96;
labLocs(:,2)=(imSz(2)-labLocs(:,2))*0.96;
labLocs(:,1)=(labLocs(:,1)-imSz(1))*1.01 + imSz(1) + 20;

if strcmp(parcType, 'combine')
    regLabsShrt(end+1) = {'UA'};
end

if strcmp(parcType, 'UA')
    home_dir = '/space/seh8/2/mdeh1-3/halgdev/projects/cgonzalez/Units';
    chmap_file = sprintf('%s/mg49_and_mg29_neuroport_channel_map.csv',home_dir);
    
    chMap = dlmread(chmap_file);
    xmin = min(locs(:,1)); xmax = max(locs(:,1)); 
    ymin = min(locs(:,2)); ymax = max(locs(:,2)); 

    xnorm = linspace(xmin,xmax,10);
    ynorm = linspace(ymin,ymax,10);
    
    [clusterNum] = gridCluster(10, 2);

    UAlocs = nan(length(roiIdx),2);

    for r = 1:size(UAlocs,1)
        [xx, yy] = find(clusterNum == r);
        UAlocs(r,:) = [mean(xnorm(xx)) mean(ynorm(yy))];
    end

    locs = UAlocs;

end

% locs= locs*0.97 + imSz(1)/2;

% labOffset=labOffset*0.97;

% locs(:,1)=locs(:,1)+40;
sclFctr=200;
locs=locs/sclFctr;
locs(end+1,:)=[3, 4.5];

labLocs=labLocs/sclFctr;
labLocs(end+1,:)=[3, 4.5];

% labLocs(end+1,:)=labLocs(end,:) - 1;

% locs(:,2)=locs(:,2)-1.11;
% labOffset=labOffset/sclFctr;
% labLocs=locs+labOffset;

scrnScl=max([1080 1920]./imSz(1:2));
axSz=flip(imSz(1:2)*scrnScl/0.9);

%legend position
lgdMaxPos=[imSz(2)/sclFctr imSz(1)/sclFctr]/8;
lgdLocs = [0 0; lgdMaxPos(1) 0;...
    0 lgdMaxPos(2)/2; lgdMaxPos(1) lgdMaxPos(2)/2;...
    0 lgdMaxPos(2); lgdMaxPos(1) lgdMaxPos(2)];
lgdLocs=lgdLocs + lgdMaxPos/3;
locs=[locs;lgdLocs];

% region labels and placeholders for legend
% plotRegs=cellfun(@(x) any(strcmp(x,regLabs)),strRoi);
plotRegs=true(numel(strRoi),1);
plotRoi=printRoi(plotRegs)';
plotRoi=[plotRoi; arrayfun(@(x) repmat(' ',1,x), 1:size(lgdLocs,1),'uni',0)'];

% comparison parameters
nChThresh=1;
close all
% stimCondSets={'word','CS';'target','CS';'target','word'};
stimCondSets=stimCond(stimTest);
% stimCondSets={'target','pretaskbaseline'};
stimCondSetIdx=cellfun(@(x) find(strcmp(x,stimCond)),stimCondSets);
%
figure(1)
colormap jet
cmap=colormap;
close(1)
gap=[0.03 0.05];
marg_h=[0.03 0.03];
marg_w=[0.02 0.05];
% pCoripCompBoot=pCoripCompBoot(:,:,1:3,:);
timeLab=arrayfun(@(x) sprintf('%d-%d ms',bins(x),bins(x+1)),testBins,'uni',0);
rippleTestLgd=cell(nTest,1);
coripTestLgd=cell(nTest,1);
for iSet=plotTests
    
    close all
    stimLab=stimCondSets(iSet,:)
    stimCondIdx=stimCondSetIdx(iSet,:);
    
    G=cell(nTestBin,2);
    Gplot=G;
    emptyNodes=G;
%     
    rejRegs=nChRegPool<nChThresh|~plotRegs|~any(pCoripCompBoot(:,:,:,iSet)<=alpha,[2 3]);
    if all(rejRegs)
%         continue
    end
%     rejRegs=nChRegPool<nChThresh|~plotRegs|all(pCoripComp(:,:,:,iSet)>pThresh(iSet),[2 3]);
    tRip=cell(nTestBin,2);
    iNode=cell(nTestBin,2);
    edgeStr=cell(nTestBin,2);
    h=cell(nTestBin,2);
    
    
    
        %legend values
    if eveIdx==1&&coEveIdx==1&&iSet==2||isempty(rippleTestLgd{iSet}) %carry over ripple scaling to any other event type
    
        %compute corip reference rates
        refRates=coripRate(:,:,testBins,stimCondIdx(1));%(:,:,:,stimCondIdx(1));
        refRates(pCoripCompBoot(:,:,:,iSet)>alpha)=nan;
        refRates=refRates(~rejRegs,~rejRegs,:);%,testBins);
        refRates(refRates==0)=nan;
        %legend co-ripple increments
%         coripRateLgdVals=linspace(min(refRates,[],'all'),max(refRates,[],'all'),3);
%         coripRateLgdVals=[0.001 0.005 0.01];
        coripRateLgdVals=[0.003 0.01 0.03];
        if zeroPh==2
            coripRateLgdVals=coripRateLgdVals*100;
        end
        
        %store to carry over
        coripTestLgd{iSet}=coripRateLgdVals;

        %ripple reference rates
        refRates=ripRate(~rejRegs,testBins,stimCondIdx(1));
        %single ripple legend increments
%         ripRateLgdVals = linspace(min(refRates,[],'all'),max(refRates,[],'all'),6)';
        ripRateLgdVals=linspace(0.02,0.1,6)';
%         ripRateLgdVals=ripRateLgdVals/densFctr;%correction if density metric
        %store to carry over
        rippleTestLgd{iSet}=ripRateLgdVals;
    else
%         coripRateLgdVals=coripTestLgd{iSet};
%         ripRateLgdVals=rippleTestLgd{iSet};
        coripRateLgdVals=coripTestLgd{2};
        ripRateLgdVals=rippleTestLgd{2};
    end
    
    %single ripple rate legend labels
    plotRoi(end-5:end)=arrayfun(@(x) num2str(x,'%.2f'),ripRateLgdVals,'uni',0);
    
    %co-ripple legend connectivity
    lgdGrph=zeros(6,6);
    lgdGrph(1,2)=1;
    lgdGrph(3,4)=1;
    lgdGrph(5,6)=1;
    lgdGrph(lgdGrph~=0)=coripRateLgdVals;
    lgdGrph=mirrormat(lgdGrph);
    lgdGrph(isnan(lgdGrph))=eps;
    
    lgdLab=arrayfun(@(x) num2str(x,'%.3f'),coripRateLgdVals,'uni',0);
    
    for c=1:2
        coripCond=coripRate(:,:,testBins,stimCondIdx(c));
        %corip(repmat(all(pCoripBase(:,:,:,stimCondIdx(1))>pThreshBase(iSet),3),[1 1 nBin]))=0;
        coripCond(~(pCoripCompBoot(:,:,:,iSet)<=alpha))=0;
        coripCond(rejRegs,:,:)=0;
        coripCond(:,rejRegs,:)=0;
        coripCond(isnan(coripCond))=0;
        
        
        rip=ripRate(:,testBins,stimCondIdx(c));
        rip(isnan(rip))=0;
        for t=1:numel(testBins)
            h{t,c}=figure('Position',[1 1 1080 1080]);
            
            idxVec=t;
            
            tCorip=coripCond(plotRegs,plotRegs,idxVec);
            tCorip(1:(size(tCorip,1)+1):end)=0;
            tCorip=diagConcat(tCorip,lgdGrph);
            tCorip(isnan(tCorip))=0;
%             tCorip(tCorip>0.05)=0.05;
            edgeLab = [repmat({''},sum(~~tCorip,'all')/2 - 3,1); lgdLab'];
            G{t,c}=graph(tCorip,plotRoi,'omitselfloops');
            
            tRip{t,c}=[nanmean(rip(plotRegs,idxVec),2);ripRateLgdVals];
            
            %emptyNodes{t,c}=cellfun(@(x) ~any(strcmp(x,G{t,c}.Edges.EndNodes(:))),strRoi);%|strcmp(G{t,c}.Nodes.Name,'rh-SPL');
            emptyNodes{t,c}=[rejRegs(plotRegs)|~any(pCoripCompBoot(:,:,t,iSet)<=alpha,2) ;false(6,1)];
            G{t,c}=rmnode(G{t,c},find(emptyNodes{t,c}));
            tRip{t,c}(emptyNodes{t,c})=0;
            
            image((1:imSz(2))/sclFctr, (1:imSz(1))/sclFctr, brainImg)
            hold on
            Gplot{t,c}=plot(G{t,c}, ...'NodeColor', nodeClr(~emptyNodes{t,c},:),...
                'EdgeColor',repmat([0.5 0.5 0.5], numel(G{t,c}.Edges.Weight),1),'EdgeAlpha',0.7,...
                'NodeFontSize',[18*ones(height(G{t,c}.Nodes)-6,1);16*ones(6,1)],...
                'EdgeFontSize',18,...
                'MarkerSize',10,...
                'NodeLabelColor',[1 1 1], 'EdgeLabel', edgeLab,'EdgeLabelColor',[1 1 1]);
            set(gca,'GridColor',[0 0 0])
        end
    end
    
    for t=1:nTestBin
        for c=1:2
            edgeStr{t,c} = arrayfun(@(x) [G{t,c}.Edges.EndNodes{x,1} G{t,c}.Edges.EndNodes{x,2}],(1:height(G{t,c}.Edges))','uni',0);
        end
        
        [~,iNode{t,1},iNode{t,2}]=intersect(edgeStr{t,1},edgeStr{t,2});
        for c=1:2
            %set line width normalized by maximum
            Gplot{t,c}.LineWidth = edgeFctr*G{t,c}.Edges.Weight/coripRateLgdVals(end);%max(cell2mat(cellfun(@(x) max(x.Edges.Weight),G(:,1),'uni',0)),[],'all');
            
            othCond=mod(c,2)+1;
            % store unnormalized edge colors
            edgeClr{t,c}=log10(G{t,c}.Edges.Weight(iNode{t,c})...
                ./G{t,othCond}.Edges.Weight(iNode{t,othCond}))*(3-2*c);
            edgeClr{t,c}=edgeClr{t,c}';
            
            % store unnormalized node colors
            nodeClr{t,c}=log10(tRip{t,c}...
                ./tRip{t,othCond})*(3-2*c);
            
            %             [~,i2,i1]=intersect(cellfun(@(x) x(4:end), Gplot{t,c}.NodeLabel,'uni',0)',regLabs);
            [~,i2,i1]=intersect(Gplot{t,c}.NodeLabel,plotRoi);
            
            Gplot{t,c}.XData(i2) = locs(i1,1);
            Gplot{t,c}.YData(i2) = locs(i1,2);
            
            Gplot{t,c}.NodeLabel(i2)=plotRoi(i1);
%             Gplot{t,c}.NodeFontSize=36;
        end
    end
    nodeSzNorm=max(abs([tRip{:}]),[],'all');
    nodeClrNorm=max(abs([nodeClr{:}]),[],'all');
%     edgeNorm=max(abs([edgeClr{:}]),[],'all');
    edgeNorm=clrRange;
    
    
    %sgtitle('co-rippling probability graphs for word vs CS')
    for t=1:nTestBin
        
        % set edge colors
        for c=1:2
            edgeClr{t,c}=127*edgeClr{t,c}/edgeNorm + 128;
            edgeClr{t,c}(edgeClr{t,c}>255)=255;
            edgeClr{t,c}(edgeClr{t,c}<0.5)=1;
            edgeClr{t,c}=cmap(round(edgeClr{t,c}),:);
            nEdge=numel(Gplot{t,c}.LineWidth);
            clrTmp=repmat(cmap(256-255*(c-1),:),[nEdge 1]);
            clrTmp(iNode{t,c},:)=edgeClr{t,c};
            Gplot{t,c}.EdgeColor(:,:)=clrTmp;
            
            
            % set node sizes
            Gplot{t,c}.MarkerSize = nodeFctr*sqrt(tRip{t,c}(~emptyNodes{t,c})/ripRateLgdVals(end)) + eps;%nodeSzNorm + eps;
            
            % set node colors
            %    nodeClr{t}= round(127*nodeClr{t}/nodeClrNorm + 128);
            nodeClr{t,c}= round(127*nodeClr{t,c}/edgeNorm + 128);
            nodeClr{t,c}(nodeClr{t,c}>255)=255;
            nodeClr{t,c}(nodeClr{t,c}<1)=1;
            nodeClr{t,c}(isnan(nodeClr{t,c}))=128;
            Gplot{t,c}.NodeColor = cmap(nodeClr{t,c}(~emptyNodes{t,c}),:);
            
            figure(h{t,c})
            figClr=0.8*ones(1,3);
            set(gca,'color',figClr)
            set(gcf,'color',figClr)
            set(gca,'XColor',figClr,'YColor',figClr,'TickDir','out')
            othCond=mod(c,2)+1;
            title(sprintf('%s, %s>%s, %s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{c},stimLab{othCond},timeLab{t}))
            cbar=colorbar;
            colormap jet
            cmap=colormap;
            ylim([0 imSz(1)/sclFctr])
            xlim([0 imSz(2)/sclFctr])
            if c==2
                colormap(flipud(cmap));
            end
            caxis([-edgeNorm edgeNorm])
            cbar.Label.String=sprintf('log_{10}  %s rate/%s',eveName,stimLab{mod(c,2)+1});
            cbar.Color=[1 1 1];
            cbar.FontSize=18;
        end
    end
    
    
    figOutDir=sprintf('%s/%sv%s%s_%s/',outDir,stimLab{:},apndTest,hemi);
    if ~exist([figOutDir '/noLgd/'],'dir');
        unix(sprintf('mkdir -p %s',[figOutDir '/noLgd/']));
    end
    if saveFigs
        for t=1:nTestBin
            %set(h{t,2},'Visible','off')
            if testDir==1
                printCond=1;
                close(h{t,2})
            elseif testDir==2
                printCond=2;
                close(h{t,1})
            elseif testDir==3
                printCond=[1 2];
            end
            for c=printCond%:2
                set(h{t,c},'InvertHardCopy','off','Color',[0 0 0])
                figure(h{t,c})
                drawnow
                set(gcf,'Units','Pixels')
                set(gcf,'Position',get(h{end,printCond},'Position'))
                drawnow
                pause(1)
                set(gcf,'Position',get(h{end,printCond},'Position'))
                drawnow
                ax=gca;
                set(ax,'XColor',[0 0 0],'YColor',[0 0 0])
                ax.Title.Color=[1 1 1];
                ax.Title.FontSize=[32];
                nodeLabs=Gplot{t,c}.NodeLabel;
                replotLabs=find(cellfun(@(x) any(strcmp(nodeLabs,x)),regLabsShrt));
                Gplot{t,c}.NodeLabel(1:end-6)=arrayfun(@(x) repmat(' ',x,1), 1:numel(Gplot{t,c}.NodeLabel)-6,'uni',0);
                for r=1:numel(replotLabs)
                    rIdx=replotLabs(r);
%                     text(labLocs(rIdx,1),labLocs(rIdx,2),regLabsShrt(rIdx),'color',[1 1 1],'FontSize',55,'FontName','Arial');
                    text(labLocs(rIdx,1),labLocs(rIdx,2),regLabsShrt(rIdx),'color',[1 1 1],'FontSize',40,'FontName','Arial');
                end
                saveGraphic(h{t,c},sprintf('%s/%s_%s_%sv%s_%sGraph',figOutDir,timeLab{t},stimLab{c},stimLab{1:2},eveName),[],600);
                
                Gplot{t,c}.XData(:,end-5:end)=-1000;
                title(timeLab{t},'FontSize',32)
                drawnow
                print(h{t,c},sprintf('%s/noLgd/%s_%s_%sv%s_%sGraph.png',figOutDir,timeLab{t},stimLab{c},stimLab{1:2},eveName),'-dpng','-r400')
            end
        end
%         if fuseBins==1
%             h6=figure('Position',[1 39 1920 520+80*uint16(resp)]);
%         else
%             h6=figure('Position',get(0,'Screensize'));
            h6=figure('Position',[1 39 1920 400]);
%         end
        set(gcf,'InvertHardCopy','off','Color',[0 0 0])
        imSzAll=[0 0];
        for t=1:nTestBin
%             if t==1
%                 img=imread(sprintf('%s/png/%s_%s_%sv%s_%sGraph.png',figOutDir,timeLab{t},stimLab{printCond(1)},stimLab{1:2},eveName));
%                 img(:,end-670:end,:)=0;
%                 img(:,[1:420 end-250:end],:)=[];
%             elseif t==nTestBin
%                 img=imread(sprintf('%s/noLgd/%s_%s_%sv%s_%sGraph.png',figOutDir,timeLab{t},stimLab{printCond(1)},stimLab{1:2},eveName));
%                 img(:,[1:420 end-125:end],:)=[];
%             else
                img=imread(sprintf('%s/noLgd/%s_%s_%sv%s_%sGraph.png',figOutDir,timeLab{t},stimLab{printCond(1)},stimLab{1:2},eveName));

                
                img(:,[1:500 end-700:end],:)=[];
%                 img(:,[1:500 end-670:end],:)=[];
                img([1:630 end-600:end],:,:)=[];
                imSzAll(1)=imSzAll(1)+size(img,1);
                imSzAll(2)=size(img,2);
%             end
            
            subtightplot(1,nTestBin,t,[0 0],[0 0],[0 0])
            figClr=[0 0 0];
            set(gcf,'color',figClr)
            set(gca,'XColor',figClr,'YColor',figClr,'color',figClr,'GridColor',figClr)
%             if t==1
% %                 title(sprintf('%s, %s>%s, %s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{testDir},stimLab{mod(testDir,2)+1},timeLab{t}))
%                 text(size(img,1)/2+250,80,...
%                     sprintf('%s, %s>%s, %s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{testDir},stimLab{mod(testDir,2)+1},timeLab{t}),...
%                     'HorizontalAlignment','center','Color',[1 1 1],'FontSize',32)
%             else
% %                 title(timeLab{t},'FontSize',32)
%                 text(size(img,1)/2+250,80,...
%                     timeLab{t},...
%                     'HorizontalAlignment','center','Color',[1 1 1],'FontSize',32)
%             end
            ax=gca;
            ax.Title.Color=[1 1 1];
            imshow(img,'Border','tight','InitialMagnification','fit')
            if fuseBins
                fontSz=16;
            else
                fontSz=10;
            end

%             fontSz=24;
            if t==1
%                 title(sprintf('%s, %s>%s, %s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{testDir},stimLab{mod(testDir,2)+1},timeLab{t}))
                text(size(img,1)/2+350,100,...
                    sprintf('%s, %s>%s\n%s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{testDir},stimLab{mod(testDir,2)+1},timeLab{t}),...
                    'HorizontalAlignment','center','Color',[1 1 1],'FontSize',fontSz,'FontWeight','bold','VerticalAlignment','bottom')
%                                 text(size(img,1)/2+350,100,...
%                     sprintf('%s, %s>%s, %s',[upper(eveName(1)) eveName(2:end) 's'], stimLab{testDir},stimLab{mod(testDir,2)+1},timeLab{t}),...
%                     'HorizontalAlignment','center','Color',[1 1 1],'FontSize',fontSz,'FontWeight','bold')%,'VerticalAlignment','bottom')
            else
%                 title(timeLab{t},'FontSize',32)
                text(size(img,1)/2+300,100,...
                    timeLab{t},...
                    'HorizontalAlignment','center','Color',[1 1 1],'FontSize',fontSz,'FontWeight','bold')%,'VerticalAlignment','bottom')
            end
        end
        drawnow
        set(h6,'Position',[1 1  imSzAll/4])
        print(h6,sprintf('%s/%sv%s_%s_allLatency.png',figOutDir,stimLab{1:2},eveName),'-dpng','-r800')
    end
    
end

    end
end

end