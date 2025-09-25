function WNcoripRatePerChan(subjIdcs,trlOpt)

% if trlOpt
% outDir='/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/coripRate/NY226FWIO/';
% subIdx=1:1880;
% else
% outDir='/space/seh8/5/halgdev/projects/jgarret/FW_data/taskPreproc/coripRate/NY226FWvar/';
% subIdx=1881:3560;
% end
subIdx=[];
task='WN';
plotChans=1;
nCorSamp=4000;
corWin=[1:nCorSamp]-nCorSamp/2;
%%
if nargin==0
%     subjIdcs=[311 297  233 226 190 133 130 111 101 -89 86 68];
    subjIdcs=[49 51];
end
subjects=genSubjList(subjIdcs,'MG');

ripXcorAll=[];
ripXcorCentAll=[];
ripXcorPkAll=[];
pCoRall=[];
pCoRtri=[];
distAll=[];
distSymAll=[];
relRateAll=[];
nCoRall=[];
PLVall=[];
close all
outFigDir=sprintf('/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/WN/rateDat/');
outStatDir='/space/seh8/6/halgdev/projects/jacob/auditory_ripples/coripDat/summaryDat/';
if ~exist(outFigDir,'dir');unix(sprintf('mkdir -p %s',outFigDir));end
if ~exist(outStatDir,'dir');unix(sprintf('mkdir -p %s',outStatDir));end
for s=1:numel(subjects)
    subj=subjects{s}

load(sprintf('/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/%s/tasks/%s_WN.mat',subj,subj),'data','label','Fs')
% ripAmpAll=[rippleStats.rippleAmp{:}];
% [ripAmpAll,idx]=sort(ripAmpAll,'descend');
nS=size(data,1);



[trl,resp,rej]=loadTrials(subj,task);
trl=trl(~rej,:);

codes={1, 2, [1 2], [3 4]};

condNames={'match','mismatch','word','noise'};
gTrl=trl(:,1);
[labParc,~,~,~,~,rejChan,coords]=retrieveChanParcels_MG(subj);
rmCh =[rejChan; true(size(labParc,1)-numel(rejChan),1)]|strcmp(labParc(:,1),'N/A')|strcmp(labParc(:,3),'rh');



labParc(rmCh,:)=[];
coords(rmCh,:)=[];
data(:,rmCh)=[];
nCh=sum(~rmCh);


[b,a]=butter(3,[70 190]/(1000/2));
HGenv=abs(hilbert(filtfilt(b,a,double(data))));
[b,a]=butter(3,[70 110]/(1000/2));
ripPh=angle(hilbert(filtfilt(b,a,double(data))));

% words=trl(any(abs(trl(:,2))==codes{1}),1);
% targs=trl(any(abs(trl(:,2))==codes{4}),1);

%     coords=retrieveChanCoords(subjects{s},'lh');
%     coords=cell2mat(cellfun(@(x) cell2mat(coords(strcmp(coords(:,1),x),2:end)),labParc(:,4),'uni',0));

%%
% load(sprintf('/space/seh8/6/halgdev/projects/jacob/NC_ripple/matFiles/MG/%s_ripple_stats_%s_0_cln.mat',subj,task),'rippleStats')
load(sprintf('/space/seh8/6/halgdev/projects/jacob/NC_ripple/matFiles/MG/%s_ripple_stats_%s_0_trlSlct.mat',subj,task),'rippleStats')
%
% rip_mask = zeros(size(data));
rip_bounds_mask = false(size(data));
rip_peak_mask=false(size(data));
gCh=find(~rmCh)';

% [b2,a2]=butter(3,[70 110]/(Fs/2),'bandpass');
% RB=filtfilt(b2,a2,double(data));
% RBamp=abs(hilbert(RB));
% baseMsk=bounds2mask([-100 50]+gTrl,nS);
% meanRB=mean(RBamp(baseMsk,:),1);
% stdRB=std(RBamp(baseMsk,:),0,1);
% 

for ch = gCh
    %     rip_mask(ch,rippleStats.locs{ch}) = true;
    for rip=1:size(rippleStats.window{ch},1)
        rip_bounds_mask(rippleStats.window{ch}(rip,1):rippleStats.window{ch}(rip,2),ch==gCh) = true;
        rip_peak_mask(rippleStats.locs{ch}(rip),ch==gCh) = true;
    end
end


trlMsk=bounds2mask([-200 1600]+round(gTrl),nS);
trlMsk(:)=true;
% trlMskMap=find(trlMsk);
nTrlSamp=sum(trlMsk);
rip_bounds_mask(~trlMsk,:)=[];
rip_peak_mask(~trlMsk,:)=[];
ripPh(~trlMsk,:)=[];
corip_ph = single(ripPh) - permute(single(ripPh),[1 3 2]);
corip_bounds_mask = rip_bounds_mask & permute(rip_bounds_mask,[1 3 2]);
coripBnds=cell(nCh,nCh);
coripPh = cell(nCh,nCh);
ripBnds=cell(nCh,1);
corip_cent_mask = false(nTrlSamp,nCh,nCh);
rip_cent_mask = false(nTrlSamp,nCh);
ripCents=cell(nCh);
ripPeaks=cell(nCh);

ripXcor=nan(nCorSamp,nCh,nCh);
ripXcorCent=nan(nCorSamp,nCh,nCh);
ripXcorPk=nan(nCorSamp,nCh,nCh);
nCoR = nan(nCh,nCh);
pCoR = nan(nCh,nCh);
relRateCoR = nan(nCh,nCh);
nRip = nan(nCh,1);
nRipCent = nan(nCh,1);
chanPLV=nan(nCh,nCh);

dist=nan(nCh,nCh);

h9=figure('Position',[1 1 1000 500])
for c1=1:nCh
    ripBnds{c1}=mask2bounds(rip_bounds_mask(:,c1));
    nRip(c1)=size(ripBnds{c1},1);
    
    rip_cent_mask(:,c1) = bounds2mask(floor(mean(ripBnds{c1},2)),nTrlSamp);
    ripCents{c1}=find(rip_cent_mask(:,c1));
    ripCents{c1}(ripCents{c1}<=nCorSamp/2|ripCents{c1}>nTrlSamp-nCorSamp/2)=[];
    nRipCent(c1)=length(ripCents{c1});
    
    ripPeaks{c1}=find(rip_peak_mask(:,c1));
    ripPeaks{c1}(ripPeaks{c1}<=nCorSamp/2|ripPeaks{c1}>nTrlSamp-nCorSamp/2)=[];
    nRipPeak(c1)=length(ripPeaks{c1});
    for c2=1:c1-1
        coripBnds{c1,c2}=mask2bounds(corip_bounds_mask(:,c1,c2));
        coripBnds{c1,c2}(diff(coripBnds{c1,c2},1,2)<25,:)=[];
        
        nCoR(c1,c2)=size(coripBnds{c1,c2},1);
        coripPh{c1,c2}=nan(nCoR(c1,c2),1);
        for cr=1:nCoR(c1,c2)
            tmpVals=corip_ph(coripBnds{c1,c2}(cr,1):coripBnds{c1,c2}(cr,2),c1,c2);
            tmpVals = cos(tmpVals) + 1i*sin(tmpVals);
            coripPh{c1,c2}(cr) = angle(mean(tmpVals));
        end
        chanPLV(c1,c2)=abs(mean(cos(coripPh{c1,c2})+1i*sin(coripPh{c1,c2})));
        
        
        corip_cent_mask(:,c1,c2) = bounds2mask(floor(mean(coripBnds{c1,c2},2)),nTrlSamp);
        if nRipCent(c1)
            ripXcor(:,c1,c2) = mean(cell2mat(arrayfun(@(x) rip_bounds_mask(corWin+ripCents{c1}(x),c2),1:nRipCent(c1),'uni',0)),2)/mean(rip_bounds_mask(:,c2));
            ripXcorCent(:,c1,c2) = mean(cell2mat(arrayfun(@(x) rip_cent_mask(corWin+ripCents{c1}(x),c2),1:nRipCent(c1),'uni',0)),2)/mean(rip_cent_mask(:,c2));
            ripXcorPk(:,c1,c2) = mean(cell2mat(arrayfun(@(x) rip_peak_mask(corWin+ripPeaks{c1}(x),c2),1:nRipPeak(c1),'uni',0)),2)/mean(rip_peak_mask(:,c2));
            
            pCoR(c1,c2) = mean(arrayfun(@(x) any(rip_bounds_mask(ripBnds{c1}(x,1):ripBnds{c1}(x,2),c2)),1:nRip(c1)));
                
        end
        if nRipCent(c2)
            pCoR(c2,c1) = mean(arrayfun(@(x) any(rip_bounds_mask(ripBnds{c2}(x,1):ripBnds{c2}(x,2),c1)),1:nRip(c2)));
        end
        
        relRateCoR(c1,c2) = sum(corip_cent_mask(:,c1,c2))/(nRip(c1)*nRip(c2));
            
%         pCoR(c1,c2) = nCoR(c1,c2)/(nRip(c1)*nRip(c2));
        
        
        dist(c1,c2)=norm(diff(coords([c1 c2],:),1,1));
    end
    
gap=[0.07 0.05];
marg_h=[0.09 0.07];
marg_w=[0.08 0.05];
    
    binSamp=200;

    
    if plotChans
    subtightplot(2,1,1,gap,marg_h,marg_w)
    hold on
    psh=-1000:2000;
    nHsmp=numel(psh);
%     rsVecRip = [binSamp/2:binSamp:nHsmp-(3/2)*binSamp; (3/2)*binSamp:binSamp:nHsmp-binSamp/2]';
    rsVecRip = [1:binSamp:nHsmp-binSamp; binSamp:binSamp:nHsmp]';
    nTrlCond = cellfun(@(x) sum(any(trl(:,2)==x,2)), codes);
%     lgdStrs=condNames(1:4);
    lgdStrs=arrayfun(@(x) [condNames{x} ' ' num2str(nTrlCond(x))], 1:4,'uni',0);
    inclCond=[];
    for cond=1:4
        condTrl=trl(any(trl(:,2)==codes{cond},2),1);
%         condTrl=arrayfun(@(x) find(trlMskMap==x),condTrl);
        condTrl(condTrl+psh(1)<1|condTrl+psh(end)>nTrlSamp)=[];
        if isempty(condTrl)
            lgdStrs(cond)=[];
            continue
        else
            inclCond=[inclCond cond];
        end
        condRip = mean(cell2mat(arrayfun(@(x) rip_cent_mask(psh+x,c1),condTrl','uni',0)),2);%/mean(rip_cent_mask(:,c1));
        condRipBin=cell2mat(arrayfun(@(x) sum(condRip(rsVecRip(x,1):rsVecRip(x,2),:),1),1:size(rsVecRip,1),'uni',0)');
        plot(psh(binSamp/2:binSamp:end-binSamp/2),condRipBin,'LineWidth',2);
    end
    xticks(psh(1:binSamp:end))
    xlim(psh([1 end]))
    title([labParc{c1,1} ' ' labParc{c1,4}])
    ylabel('N ripple/N trial')
    xlabel('ms relative to stimulus')
    line([0 0],ylim,'Color','k','LineStyle','--')
    line([1 1]*-500,ylim,'Color','k','LineStyle','--')
    lgd=legend(lgdStrs);
    lgd.FontSize=7;
    
    
    subtightplot(2,1,2,gap,marg_h,marg_w)
    hold on
    for cond=inclCond
        condTrl=trl(any(trl(:,2)==codes{cond},2),1);
        condHG = mean(cell2mat(arrayfun(@(x) HGenv(psh+x,c1),condTrl','uni',0)),2);
%         condHG=smoothdata(condHG,'gaussian',50);
        plot(psh,condHG,'LineWidth',2);
%         condHGbin = cell2mat(arrayfun(@(x) mean(condHG(rsVecRip(x,1):rsVecRip(x,2),:),1),1:size(rsVecRip,1),'uni',0)');
%         plot(psh(binSamp/2:binSamp:end-binSamp/2),condHGbin,'LineWidth',2);
    end
    xticks(psh(1:binSamp:end))
    xlim(psh([1 end]))
%     title([labParc{c1,1} ' ' labParc{c1,4}])
    ylabel('uV 70-190Hz')
    xlabel('ms relative to stimulus')
    line([0 0],ylim,'Color','k','LineStyle','--')
    line([1 1]*-500,ylim,'Color','k','LineStyle','--')
    legend(lgdStrs);
    lgd.FontSize=7;
    
%     subtightplot(3,1,3,gap,marg_h,marg_w)
%     hold on
%     psh=-800:400;
%     nHsmp=numel(psh);
% %     rsVecRip = [binSamp/2:binSamp:nHsmp-(3/2)*binSamp; (3/2)*binSamp:binSamp:nHsmp-binSamp/2]';
%     rsVecRip = [1:binSamp:nHsmp-binSamp; binSamp:binSamp:nHsmp]';
%     for cond=[2 -2]
%         condTrl=resp(any(resp(:,2)==cond,2),1);
%         if isempty(condTrl)
%             continue
%         end
% %         condTrl=arrayfun(@(x) find(trlMskMap==x),condTrl);
%         condTrl(condTrl+psh(1)<1|condTrl+psh(end)>nTrlSamp)=[];
%         condRip = mean(cell2mat(arrayfun(@(x) rip_cent_mask(psh+x,c1),condTrl','uni',0)),2);%/mean(rip_cent_mask(:,c1));
%         condRipBin=cell2mat(arrayfun(@(x) sum(condRip(rsVecRip(x,1):rsVecRip(x,2),:),1),1:size(rsVecRip,1),'uni',0)');
% % %         condHG = mean(cell2mat(arrayfun(@(x) HGenv(psh+x,c1),condTrl','uni',0)),2);
% %         condHGbin = cell2mat(arrayfun(@(x) mean(condHG(rsVecRip(x,1):rsVecRip(x,2),:),1),1:size(rsVecRip,1),'uni',0)');
%         plot(psh(binSamp/2:binSamp:end-binSamp/2),condRipBin,'LineWidth',2);
% %         yyaxis right
% %         plot(psh(binSamp/2:binSamp:end-binSamp/2),condHGbin,'LineWidth',2);
%         
%     end
%     xticks(psh(1:binSamp:end))
%     ylabel('N ripple/N trial')
%     xlabel('ms relative to response')
%     line([0 0],ylim,'Color','k','LineStyle','--')
%     legend({'TP','FP'})
    saveGraphic(h9,sprintf('%s/%s/%s_ripMean',outFigDir,subj,labParc{c1,4}))
    clf
    end
end

ripXcorPr = reshape(ripXcor(repmat(permute(tril(true(nCh),-1),[3 1 2]),nCorSamp,1)),[nCorSamp nCh*((nCh-1)/2)]);
ripXcorCentPr = reshape(ripXcorCent(repmat(permute(tril(true(nCh),-1),[3 1 2]),nCorSamp,1)),[nCorSamp nCh*((nCh-1)/2)]);
ripXcorPkPr = reshape(ripXcorPk(repmat(permute(tril(true(nCh),-1),[3 1 2]),nCorSamp,1)),[nCorSamp nCh*((nCh-1)/2)]);
ripXcorAll=[ripXcorAll ripXcorPr];
ripXcorCentAll=[ripXcorCentAll ripXcorCentPr];
ripXcorPkAll=[ripXcorPkAll ripXcorPkPr];
% distPr=dist(tril(true(nCh),-1));
% distAll=[distAll;distPr];


nanMskTri=isnan(dist);
distAll=[distAll; dist(~nanMskTri)];
dist(isnan(dist))=0;
distSym =mirrormat(dist);
dist(dist==0)=nan;
distSym(distSym==0)=nan;

distSymPair = distSym(~isnan(pCoR));
distSymAll = [distSymAll; distSymPair];
pCoRpair=pCoR(~isnan(pCoR));
pCoRall=[pCoRall;pCoRpair];

nCoRpair=nCoR(~nanMskTri);
nCoRall=[nCoRall; nCoRpair];
PLVpair=chanPLV(~nanMskTri);
PLVall=[PLVall;PLVpair];
pCoRtri=[pCoRtri;pCoR(~nanMskTri)];


relRateAll=[relRateAll; relRateCoR(~nanMskTri)];


end
save(sprintf('%s/coripRateDat.mat',outStatDir),'ripXcorAll','ripXcorCentAll','ripXcorPkAll','pCoRall','distAll','relRateAll','distSymAll','nCoRall','PLVall','pCoRtri')

%%
nanMsk=any(isnan(ripXcorAll),1);
ripXcorAll(:,nanMsk)=[];
nanMskSym=isnan(pCoRall);
pCoRall(nanMskSym)=[];
distAll(nanMskSym)=[];
nPr=size(ripXcorAll,1);

ripXcorMean=nanmean(ripXcorAll,2);
ripXcorSEM = nanstd(ripXcorAll,1,2)/sqrt(nPr);
% ripXcorCentSmth=smoothdata(ripXcorCentAll,'gaussian',20);
binSamp=20;
rsVec = [binSamp/2:binSamp:nCorSamp-(3/2)*binSamp; (3/2)*binSamp:binSamp:nCorSamp-binSamp/2]';

ripXcorCentSmth=cell2mat(arrayfun(@(x) mean(ripXcorCentAll(rsVec(x,1):rsVec(x,2),:),1),1:size(rsVec,1),'uni',0)');
% ripXcorCentSmth=nandownsample(ripXcorCentAll,1,10);
ripXcorCentMean=nanmean(ripXcorCentSmth,2);
ripXcorCentSEM = nanstd(ripXcorCentSmth,1,2)/sqrt(nPr);
%%
% ripXcorPkSmth=nandownsample(ripXcorPkAll(10:end-10,:),1,20);
ripXcorPkSmth=cell2mat(arrayfun(@(x) mean(ripXcorPkAll(rsVec(x,1):rsVec(x,2),:),1),1:size(rsVec,1),'uni',0)');
ripXcorPkMean=nanmean(ripXcorPkSmth,2);
ripXcorPkSEM = nanstd(ripXcorPkSmth,1,2)/sqrt(nPr);
%
%% plot
h1=figure();
shadedErrorBar(corWin(1501:2500),ripXcorMean(1501:2500),ripXcorSEM(1501:2500));
ylabel('co-ripple density/chance')
xlabel('ms vs ripple center')
xlim([corWin(1501)-1 corWin(end-1500)])
xticks(-nCorSamp/2:100:nCorSamp/2)
prettifyPlot(h1,gca)
saveGraphic(h1,sprintf('%s/ripCorrelogram',outFigDir))
%%
%{
h0=figure();
shadedErrorBar(corWin(1:10:end),ripXcorCentMean,ripXcorCentSEM);
ylabel('channel 2 ripple centers/chance')
xlabel('ms vs channel 1 ripple centers')
% xlim([corWin(1)-1 corWin(end)])
xlim([-500 500])
xticks(-nCorSamp/2:250:nCorSamp/2)
% xticklabels(xticks*10)
prettifyPlot(h0,gca)
saveGraphic(h0,sprintf('%s/ripCentCorrelogram',outFigDir))
%}
%%
%{
h9=figure();
shadedErrorBar(corWin(binSamp:binSamp:end-binSamp),ripXcorPkMean,ripXcorPkSEM);
ylabel('channel 2 ripple peaks/chance')
xlabel('ms vs channel 1 ripple peaks')
% xlim([corWin(1)-1 corWin(end)])
xlim([-500 500])
% xticks(
xticks(-nCorSamp/2:250:nCorSamp/2)
% xticklabels(xticks*10)
prettifyPlot(h9,gca)
saveGraphic(h9,sprintf('%s/ripPkCorrelogram',outFigDir))
%}
%%
%{
h2=figure();
s=scatter(distSymAll/10,pCoRall,80,'.')
s.MarkerEdgeAlpha=0.5;
% scatter(distSymAll/10,pCoRall,2)
hold on
mdl=fitlm(distSymAll/10,pCoRall);
m=mdl.Coefficients.Estimate(2);
b=mdl.Coefficients.Estimate(1);
p=mdl.Coefficients.pValue(2);
rsq=mdl.Rsquared.Adjusted;
title(sprintf('pairs, b=%.5f R^{2}=%.4f p=%.4f',m,rsq,p))
line(xlim,xlim*m + b,'color','r')
xlabel('euclidian distance [cm]')
ylabel('p(corip|rip)')
prettifyPlot(h2,gca)
saveGraphic(h2,'/space/seh8/5/halgdev/projects/jgarret/taskAnalysis/figs/FW/coripRateAnalysis/coripVdist')
%}
%%
%
h21=figure();
% pCoRtri
PLVmagAll=abs(PLVall);
s=scatter(distAll(nCoRall>20)/10,PLVmagAll(nCoRall>20),8,'.')
% legend('pCoR<0.2','pCoR>0.2')
s.MarkerEdgeAlpha=0.5;
% scatter(distSymAll/10,pCoRall,2)
hold on
% mdl=fitlm(distSymAll/10,pCoRall);
mdl=fitlm(distAll(nCoRall>20)/10,PLVmagAll(nCoRall>20))
m=mdl.Coefficients.Estimate(2);
b=mdl.Coefficients.Estimate(1);
p=mdl.Coefficients.pValue(2);
rsq=mdl.Rsquared.Adjusted;
title(sprintf('pairs, b=%.5f R^{2}=%.4f p=%.4f',m,rsq,p))
line(xlim,xlim*m + b,'color','r')
xlabel('euclidian distance [cm]')
% ylabel('p(corip|rip)')
ylabel('PLV')
ylim([0 1])
prettifyPlot(h21,gca)
saveGraphic(h21,sprintf('%s/plvVdist',outFigDir))
%}
%%
h22=figure();
% pCoRtri
PLVmagAll=abs(PLVall);
s=scatter(pCoRtri(nCoRall>20)/10,PLVmagAll(nCoRall>20),8,'.')
s.MarkerEdgeAlpha=0.5;
% scatter(distSymAll/10,pCoRall,2)
hold on
% mdl=fitlm(distSymAll/10,pCoRall);
mdl=fitlm(pCoRtri(nCoRall>20)/10,PLVmagAll(nCoRall>20))
m=mdl.Coefficients.Estimate(2);
b=mdl.Coefficients.Estimate(1);
p=mdl.Coefficients.pValue(2);
rsq=mdl.Rsquared.Adjusted;
title(sprintf('pairs, b=%.5f R^{2}=%.4f p=%.4f',m,rsq,p))
line(xlim,xlim*m + b,'color','r')
xlabel('p(coR|R)')
% ylabel('p(corip|rip)')
ylabel('PLV')
ylim([0 1])
prettifyPlot(h22,gca)
saveGraphic(h22,sprintf('%s/PLVvPcoR',outFigDir))
%%
h3=figure();
scatter(distAll/10,relRateAll,2.5,'filled')
% scatter(distSymAll/10,pCoRall,2)
hold on
mdl=fitlm(distAll/10,relRateAll);
m=mdl.Coefficients.Estimate(2);
b=mdl.Coefficients.Estimate(1);
p=mdl.Coefficients.pValue(2);
rsq=mdl.Rsquared.Adjusted;
title(sprintf('pairs, b=%.5f R^{2}=%.4f p=%.4f',m,rsq,p))
line(xlim,xlim*m + b,'color','r')
xlabel('euclidian distance [cm]')
ylabel('N_{co-ripple(ch1,ch2){/(N_{ripple(ch1)}*N_{ripple(ch2)})')
prettifyPlot(h3,gca)
saveGraphic(h3,sprintf('%s/relCoripRateVdist',outFigDir))
%}

%%
end