function normalizedSpectrumPlot(spectra,fVec,oof,xLabel,yLabel,subTitles,nSbpltHrz,legendText,colors,plotOof,skipChans, specErr, lines)

gap=[0.04 0.030];%[0.02 0.015]; %vertical, horizontal
h_marg=[0.05 0.09]; %bottom, top
w_marg=[0.03 0.02]; %left, right

nChans=size(spectra,2);

plotErr=false;
if exist('specErr','var')
    plotErr=true;
end
if ~exist('lines','var')
    lines=[];
end
if ~exist('plotOof','var')
    plotOof=false;
end
if ~exist('skipChans','var')
    skipChans=false(nChans,1);
end
if ~exist('colors','var')||isempty(colors)
    colors=[0, 0.4470, 0.7410;...
        0.8500, 0.3250, 0.0980;...
        0.9290, 0.6940, 0.1250;...
        0.4940, 0.1840, 0.5560;...
        0.4660, 0.6740, 0.1880;...
        0.3010, 0.7450, 0.9330;...
        0.6350, 0.0780, 0.1840];
end


% if any(isnan(spectra(:,~skipChans,:)),'all')
%     error('meanSpectra matrix has nan elements')
% end
if size(oof,2)~=nChans
    error('number of OOF fits doesn''t match number of spectra')
end
if size(spectra,1)~=size(fVec,1)
    error('number of frequency values in spectrum does not match fVec')
end

subtightplot(ceil((nChans+1)/nSbpltHrz), nSbpltHrz, 1, gap, h_marg, w_marg)
plot(nan,nan);hold on
lgndClrs=unique(colors,'rows','stable');
for p=1:size(lgndClrs,1)
    plot(nan,nan,'Color',lgndClrs(p,:))
end
hold off

ax=gca;
ax.Visible=false;
[lgd,objH]=legend(legendText,...
    'Location','none','Orientation','vertical');
set(findobj(objH,'Tag',legendText{1}),'Vis','off')
pos=get(objH(1),'Pos');
set(objH(1),'Pos',[0.1 pos(2:3)], 'String', legendText{1})
subPltDim=ax.Position;
subPltDim(1:2)=subPltDim(1:2)+subPltDim(3:4)/4;
subPltDim(3:4)=subPltDim(3:4)/2;
lgd.Position=subPltDim;


if ~plotOof %~plotOof means normalize by oof
    oofTmp=repmat(oof,[1 1 size(spectra,3)]);
    if plotErr
        specErr=cat(4,(specErr+spectra)./oofTmp,(specErr-spectra)./oofTmp); %top and bottom error bars are now different; note that specErr is now in absolute terms
        spectra=spectra./oofTmp;
        specErr=specErr-repmat(spectra,[1 1 1 2]); %shift specErr to relative error
    end
    %if plotErr;specErr=specErr.*spectra;end
end

for n=1:nChans
    subtightplot(ceil((nChans+1)/nSbpltHrz), nSbpltHrz, n+1, gap, h_marg, w_marg);
    if skipChans(n)||all(isnan(spectra(:,n,:)),[1 3])
        continue
    end
    

    for iCond=1:size(spectra,3)
        if all(isnan(spectra(:,n,iCond)));continue;end
        if ~plotErr
            loglogjg(fVec, spectra(:,n,iCond),[],'Color',colors(iCond,:));
        else
            loglogjg(fVec, spectra(:,n,iCond),squeeze(specErr(:,n,iCond,:)),'Color',colors(iCond,:,:));
        end
        hold on
    end
    if plotOof;loglogjg(fVec, oof(:,n), [], 'LineStyle', '--', 'Color','k');end
    ylim(log10([min(spectra(:,n,:),[],'all') max(spectra(:,n,:),[],'all')]))
    for iLine=1:numel(lines)
         line(log10([lines(iLine) lines(iLine)]),ylim,'Color',[0 1 0])
    end
    
    ax=gca;ax.FontSize=7;
    ax.Children=flip(ax.Children);
    xlim(log10([fVec(1) fVec(end)+eps]))
    
    if n+1<=nChans-(nSbpltHrz-1)
        %ax.XTickLabel=[];
    else
        xlabel(xLabel)
    end
    if n==1||mod(n+1,nSbpltHrz)==1
        ylabel(yLabel);
    end
    title(subTitles{n})
end
end

