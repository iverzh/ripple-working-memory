function multiChanPlot(dat,xVec,xLabel,yLabel,subTitles,nSbpltHrz,legendText,colors,skipChans, err, lines)

gap=[0.04 0.030];%[0.02 0.015]; %vertical, horizontal
h_marg=[0.05 0.09]; %bottom, top
w_marg=[0.03 0.02]; %left, right

nChans=size(dat,2);

plotErr=false;
if exist('err','var')
    plotErr=true;
end
if ~exist('lines','var')
    lines=[];
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
if size(dat,1)~=size(xVec,1)
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


for n=1:nChans
    subtightplot(ceil((nChans+1)/nSbpltHrz), nSbpltHrz, n+1, gap, h_marg, w_marg);
    if skipChans(n)||all(isnan(dat(:,n,:)),[1 3])
        continue
    end
    

    for iCond=1:size(dat,3)
        if all(isnan(dat(:,n,iCond)));continue;end
        if ~plotErr
            plot(xVec, dat(:,n,iCond),[],'Color',colors(iCond,:));
        else
            shadedErrorBar(xVec, dat(:,n,iCond),squeeze(err(:,n,iCond,:)),'lineProps',{'Color',colors(iCond,:,:)});
        end
        hold on
    end
    %ylim([min(dat(:,n,:),[],'all') max(dat(:,n,:),[],'all')])
    for iLine=1:numel(lines)
         line([lines(iLine) lines(iLine)],ylim,'Color',[0 0 0])
    end
    
    ax=gca;ax.FontSize=7;
    ax.Children=flip(ax.Children);
    
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

