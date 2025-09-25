function out = loglogjg(x,y,err,varargin)
%LOGLOGJG I'm at my wits' end with loglog displaying incomprehensible tick
%arrangements so I am programming this to automatically log-transform data,
%use matlab's linear plot function, and make sensical ticks.
logx=log10(x);
logy=log10(y);

if isempty(err)
    out=plot(logx, logy, varargin{:});
else
    if size(err,2)==1
        logErrPlus=log10(y+err)-logy;
        logErrMinus=log10(y-err)-logy;
    else
        logErrPlus=log10(y+err(:,1))-logy;
        logErrMinus=log10(y+err(:,2))-logy;
    end
    out=shadedErrorBar(logx,logy,[logErrPlus -logErrMinus],'lineProps',varargin);
end
xTicks=floor(min(xlim)):ceil(max(xlim));
xMinorTicks=nan((size(xTicks,2)-1)*9+1,1);
xTickLabels=repmat({''},[size(xMinorTicks,1),1]);
xTickLabels(1:9:end)=arrayfun(@(x) sprintf('10^{%d}',x),xTicks,'UniformOutput',false);
minorTickProps=log10(2:9);
for x=1:numel(xTicks)
    xMinorTicks(1+(x-1)*9)=xTicks(x);
    if x~=numel(xTicks)
        xMinorTicks((2:9)+(x-1)*9)=xTicks(x)+minorTickProps;
    end
end
xticks(xMinorTicks)
xticklabels(xTickLabels)
%xlim([min(logx) max(logx)])


yTicks=floor(ylim):ceil(max(ylim));
yMinorTicks=nan((size(yTicks,2)-1)*9+1,1);
yTickLabels=repmat({''},[size(yMinorTicks,1),1]);
yTickLabels(1:9:end)=arrayfun(@(y) sprintf('10^{%d}',y),yTicks,'UniformOutput',false);
minorTickProps=log10(2:9);
for y=1:numel(yTicks)
    yMinorTicks(1+(y-1)*9)=yTicks(y);
    if y~=numel(yTicks)
        yMinorTicks((2:9)+(y-1)*9)=yTicks(y)+minorTickProps;
    end
end
yticks(yMinorTicks)
yticklabels(yTickLabels)
%ylim([min(logy) max(logy)])

end

