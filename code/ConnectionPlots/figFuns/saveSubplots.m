function saveSubplots(h,baseFileName,saveType,res)
if nargin<4
    res=0;
end

lgnd=findobj(h,'Type','Legend');
ax=findobj(h,'Type','Axes');
allAxPos= reshape(cell2mat({ax.Position}'),numel(ax),4);
axPos = [abs(min(diff(allAxPos(:,1)))) max(diff(allAxPos(:,2))) allAxPos(1,3) allAxPos(1,4)];
absPos = [(axPos(1:2)/2).*h.Position(3:4) axPos(3:4).*h.Position(3:4)];
newFigPos = [0 1 absPos(1)+absPos(3) absPos(2)+absPos(4)];
if ~isempty(lgnd)
    
    allLgndPos=reshape(cell2mat({lgnd.Position}'),numel(lgnd),4);
    relPos = allLgndPos-allAxPos;
    absLgndPos = [absPos(1:2)/2+relPos(1,1:2).*h.Position(3:4) allLgndPos(1,3:4).*h.Position(3:4)];
end
%axPos = ax(1).Position;
for iAx = 1:numel(ax)
    hNew = figure;
    if ~isempty(lgnd)
        [axNew] = copyobj([ax(iAx) lgnd(iAx)],hNew);
        lgndNew = axNew(2);
        axNew=axNew(1);
        set(lgndNew,'Position',[absLgndPos(1:2)./newFigPos(3:4) absLgndPos(3:4)./newFigPos(3:4)]);
    else
        axNew = copyobj(ax(iAx),hNew);
    end
    set(hNew,'Position',newFigPos)
    set(axNew,'Position',[(absPos(1:2)/2)./newFigPos(3:4) absPos(3:4)./newFigPos(3:4)])
    if strcmp(saveType,'png')
        savepng(hNew,sprintf('%s_%d.%s',baseFileName,iAx,saveType),[1 1 1],res)
    elseif strcmp(saveType,'pdf')
        savepdf(hNew,sprintf('%s_%d.%s',baseFileName,iAx,saveType))
    elseif strcmp(saveType,'both')
        savepng(hNew,sprintf('%s_%d.%s',baseFileName,iAx,'png'),[1 1 1],res)
        savepdf(hNew,sprintf('%s_%d.%s',baseFileName,iAx,'pdf'))
    end
    close(hNew)
end
end

