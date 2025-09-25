function savepng(h,filename,clr,res)

if nargin<4
    res=0;
end
if nargin<3
    clr='w';
end
% figure and axes settings
set(h,'color',clr);% make figure white
ax = findall(h,'type','axes');
for iAx=1:numel(ax)
%axis(ax, 'square');
box(ax(iAx),'off')
set(ax(iAx),'linewidth',1,'FontSize',6);
end

set(h,'PaperUnits','inches')
PP = get(h,'paperposition');
PP(1:2) = 0.5;
set(h,'papersize',PP(3:4)+.5);
set(h,'InvertHardCopy','off') % preserve background color
print(h,filename,'-dpng',sprintf('-r%d',res)); % -r###  is the resolution of the image in ppi


end