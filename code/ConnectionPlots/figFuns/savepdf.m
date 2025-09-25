function savepdf(h, save_path,res)
% saves PDFs that don't suck
if nargin<3||isempty(res)
    res=[];
end
set(h, 'Units', 'Inches');
pos = get(h, 'Position');
set(h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);
% if exist(save_path,'file')
%     unix(sprintf('rm %s',save_path));
% end
% if isnumeric(res)&&~isempty(res)
%     res=round(res);
%     res=sprintf('-r%d',res);
%     print(h, save_path, '-dpdf', res)
% elseif ischar(res)
%     print(h,save_path,'-dpdf',res)
% elseif isempty(res)
%     set(0,'defaultFigureRenderer','painters')
    print(h,save_path,'-dpdf')
% else
%     error('resolution option not recognized')
% end
    
end

