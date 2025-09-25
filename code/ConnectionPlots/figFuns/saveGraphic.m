function saveGraphic(fig,path,name,res)
%SAVEBOTH Summary of this function goes here
%   Detailed explanation goes here

if nargin<4||isempty(res)
    res='-r0';
    resPDF=[];
else
    res=sprintf('-r%d',round(res));
    resPDF=sprintf('-r%d',res);
end
if nargin<3||isempty(name)
    if strcmp(path(end),'/')
        error('must specify file name')
    else
        slashes=regexp(path,'/');
        name=path(slashes(end)+1:end);
        path=path(1:slashes(end));
    end
end
if ~exist(path,'dir');unix(sprintf('mkdir -p %s',path));end
if ~exist([path '/png/'],'dir');mkdir([path '/png/']);end
if ~exist([path '/pdf/'],'dir');mkdir([path '/pdf/']);end
print(fig,[path '/png/' name '.png'],'-dpng',res);
if ~isempty(resPDF)
    savepdf(fig,[path '/pdf/' name '.pdf'],resPDF);
else
    savepdf(fig,[path '/pdf/' name '.pdf']);
end
    
end

