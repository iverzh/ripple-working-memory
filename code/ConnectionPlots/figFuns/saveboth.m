function saveboth(fig,path,clr)
%SAVEBOTH Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    clr=[1 1 1];
end
    savepng(fig,[path '.png'], clr);
    savepdf(fig,[path '.pdf']);
end

