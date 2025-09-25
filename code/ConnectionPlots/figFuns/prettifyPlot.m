function prettifyPlot(h,ax)
%% pretty plot
addpath('/home/bqrosen/matlab')
% figure and axes settings
set(h,'color','w');% make figure white
%axis('square')
box('off')
set(ax,'linewidth',1.5,'FontSize',15);

% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');


% wysiwyg export
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved

end

