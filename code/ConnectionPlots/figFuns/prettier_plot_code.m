% Burke's example code to prettify plots

%% basic plot
rng(0);
x = [1:100]';
y = -x+50*randn(100,1);

figure(1);clf
sH = scatter(x,y);hold on;
[F,gof] = fit(x,y,'poly1');%% basic plot
rng(0);
x = [1:100]';
y = -x+50*randn(100,1);

figure(1);clf
sH = scatter(x,y);hold on;
[F,gof] = fit(x,y,'poly1');
lH = plot(xlim,F(xlim));
xlabel('x-var')
ylabel('y-var')
title(sprintf('r^2 = %.2G',gof.rsquare))

%% pretty plot
% make plotted objects pop
sH.MarkerFaceColor = sH.MarkerEdgeColor ;
lH.LineWidth = 2;

% figure and axes settings
set(gcf,'color','w');% make figure white
axis('square');box('off')
set(gca,'linewidth',1.5,'FontSize',15);

% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');


% wysiwyg export
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved

 print(1,'./prettyfig.png','-dpng','-r200'); % -r###  is the resolution of the image in ppi


lH = plot(xlim,F(xlim));
xlabel('x-var')
ylabel('y-var')
title(sprintf('r^2 = %.2G',gof.rsquare))

%% pretty plot
% make plotted objects pop
sH.MarkerFaceColor = sH.MarkerEdgeColor ;
lH.LineWidth = 2;

% figure and axes settings
set(gcf,'color','w');% make figure white
axis('square');box('off')
set(gca,'linewidth',1.5,'FontSize',15);

% Turn on Hardware-acceleration features (if possible)
set(findall(0,'-property','AlignVertexCenters'),'AlignVertexCenters','On');
set(findall(0,'-property','GraphicsSmoothing'),'GraphicsSmoothing','On');


% wysiwyg export
setpaper; % /home/bqrosen/matlab/setpaper.m
          % makes the papersize equal to the image size
          % also turns InvertHardCopy off, so the image background color is preserved

 print(1,'./prettyfig.png','-dpng','-r200'); % -r###  is the resolution of the image in ppi

