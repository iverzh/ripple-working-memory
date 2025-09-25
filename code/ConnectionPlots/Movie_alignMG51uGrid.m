load('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/coripDat/MG51/phRip_MG51.mat', 'phRip','labParc')
phRip=cos(phRip)+1i*sin(phRip);

load('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/altnalysis_1kHz_single/MG/MG51/tasks/MG51_WN.mat', 'data','label')

%%

dat=phRip;
% dat=data;
trl=loadTrials('MG51','WN');
dat=dat(trl(1,1):trl(end,1)+2000,:);

n=200;

%phRip
% uCh=113;
% mCh=45;
% uCh=109:117;
% mCh=43:45;

uCh=105:126;
mCh=[45 46 53 54];

%full data
% uCh=135:143;
% mCh=45;

% uCh=125:152;
% mCh=[45 46 53 54];

uDat=nanmean(dat(:,uCh),2);
mDat=nanmean(dat(:,mCh),2);

covVec=cellfun(@(y) abs(y(2)),...
    arrayfun(@(x) cov(mDat(1+n+x:end-n+x),uDat(1+n-x:end-n-x),'omitrows'),...
    -n:n,'UniformOutput',false));

figure();plot(-n*2:2:n*2,covVec)
