function [msk,bnds] = getVHFnoisePeriods(subjIdx,site,task,thresh,pad)
if numel(subjIdx)>1
    error('one subject at a time')
end
subj=genSubjList(subjIdx,site);
subj=subj{1};
if nargin<4
    thresh=4;
end
if nargin<5
    pad=0.4;% seconds
end

load(sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/altnalysis_1kHz_single/%s/%s/tasks/%s_%s.mat',site,subj,subj,task),...
    'data','Fs')
% load(sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/ripples/%s/%s_ripple_stats_%s_0_cln.mat',site,subj,subj,task),...
%     'rippleStats')
[b,a]=butter(3,200/(Fs/2),'high');
VHFamp=abs(hilbert(filtfilt(b,a,double(data))));
VHFscore=zscore(smoothdata(median(VHFamp,2),1,'gaussian',200),1);
bnds=mask2bounds(VHFscore>thresh);
msk=bounds2mask(bnds,numel(VHFscore),pad*Fs);
end

