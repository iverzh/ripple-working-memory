function fuseMG51dat
noiseNotches_MG(51,1)
mMobj=matfile('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/altnalysis_1kHz_single/MG/MG51/tasks/MG51_WN.mat');
mMobj.Properties.Writable=true;
nSm=size(mMobj,'data',1);
uMobj=matfile('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/altnalysis_1kHz_single/MG/MG51v2/tasks/MG51v2_WN.mat');
mTrl=loadTrials('MG51','WN');
uTrl=loadTrials('MG51v2','WN');
smpCorr=-6;
polCorr=-1;
baseOffset=round(uTrl(1,1)-mTrl(1,1));

oldDat=uMobj.data;
uDat=oldDat((1:nSm)+baseOffset,:);

for t=1:1000 %align
    newIdx = (mTrl(t,1)-200):nSm;
    oldIdx = (uTrl(t,1)-200)+smpCorr + (0:numel(newIdx)-1);
    uDat(newIdx,:)=oldDat(oldIdx,:)*polCorr;
end

% save
newLab=arrayfun(@(x) sprintf('UG%d',x),1:28,'uni',0)';
mMobj.label=[mMobj.label; newLab];
mMobj.data=[mMobj.data uDat];


%% localization
mMobj2=matfile('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/matlab/MG/parcs/MG51_lh_iEEG_Destrieux_v3_split.mat');
mMobj2.Properties.Writable=true;
uMobj2=matfile('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/matlab/MG/parcs/MG51v2_lh_iEEG_Destrieux_v3_split.mat');

MTstr={'G_temporal_middle'};
STstr={'G_temp_sup-Lateral'};
uLocStr=repmat(STstr,[28 1]);
uLocStr([1:8 10])=MTstr;

mMobj2.eegParcs=[mMobj2.eegParcs; [uLocStr newLab]];
mMobj2.siteCoords=[mMobj2.siteCoords; uMobj2.siteCoords];

end

