
    load('/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/MG51/tasks/MG51_WN_raw.mat', 'data')
    trl=loadTrials('MG51','WN');
    trlMsk=bounds2mask(trl(:,1)+[-1000 2000],length(data));
    data=data-median(data,1);
%     data(~trlMsk,:)=0;
    save('/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/MG51/tasks/MG51_WN_preICA.mat', 'data','-append')
    