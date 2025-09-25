function WN_save_rip_phase(subjIdcs)


if nargin==0
    subjIdcs=[49 51];
end
subjects=genSubjList(subjIdcs,'MG');

HGapnd='';
for s=1:numel(subjects)
    subj=subjects{s}
    outDir=sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/coripDat/%s/',subj);
    if ~exist(outDir,'dir');unix(sprintf('mkdir -p %s',outDir));end
    
    tic;fprintf('loading data...')
    load(sprintf('/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/MG/%s/tasks/%s_WN.mat',subj,subj),'data','label','Fs')
    nS=size(data,1);
%     load(sprintf('/space/seh8/6/halgdev/projects/jacob/NC_ripple/matFiles/MG/%s_ripple_stats_WN_0_cln%s.mat',subj,HGapnd),'rippleStats')
    load(sprintf('/space/seh8/6/halgdev/projects/jacob/NC_ripple/matFiles/MG/%s_ripple_stats_WN_0_trlSlct%s.mat',subj,HGapnd),'rippleStats')
    
    [labParc,~,~,~,~,rejChan,coords]=retrieveChanParcels_MG(subj);
    toc
    
    tic;fprintf('channel bookkeeping...')
    
    rmCh =[rejChan; true(size(labParc,1)-numel(rejChan),1)]|strcmp(labParc(:,1),'N/A');
    labParc(rmCh,:)=[];
    data(:,rmCh)=[];
%% analysis


rip_bounds_mask = false(size(data));
gCh=find(~rmCh)';

[b,a]=butter(3,[70 110]/(Fs/2),'bandpass');
phRip=angle(hilbert(filtfilt(b,a,double(data))));


tic;fprintf('populating ripples...')
for ch = gCh
    %     rip_mask(ch,rippleStats.locs{ch}) = true;
    for rip=1:size(rippleStats.window{ch},1)
        rip_bounds_mask(rippleStats.window{ch}(rip,1):rippleStats.window{ch}(rip,2),ch==gCh) = true;
%         rip_peak_mask(rippleStats.locs{ch}(rip),ch==gCh) = true;
    end
end

phRip(~rip_bounds_mask)=nan;


tic;fprintf('saving data...')
save(sprintf('%s/phRip_%s.mat',outDir,subj),...
    'gCh','labParc','phRip','subj','-v7.3','-nocompression')

end
end