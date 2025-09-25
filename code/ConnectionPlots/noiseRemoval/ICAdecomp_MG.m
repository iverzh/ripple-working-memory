function ICAdecomp_MG(icaMethod,subjIdcs,trlSlct,filtDat,avgRef)

%% Script to generate a fieldtrip structure for subjects that do not have it
%relies on analysis files already having been generated, and the trigger
%data from loadTrialInfo_v2

%v2 created for the expanded pipeline w full subject pool

if ~exist('subjIdcs','var')||isempty(subjIdcs)
    subjIdcs=190
end

if nargin<3
    trlSlct=false
end
if trlSlct
    trlSlctApnd='_trlSlct';
else
    trlSlctApnd='';
end
if nargin<4
    filtDat=false
end
if filtDat
    filtApnd='_filtDat';
else
    filtApnd='';
end

if nargin<5
    avgRef=false;
end

if avgRef
    apndAvgRef='_avgRef';
else
    apndAvgRef='';
end

addpath(genpath('/home/jgarret/fieldtrip-20210529'))

%DONE: 70
%subjIdcs=[-240 -89 442 321 313 311 297 293 259 255 233 226 190 170 133 130 ...
%   111 101 90 87 86 77 71 68 67 50]; %59 has messed up triggers
site = 'MG';

if filtDat
    datRoot='/space/seh8/5/halgdev/projects/jgarret/altnalysis_1kHz_single/';
else
    datRoot='/space/seh8/5/halgdev/projects/jgarret/altnalysis_native_fs/';
end

[subjects,site]=genSubjList(subjIdcs,site);

%clear data cfg
for s=1:numel(subjects)
    subj=subjects{s}
    mobj=matfile(sprintf('%s/%s/%s/tasks/%s_WN_preICA%s.mat',datRoot,site,subj,subj,apndAvgRef));
    
%     data.label=mobj.label';
%     if strcmp(subj,'NY233')
%         for n=205:224
%             data.label{n}=[data.label{n}(1:4) '2-' data.label{n}(5:end)];
%         end
%     end
    
    
    %emptyChans=strcmpi(data.label,'Empty')|strcmp(data.label,'a')|strcmp(data.label,'-');
    

    
%     [labParcs,sz]=retrieveChanParcels(subj);
    
%     if strcmp(subj,'NY233')
%         labParcs(205:224,4)=strrep(labParcs(205:224,4),'EEG','EEGG');
%     elseif strcmp(subj,'NY170')
%         labParcs(125:126,4)={'X01','XX01'};
%     end
    
%     emptyChans=cellfun(@(x) ~str2num(x(regexp(x,'[0-9]')))&&~any(strcmp(x(regexp(x,'[A-Z]')),{'LC','RC'}),2),labParcs(:,4));
%     if any(strcmp(subj,{'NY255','NY190'})) %not all at end, but they're all different names so not a problem
%         emptyChans=false(size(labParcs,1),1)
%     end 
    

    dat=mobj.data;
    if trlSlct
        dat=dat-median(dat,1);
        trl=loadTrials(subj,'WN');
        trlMsk=bounds2mask(trl(:,1)+[-500 2500],length(dat));
        dat(~trlMsk,:)=0;
    end
%     if strcmp(subj,'NY190')
%         deadChans=std(dat)<20;
%     else%if strcmp(subj,'NY133')
%         deadChans=std(dat)<2;
%     end

    subjIdx=str2num(subj(regexp(subj,'[0-9]')));
    deadThresh=2;
    if any(subjIdx==[311 297 255 2402 233 226 190 130 892])
        deadThresh=deadThresh*10;
    end
    if subjIdx==255
        deadThresh=deadThresh*100;
    end
    deadChans=std(dat)<deadThresh;
    dat=dat(:,~deadChans);
    

%         data.label=labParcs(~deadChans,4);
    labAll=mobj.label;
    data.label=labAll(~deadChans);
    %data.label=data.label(~emptyChans);
    
    data.fsample=mobj.Fs;
    data.trial={dat'};
    data.time={1:size(data.trial{1},2)};
    
    
%     data.label=arrayfun(@(x) [labParcs{x,1} '-' labParcs{x,4}],1:size(labParcs,1),'uni',0)';
%     %data.label=data.label(~emptyChans);
%     data.label=data.label(~deadChans);
    
%     data.fsample=mobj.Fs;
%     data.trial={dat'};
%     data.trial={data.trial{1}(~deadChans,:)};
    
    
    
    cfg=[];
    cfg.method=icaMethod;
    if strcmp(icaMethod,'fastica')
        cfg.fastica.stabilization='on';
        cfg.fastica.approach='symm';
    else
        cfg.(icaMethod).extended=1;
    end
    comp = ft_componentanalysis(cfg,data);
    
    
%     save(sprintf('/space/seh8/5/halgdev/projects/jgarret/FW_data/ica/%s_ICA_%s.mat',subj,icaMethod),'comp','labParcs','deadChans','data','-v7.3','-nocompression')
    save(sprintf('/space/seh8/6/halgdev/projects/jacob/auditory_ripples/WN/ica/%s_ICA_%s%s%s.mat',subj,icaMethod,trlSlctApnd,filtApnd),'comp','labAll','deadChans','data','-v7.3','-nocompression')
end

end