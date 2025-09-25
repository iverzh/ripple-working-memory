function Movie_genCoripDist(subjIdcs,eveAll,runStim,runResp,zeroPh,ctOpt,HGburst,bilat)


if nargin==0
    subjIdcs=[49 51];
end

try
parpool(16)
catch
    delete(gcp('nocreate'))
    parpool(16)
end



if nargin<2
    eveAll=[1];
end

if nargin<3
    runStim=1;
end
if nargin<4
    runResp=0;
end

if nargin<5
    zeroPh=0;
end

if nargin<6
    ctOpt=3;
end

if nargin<7
    HGburst=0;
end

if nargin<8
    bilat=0;
end

if zeroPh==1
    zeroPhApnd='_zeroPh';
elseif zeroPh==2
    zeroPhApnd='_zeroPhProp';
elseif zeroPh==0
    zeroPhApnd='';
end

if ctOpt==1
    ctApnd=''; %count co-ripple centers
elseif ctOpt==2
    ctApnd='_onset'; %count co-ripple starts
elseif ctOpt==3
    ctApnd='_dens'; %count co-ripple durations (density)
end

if HGburst
    eveAll=2;
end

if bilat
    rhApnd='_bilat';
else
    rhApnd='';
end

coEveAll=1;
subjects =  {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};

fRoot='/space/seh10/6/halgdev/projects/iverzh/ripples/code/bmovie-release-NWB-BIDS/code/ConnectionPlots/out';
% if ~exist([fRoot '/coripDist/'],'dir');unix(sprintf('mkdir -p %s/coripDist/',fRoot));end
opt={'-v7.3','-nocompression'};

distApnd={'', '_ripHG','_ripNoHG';...
    '_HGburst','_HGrip','_HGnoRip';...
    };
[nEve,nCoEve]=size(distApnd);

nBoot=1000;


for eve=eveAll
for coEve=coEveAll

apnd=distApnd{eve,coEve};

for s=2:numel(subjects)
    subj=subjects{s};
    tic;disp('loading trial data...')
%load co-rippling data
fprintf('Loading trial data for event %d of %d (%s)...', coEve + (eve-1)*nCoEve, nEve*nCoEve, apnd);tic
%%
%
if runStim
    tic
load(sprintf('%s/coripChanDat_%s%s%s%s%s.mat',fRoot,subj,rhApnd,apnd,zeroPhApnd,ctApnd),'coripCond');toc
outFil=sprintf('%s/coripDist_%s%s%s%s%s.mat',fRoot,subj,rhApnd,apnd,zeroPhApnd,ctApnd);

% tests=[1 2; 2 4; 3 4; 1 6; 2 6; 3 6; 4 6];
tests=[1 2; 3 4; 5 6; 5 7];% 1 4; 2 4; 3 4];

% tests=[tests; flip(tests,2)];
nTest=size(tests,1);

% subjects=fieldnames(coripDat);
nCh=size(coripCond{1},1);
nBin=size(coripCond{1},3);
if zeroPh==2
    coripDistSubj=zeros(nCh,nCh,nBin,nBoot,nTest,2,'single');
else
    coripDistSubj=zeros(nCh,nCh,nBin,nBoot,nTest,2,'uint16');
end
nTrlSubj=cellfun(@(x) size(x,4), coripCond);
nTrlSubj(cellfun(@(x) isempty(x), coripCond)) = 0;
% parpool
rng('default');
% rng(1);
for t=1:nTest
    rng(1)
    if zeroPh==2
        trlSamp=cat(4,coripCond{tests(t,:)});
    else
        trlSamp=uint8(cat(4,coripCond{tests(t,:)}));
    end
    
    nPoolSamp=sum(nTrlSubj(tests(t,:)));
    nTrl=nTrlSubj(tests(t,:));
    fprintf('stim (1/2) | %s (event %d/%d} | %s | test: %d/%d \n',...
        apnd, coEve + (eve-1)*nCoEve, nEve*nCoEve,...
        subj,... s, nSubj,...
        t, nTest)
    parfor b=1:nBoot %randomize trial identity nBoot times
%         fprintf('%d\n',b)
        condTrls=randperm(nPoolSamp);
        
        condIdx=false(nPoolSamp,1);
        condIdx(condTrls(1:nTrl(1)))=true;
        subSamp1=sum(trlSamp(:,:,:,condIdx),4);
        
        condIdx=false(nPoolSamp,1);
        condIdx(condTrls(nTrl(1)+1:end))=true;
        if sum(condIdx)~=nTrl(2)
            error('sanity check fail on trial count for bootstrapping condition 2')
        end
        subSamp2=sum(trlSamp(:,:,:,condIdx),4);
        
        coripDistSubj(:,:,:,b,t,:)=cat(6,subSamp1,subSamp2);
    end
end
nTrlSubj=cellfun(@(x) size(x,4), coripCond);
clear subjDat trlSamp trlSampTmp subSamp

toc
fprintf('Saving stimulus-locked bootstrapped testing distributions for event %d of %d (%s)...', coEve + (eve-1)*nCoEve, nEve*nCoEve, apnd);tic
save(outFil,'coripDistSubj','nTrlSubj','tests',opt{:});toc

end
%}


%%
%order sample co-ripple counts

%
if runResp
tic
load(sprintf('%s/%s/coripChanDat_resp_%s%s%s%s%s.mat',fRoot,subj,subj,rhApnd,apnd,zeroPhApnd,ctApnd),'coripCond');toc
outFil=sprintf('%s/%s/coripDist_resp_%s%s%s%s%s.mat',fRoot,subj,subj,rhApnd,apnd,zeroPhApnd,ctApnd);
condNames={'TP'};
% conds=[1 2 4 3 9 10 7];
nCond=numel(coripCond);
% tests=[1 2; 3 2; 3 1; 1 4; 2 4; 5 6; 3 6; 2 6];
% tests=[6 2; 3 6; 2 4; 6 5; 6 4; 3 2; 3 7; 1 2; 1 4];%1 2; 3 2; 3 1; 1 4; 2 4; 5 6;  2 6];
% tests=[1 2; 1 3; 1 5; 4 1; 4 7; 6 2; 6 3; 2 3];
tests=[1 2; 1 3; 1 4];

% tests=[tests; flip(tests,2)];
nTest=size(tests,1);

% subjects=fieldnames(coripDat);
nCh=size(coripCond{1},1);
nBin=size(coripCond{1},3);
coripDistSubj=zeros(nCh,nCh,nBin,nBoot,nTest,2,'uint32');
nTrlSubj=cellfun(@(x) (~~size(x,1))*size(x,4), coripCond);
% parpool

rng('default');
% rng(1);
for t=1:nTest
%     rng(t)
%     trlSamp=uint8(cat(4,coripCond{tests(t,:)}));
    trlSamp=cat(4,coripCond{tests(t,:)});
    
    nPoolSamp=sum(nTrlSubj(tests(t,:)));
    nTrl=nTrlSubj(tests(t,:));
    fprintf('resp (2/2) | %s (event %d/%d} | %s | test: %d/%d \n',...
        apnd, coEve + (eve-1)*nCoEve, nEve*nCoEve,...
        subj,... s, nSubj,...
        t, nTest)
    if any(~nTrl);continue;end
    parfor b=1:nBoot %randomize trial identity nBoot times
        condTrls=randperm(nPoolSamp);
        
        condIdx=false(nPoolSamp,1);
        condIdx(condTrls(1:nTrl(1)))=true;
        subSamp1=sum(trlSamp(:,:,:,condIdx),4);
        
        condIdx=false(nPoolSamp,1);
        condIdx(condTrls(nTrl(1)+1:end))=true;
        if sum(condIdx)~=nTrl(2)
            error('sanity check fail on trial count for bootstrapping condition 2')
        end
        subSamp2=sum(trlSamp(:,:,:,condIdx),4);
        
        coripDistSubj(:,:,:,b,t,:)=cat(6,subSamp1,subSamp2);
    end
end
nTrlSubj=cellfun(@(x) size(x,4), coripCond);
clear subjDat trlSamp trlSampTmp subSamp

toc
fprintf('Saving stimulus-locked bootstrapped testing distributions for event %d of %d (%s)...', coEve + (eve-1)*nCoEve, nEve*nCoEve, apnd);tic
save(outFil,'coripDistSubj','nTrlSubj','tests',opt{:});toc
end

end    
end
end