function genWNanalysisFiles_nativeFs
%%
addpath('~/MULTI')
addpath('~/mat')
    subjFils=readcell('/home/jgarret/TaskAnalysis/lang/WN/WNsubjectRawData.xlsx');

resamp=0;
runNotch=0; % the purpose of these files is careful characterization of noise components, so leave the data untouched for now
subjects=subjFils(:,1);%subjFils(:,1)=[];


outRoot='/space/seh8/5/halgdev/projects/jgarret/altnalysis_native_fs/MG/';

for s=1:numel(subjects)
    subj=subjects{s};
%     if strcmp(subj,'MG51')
%         runNotch=1;
%     else
%         runNotch=1;
%     end
    %if ~strcmp(subj,'NY255');continue;end
    %if strcmp(subj,'NY255');warning('skipping NY255 due to error, fix this');continue;end
    nFils=0;
    if strcmp(subj(3),'0');subj(3)=[];end
    outDir=[outRoot subj '/tasks/'];
    if ~exist(outDir,'dir')
        unix(sprintf('mkdir -p %s',outDir))
    end
    outFil=[outDir subj '_WN.mat'];
    if exist(outFil,'file')
        warning('%s file exists, skipping...',subj)
        continue
    end
    clear partStrct
    for f=1:size(subjFils,2)-1
        currInFil=subjFils{strcmp(subjFils(:,1),subj),f+1};
        currOutFil=sprintf('%s_part%d.mat',outFil(1:end-4),f);
        if ~ismissing(currInFil)
            nFils=f;
            if strcmp(currInFil(end-3:end),'.set')
                MULTI_set2mat_task_highPrecision(currInFil,currOutFil,[],runNotch,resamp);
            elseif strcmp(currInFil(end-3:end),'.dat')
                MULTI_dat2mat_task_highPrecision(currInFil,currOutFil);
            elseif strcmp(currInFil(end-3:end),'.edf')
                MULTI_edf2mat_task_highPrecision(currInFil,currOutFil,[],runNotch,resamp);
            elseif strcmp(currInFil(end-3:end),'.eeg')
                MULTI_eeg2mat_task_highPrecision(currInFil,currOutFil);
            elseif strcmp(currInFil(end-1:end),'.e')
                MULTI_e2mat_task_highPrecision(currInFil,currOutFil);
            end
            ctr=0;
            while ~exist(currOutFil,'file')
                ctr=ctr+1;
                if ctr>10
                    error('could not load part %d file',f)
                end
                pause(1)
            end
            partStrct(f)=load(currOutFil);
        else
            if f==1
                warning('skipping subject %s...',subj)
            end
            break
        end
    end
    
    %fuse part files
    if nFils>0
        fusedStrct = partStrct(1);

        if nFils>1
            fusedStrct.data=cell2mat({partStrct.data}');
            fusedStrct.trigData=cell2mat({partStrct.trigData}');
            fusedStrct.clockT=1:length(fusedStrct.data);
            fusedStrct.edfFiles=[partStrct.edfFiles];
            checkLabs=[partStrct.label];
            if ~all(arrayfun(@(x) all(strcmp(checkLabs(x,1),checkLabs(x,2:end))), 1:size(checkLabs,1)))
                error('not all labels match across files')
            end
        end

        save(outFil,'-struct','fusedStrct','-v7.3','-nocompression')
    end
end

end
