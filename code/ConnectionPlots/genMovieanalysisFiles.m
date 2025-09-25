function genWNanalysisFiles
%%
addpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/MULTI')
% addpath('/home/jgarret/mat')
    subjFils=readcell('/home/jgarret/TaskAnalysis/lang/WN/WNsubjectRawData.xlsx');


subjects=subjFils(:,1);%subjFils(:,1)=[];
subjects = {'MG49'};

outRoot='/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/MG49/data_1kHz/';

for s=1 %:numel(subjects)
    subj=subjects{s};
    if strcmp(subj,'MG51v2')
        lowPass=1;
    else
        lowPass=0;
    end
    if strcmp(subj,'MG51')
        runNotch=1;
    else
        runNotch=1;
    end
    %if ~strcmp(subj,'NY255');continue;end
    %if strcmp(subj,'NY255');warning('skipping NY255 due to error, fix this');continue;end
    nFils=0;
    if strcmp(subj(3),'0');subj(3)=[];end
    outDir= outRoot; %[outRoot subj '/tasks/'];
    if ~exist(outDir,'dir')
        unix(sprintf('mkdir -p %s',outDir))
    end
    outFil=[outDir subj '_SA_v2.mat'];
    if exist(outFil,'file')
        warning('%s file exists, skipping...',subj)
        continue
    end
    clear partStrct
    for f=1%:size(subjFils,2)-1
%         currInFil='/space/seh8/6/halgdev/projects/jacob/auditory_ripples/raw_data_backup/MG49/MG49_SA/MG4909SAtele.set';
        currInFil='/space/seh8/6/halgdev/projects/jacob/auditory_ripples/raw_data_backup/MG49_1/MG4909SAtele.edf';

        currOutFil=sprintf('%s_part%d.mat',outFil(1:end-4),f);
        if ~ismissing(currInFil)
            nFils=f;
            if strcmp(currInFil(end-3:end),'.set')
                MULTI_set2mat_task_highPrecision(currInFil,currOutFil,[],1,1,lowPass);
            elseif strcmp(currInFil(end-3:end),'.dat')
                MULTI_dat2mat_task_highPrecision(currInFil,currOutFil);
            elseif strcmp(currInFil(end-3:end),'.edf')
                MULTI_edf2mat_task_highPrecision(currInFil,currOutFil,[],runNotch);
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
