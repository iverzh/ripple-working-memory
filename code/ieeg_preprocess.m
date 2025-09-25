
close all
clc
clear



%% INPUTS

dataDirectory =     '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/preprocess';
exportDirectory = fullfile(dataDirectory, '../out');
if ~isfolder(exportDirectory); mkdir(exportDirectory); end

processMacro = false;
processMicro = true;

%% macro contacts

files = dir(fullfile(dataDirectory, '*macro*'));
files = {files.name}';


switch processMacro
    case true

        figure; 
        for f = 1:length(files)
            fprintf('... loading %s ... \n', files{f})
            mObj = load(fullfile(dataDirectory, files{f}));
            lab = cellstr(mObj.chan_labels);
            dat = mObj.lfp_data;
            subject = files{f}(1:8);
        
        
        
        
        
            nL = numel(lab);
            [alph,numr] = deal(cell(nL,1));
            for iL = 1:nL
              nMsk = isstrprop(lab{iL},'digit');
              if ~nMsk(end) % if label does not end in a number
                alph{iL} = lab{iL};
                numr{iL} = NaN;
              else% if label ends in a number
                % consider everything before the number at the end the name aka alph
                % (no matter how many digits that number is)
                idx = find(~nMsk,1,'last');
                alph{iL} = lab{iL}(1:idx);
                numr{iL} = str2double(lab{iL}(idx+1:end));
              end
            end
        
            bpLabel = cell(nL-1,1);
            bpCoords = nan(nL-1, 3);
            dIdx = false(nL-1,1);
            for iL = 1:nL-1
              if ~isequal(alph(iL),alph(iL+1))
                dIdx(iL) = true;
              else
                bpLabel{iL} = sprintf('%s%.2i-%.2i',...
                  alph{iL},numr{iL+1},numr{iL});
                bpCoords(iL, :) = mean([mObj.chan_coords(iL,:); mObj.chan_coords(iL+1,:)]);
              end
            end
            bpLabel(cellfun(@(X) isempty(X), bpLabel)) = [];
            bpCoords(dIdx, :) = [];
            seegMsk = true(1, nL);
            
        
        
            bpV = double(seegMsk);
            bpOp = diag(bpV(1:end-1),-1) + diag(-bpV); 
            bpOp(:,~seegMsk) = [];
            bpOp(:,[find(~arrayfun(@isequal,alph(1:end-1),alph(2:end)));size(bpOp,2)]) = [];%remove intershaft derived channels
        
            data = 1e6 * dat'*bpOp; %bipolarize and convert from V to microV
            for ch = 1:size(data,2)
                [b,a] = butter(3,0.1/(mObj.fs/2), "high");
                data(:,ch) = filtfilt(b,a,data(:,ch));
        
                for nf = 60:60:mObj.fs/2 %240
                    Wo = nf/(mObj.fs/2);
                    BW = Wo/35;
                    [b,a] = iirnotch(Wo, BW);
                    data(:,ch) = filtfilt(b,a,data(:,ch));
                end
        
            end
        
        
            chan_labels = bpLabel;
            fs = mObj.fs;
            filename = sprintf('%s_LFP_macro_bp.mat', subject);
            fprintf('... saving  %s ... \n', filename)
            save(fullfile(exportDirectory, filename), 'data', 'chan_labels', 'fs', 'subject', 'bpCoords', '-v7.3')
        
            pl = plot3(mObj.chan_coords(:,1), mObj.chan_coords(:,2), mObj.chan_coords(:,3), '.'); hold on;
            pl.Color = [0.7 0.7 0.7];
        
            pl = plot3(bpCoords(:,1), bpCoords(:,2), bpCoords(:,3), 'r.'); hold on;
        
        
          
        
        end
        % dataBp = lfp' * bpOp;
        
        axis equal
end

%% micro data 

switch processMicro
    case true


        LFPfiles = dir(fullfile(dataDirectory, '*micro*'));
        LFPfiles = {LFPfiles.name}';
        
        unitfiles = dir(fullfile(dataDirectory, '*unit*'));
        unitfiles = {unitfiles.name}';

        
        
        for f = 1:length(files)
            fprintf('... loading %s ... \n', LFPfiles{f})
            lObj = load(fullfile(dataDirectory, LFPfiles{f}));
            uObj = load(fullfile(dataDirectory, unitfiles{f}));
            
            lab = cellstr(lObj.chan_labels);
            dat = lObj.lfp_data;
            subject = LFPfiles{f}(1:8);
            uElec = double(uObj.electrode_id);

            units = LoadSpikeTimes(subject, 'RutishauserLab');
        
            nL = numel(lab);
            [alph,numr] = deal(cell(nL,1));
            for iL = 1:nL
              nMsk = isstrprop(lab{iL},'digit');
              if ~nMsk(end) % if label does not end in a number
                alph{iL} = lab{iL};
                numr{iL} = NaN;
              else% if label ends in a number
                % consider everything before the number at the end the name aka alph
                % (no matter how many digits that number is)
                idx = find(~nMsk,1,'last');
                alph{iL} = lab{iL}(1:idx);
                numr{iL} = str2double(lab{iL}(idx+1:end));
              end
            end
        
        
        
        %     bundleNames = unique(alph);
        %     for iB = 1:length(bundleNames)
        %         iiB = find(strcmp(alph, bundleNames{iB}));
        %         ref = find(~ismember(iiB, uElec), 1, 'first'); % find electrode without unit ativity 
        % 
        %         data
        % 
        % 
        % 
        %     end
        
            data = dat'; %bipolarize and convert from V to microV
            for ch = 1:size(data,2)
                [b,a] = butter(3,0.1/(lObj.fs/2), "high");
                data(:,ch) = filtfilt(b,a,data(:,ch));
        
                for nf = 60:60:lObj.fs/2 %240
                    Wo = nf/(lObj.fs/2);
                    BW = Wo/35;
                    [b,a] = iirnotch(Wo, BW);
                    data(:,ch) = filtfilt(b,a,data(:,ch));
                end
        
            end
        
        
            chan_labels = lab;
            fs = lObj.fs;
            chan_coords = lObj.chan_coords;
            filename = sprintf('%s_LFP_micro.mat', subject);
            fprintf('... saving  %s ... \n', filename)
            save(fullfile(exportDirectory, filename), 'data', 'chan_labels', 'fs', 'subject', 'chan_coords', '-v7.3')
        
        
        end

end





































