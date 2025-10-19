
% subjects suported: MG29


function unitsFinal = LoadSpikeTimes(subject,varargin)

if nargin > 1
    formatType = varargin{1};
else
    formatType = 'Anna';
end

if nargin > 2 && strcmp(formatType, 'RutishauserLab')
    study = varargin{2};    
elseif nargin > 2
    arrayName = varargin{2};    
end

switch formatType
    case 'Anna'
        switch subject
            case 'MG29'
                spikeDirectory = '/space/seh8/2/mdeh1-3/halgdev/projects/mmilanguage/unit_sorting/mg29_7-15/Sleep';
                spikeFoldersStruct = dir(spikeDirectory);
                spikeFolders = {spikeFoldersStruct.name};
                checkIfDir = [spikeFoldersStruct.isdir];
                spikeFolders = spikeFolders(checkIfDir); 
                
                keep = true(1,length(spikeFolders));
                for f = 1:length(keep); if any(strcmp(spikeFolders{f},{'.','..'})); keep(f) = false; end; end
                spikeFolders = spikeFolders(keep);
                
                pyr = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/MG29/MG29_pyr_unit_index.mat');
                pyrChannels = unique(pyr.unit_index(1,:));
                
                int = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Int_Units/MG29/MG29_int_unit_index.mat');
                intChannels = unique(int.unit_index(1,:));
                
                mult = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Multi_Units/MG29/MG29_multi_unit_index.mat');
                multChannels = unique(mult.unit_index(1,:));

                spin_offset = 0;

            case 'MG29wc'
                spikeDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/MG29/preprocess/waveclus/out';
                spikeFoldersStruct = dir(fullfile(spikeDirectory,'*.mat'));
                spikeFolders = {spikeFoldersStruct.name};
        
                detectedChannels = [];
                for f = 1:length(spikeFolders)
                    file = spikeFolders{f}; 
                    detectedChannels(f) = str2double(file(end-7:end-4)); 
                end
                
        %         pyr = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/MG29/MG29_pyr_unit_index.mat');
                pyrChannels = []; %unique(pyr.unit_index(1,:));
                
        %         int = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Int_Units/MG29/MG29_int_unit_index.mat');
                intChannels = []; %unique(int.unit_index(1,:));
                
        %         mult = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Multi_Units/MG29/MG29_multi_unit_index.mat');
                multChannels = unique(detectedChannels);
        
                
                spin_offset = 0;
            case 'MG49'
                spikeDirectory = '/space/seh8/2/mdeh1-3/halgdev/projects/mmilanguage/unit_sorting/mg49_sleep/Sleep-Anna/';
                spikeFoldersStruct = dir(spikeDirectory);
                spikeFolders = {spikeFoldersStruct.name};
                checkIfDir = [spikeFoldersStruct.isdir];
                spikeFolders = spikeFolders(checkIfDir); 
                
                keep = true(1,length(spikeFolders));
                for f = 1:length(keep); if any(strcmp(spikeFolders{f},{'.','..'})); keep(f) = false; end; end
                spikeFolders = spikeFolders(keep);
                
                pyr = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/MG49/MG49_pyr_unit_index.mat');
                pyrChannels = unique(pyr.unit_index(1,:));
                
                int = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Int_Units/MG49/MG49_int_unit_index.mat');
                intChannels = unique(int.unit_index(1,:));
                
                mult = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Multi_Units/MG49/MG49_multi_unit_index.mat');
                multChannels = unique(mult.unit_index(1,:));
        
                
                spin_offset = 7203900;
                
            case 'MG67'
                spikeDirectory = '/space/seh8/2/mdeh1-3/halgdev/projects/mmilanguage/unit_sorting/mg67/Sleep/';
                spikeFoldersStruct = dir(spikeDirectory);
                spikeFolders = {spikeFoldersStruct.name};
                checkIfDir = [spikeFoldersStruct.isdir];
                spikeFolders = spikeFolders(checkIfDir); 
                
                keep = true(1,length(spikeFolders));
                for f = 1:length(keep); if any(strcmp(spikeFolders{f},{'.','..'})); keep(f) = false; end; end
                spikeFolders = spikeFolders(keep);
                
                pyr = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/MG67/MG67_pyr_unit_index.mat');
                pyrChannels = unique(pyr.unit_index(1,:));
                
                int = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Int_Units/MG67/MG67_int_unit_index.mat');
                intChannels = unique(int.unit_index(1,:));
                
                mult = load('/space/mdeh4/1/halgdev/projects/sargsyan/SpindleDetector/STDP_Analysis/Multi_Units/MG67/MG67_multi_unit_index.mat');
                multChannels = unique(mult.unit_index(1,:));
        
                
                spin_offset = 0; 
        
        end

    case 'CellExplorer'
        switch subject
            case 'MG29'
                spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/preprocess/waveclus_cOV_90_batch2/out',subject);
                baseDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/CellExplorerOut/', subject);

            case 'T11'
                spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/preprocess/waveclus/%s/out',subject, arrayName);
                baseDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/CellExplorerOut/%s', subject, arrayName);
            case 'MG49'
                spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/preprocess/waveclus/%s/out',subject, arrayName);
                baseDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/CellExplorerOut/%s', subject, arrayName);

            otherwise 
                spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/preprocess/waveclus/out',subject);
                baseDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/CellExplorerOut/', subject);

        end

        spikeFoldersStruct = dir(fullfile(spikeDirectory,'*.mat'));
        spikeFolders = {spikeFoldersStruct.name};

%         detectedChannels = [];
%         for f = 1:length(spikeFolders)
%             file = spikeFolders{f}; 
%             detectedChannels(f) = str2double(file(end-7:end-4)); 
%         end
        
        load(fullfile(baseDir,"cell_metrics.mat"), 'cell_metrics')
        
        detectedChannels = unique(cell_metrics.maxWaveformCh1);
        pyrInd = strcmp(cell_metrics.putativeCellType,'Pyramidal Cell');
        pyrChannels = unique(cell_metrics.maxWaveformCh1(pyrInd));

        intInd = strcmp(cell_metrics.putativeCellType,'Narrow Interneuron') | strcmp(cell_metrics.putativeCellType,'Wide Interneuron');
        intChannels = unique(cell_metrics.maxWaveformCh1(intInd));
       
        multInd = strcmp(cell_metrics.putativeCellType,'multi') | strcmp(cell_metrics.putativeCellType,'mult')| strcmp(cell_metrics.putativeCellType,'Unknown');
        multChannels = unique(cell_metrics.maxWaveformCh1(multInd));
        SU = true(1, length(multInd));

    case 'RutishauserLab'
        spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/%s/preprocess/units', study);
        files = dir(sprintf('%s/*%s*_units.mat', spikeDirectory, subject));
        files = {files.name};
        filename = files{1};
        uObj = load(fullfile(spikeDirectory, filename));
       
        detectedChannels = unique(uObj.electrode_id);
        pyrChannels = detectedChannels;
        if strcmp(study,'bmovie')
            channelLocations = cellfun(@(X) X{6}, uObj.electrodes, 'UniformOutput', false);
            unitAmplitudes = cellfun(@(X) abs(X(:, 95)), uObj.waveforms);
%             SU = true(1, length(unitAmplitudes));
            percentRefract = cellfun(@(X) sum(diff(X) < 3*1e-3) / length(X) * 100, uObj.spike_times);
            ISI = cellfun(@(X) diff(X)*1e3, uObj.spike_times, 'UniformOutput', false);
            N = cellfun(@(X) histcounts(X, -0.5:1:100.5), ISI, 'UniformOutput', false);
            burstIndex = cellfun(@(X) (mean(X(1:11))-mean(X(41:51)))/mean(X(41:51)), N);
            SU =  percentRefract < 3;
            recov = cellfun(@(X) diff(smoothdata(X(:, 105:end), 'gaussian', 10)), uObj.waveforms, 'UniformOutput', false);
            zc = cellfun(@(X) find(X(1:end-1) .* X(2:end) < 0, 1, 'first'), recov, 'UniformOutput', false);
            iZ = cellfun(@(X) isempty(X), zc); 
            if sum(iZ) > 0; zc(iZ) = {nan}; end
            troughToPeak = cell2mat(zc) + 10;

        elseif strcmp(study,'Sternberg')
            channelLocations = cellfun(@(X) X{4}, uObj.electrodes, 'UniformOutput', false);
            unitSNR = cellfun(@(X) [max(mean(squeeze(X)))-min(mean(squeeze(X)))]/(2*mean(std(squeeze(X)))), uObj.waveforms);
            percentRefract = cellfun(@(X) sum(diff(X) < 3*1e-3) / length(X) * 100, uObj.spike_times);
            ISI = cellfun(@(X) diff(X)*1e3, uObj.spike_times, 'UniformOutput', false);
            N = cellfun(@(X) histcounts(X, -0.5:1:100.5), ISI, 'UniformOutput', false);
            burstIndex = cellfun(@(X) (mean(X(1:11))-mean(X(41:51)))/mean(X(41:51)), N);
            nAP = cellfun(@(X) length(X), uObj.spike_times);
            for iU = 1:length(uObj.waveforms)
                if nAP(iU) < 2
                    uObj.waveforms{iU} = [uObj.waveforms{iU},uObj.waveforms{iU}];
                end
            end
            
            if contains(subject, 'JHU')
                unitAmplitudes = cellfun(@(X) abs(mean(squeeze(X(:, :, 122)))), uObj.waveforms);
                recov = cellfun(@(X) diff(smoothdata(mean(squeeze(X(:, :, 142:end))), 'gaussian', 10)), uObj.waveforms, 'UniformOutput', false);
                zc = cellfun(@(X) find(X(1:end-1) .* X(2:end) < 0, 1, 'first'), recov);
                troughToPeak = zc + 20; 

            else
                unitAmplitudes = cellfun(@(X) abs(mean(squeeze(X(:, :, 95)))), uObj.waveforms);
                recov = cellfun(@(X) diff(smoothdata(mean(squeeze(X(:, :, 105:end))), 'gaussian', 10)), uObj.waveforms, 'UniformOutput', false);
                zc = cellfun(@(X) find(X(1:end-1) .* X(2:end) < 0, 1, 'first'), recov);
                troughToPeak = zc + 10; 

            end
            
            SU = ~(unitSNR < 2.5 & percentRefract > 3);
%             SU = true(1, length(unitAmplitudes));




        end 
        intChannels = [];
        multChannels = [];
end
%%


switch formatType
    case 'RutishauserLab'
        try
                fldr =  sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/artifactSearch/%s', study);
                filename = sprintf('%s_ccgPeaks.mat', subject);
                load(fullfile( fldr, filename), 'ccg', 'Nab', 'pctOverlap')
                ccg = pctOverlap > 30 & cellfun(@(X) sum(X) > 50, Nab);
                fprintf('loading spikes for %s\n ccg matrix found\n', subject)
            catch
                ccg = false(length(uObj.electrodes));
                fprintf('loading spikes for %s\n NO ccg matrix found\n', subject)
        end
        
       G = graph(ccg, 'upper');
       % Find connected components
       [component_idx, ~] = conncomp(G);
       num_groups = max(component_idx);
       
       keepUnits = false(1,length(component_idx));
       for ii = 1:num_groups
           uGroup = find(component_idx == ii);
           uGroupAmp = unitAmplitudes(uGroup);
           keepUnits(uGroup(uGroupAmp == max(uGroupAmp))) = true;
           
       end
       
    otherwise
       keepUnits = false(1,length(SU));

       
        
        
        
        
end

nChannel =  length(detectedChannels);
% spikeFolderChannels = cellfun(@(x) str2double(x(end-7:end-4)), spikeFolders);
% cell_metrics.general.acgs.log10'
unitsFinal = cell(sum(keepUnits & SU),3);

count = 1;
for nF = 1:nChannel
      
    switch formatType
        case 'CellExplorer'
            ch = detectedChannels(nF);

            unitii = cell_metrics.maxWaveformCh1 == ch;
            types = cell_metrics.putativeCellType(unitii);
            timesAll = cell_metrics.spikes.times(unitii);
            troughToPeakAll = cell_metrics.troughToPeak(unitii);
%             burstIndexAll = cell_metrics.acg_tau_rise(unitii);
%             burstIndexAll = cell_metrics.burstIndex_Mizuseki2012(unitii);
            burstIndexAll = cell_metrics.burstIndex_Royer2012(unitii);
            rawWavAll = cell_metrics.waveforms.raw(unitii);
%             acgNarrow = cell_metrics.acg.narrow_normalized(:, unitii);
            acgNarrow = cell_metrics.acg.log10(:, unitii);
            stdWavAll = cell_metrics.waveforms.raw_std(unitii);

            nUnit = length(types);

            for unii = 1:nUnit
        
                times = timesAll{unii}; 
                times(times < 0 ) = [];
                spikeTimes = times;
                unitsFinal{count,8} = stdWavAll{unii};
                unitsFinal{count,7} = rawWavAll{unii};
                unitsFinal{count,6} = acgNarrow(:, unii)';
                unitsFinal{count,5} = burstIndexAll(unii);
                unitsFinal{count,4} = troughToPeakAll(unii);
                unitsFinal{count,2} = spikeTimes;
                unitsFinal{count,1} = ch;

                type = types{unii};
                if strcmp(type,'Pyramidal Cell')
                    unitsFinal{count,3} = 'pyr';
                elseif contains(type,'Interneuron')
                    unitsFinal{count,3} = 'int';
                elseif contains(type,'mult') 
                    unitsFinal{count,3} = 'mult';
                end
         
                count = count + 1;
                
            end

        case 'RutishauserLab'
            
                    
            ch = detectedChannels(nF);

            unitii = uObj.electrode_id == ch & keepUnits & SU;
            
            if size(ccg, 1) ~= length(unitii); error('ccg matrix not formatted correctly'); end
            
            
            types = repmat({'Pyramidal Cell'}, [1 sum(unitii)]);
            timesAll = uObj.spike_times(unitii);
            
            rawWavAll = uObj.waveforms(unitii);
            uLocation = channelLocations(unitii);
            uBurst = burstIndex(unitii);
            uT2P = troughToPeak(unitii);

            nUnit = length(types);
            
            

            for unii = 1:nUnit
        
                times = timesAll{unii}; 
                times(times < 0 ) = [];
                spikeTimes = times;
                unitsFinal{count,7} = uBurst(unii);
                unitsFinal{count,6} = uT2P(unii);
                unitsFinal{count,5} = uLocation{unii};
                unitsFinal{count,4} = rawWavAll{unii} - rawWavAll{unii}(1);
                unitsFinal{count,2} = spikeTimes;
                unitsFinal{count,1} = ch;

                type = types{unii};
                if uT2P(unii) >= 45
                    unitsFinal{count,3} = 'pyr';
                elseif uT2P(unii) < 45
                    unitsFinal{count,3} = 'int';
                elseif contains(type,'mult') 
                    unitsFinal{count,3} = 'mult';
                end
         
                count = count + 1;
                
            end

        case 'Anna'

            fldr = spikeFolders{nF};
            findCh = find(ismember(fldr,'Ch'),1,'last');
            ch = str2double(fldr(findCh+1:findCh+3));
        
            files = dir(fullfile(spikeDirectory,fldr,'*.mat'));
            files = {files.name};
            
            if length(files) > 1
                error('%s has more than one mat file\n', fullfile(spikeDirectory,fldr))
            end
            
            matFile = load(fullfile(spikeDirectory,fldr,files{1}));
            
        
            if ismember(ch,pyrChannels)
               
                ch_ind = find(ismember(pyrChannels,ch));
        %         spikeTimes = [];
        
                for unii = 1:length(pyr.unit{ch_ind})
        %                 spikeTimes = [];
        
                        un = pyr.unit{ch_ind}{unii};
                        times = matFile.spike.timestamp{un}-spin_offset; % need to decrease this by the spin_offset
                        times(times < 0 ) = [];
                        times = times / 1e3;
        %                 spikeTimes = [spikeTimes times];
                        spikeTimes = times;
        %         end
                
                        unitsFinal{count,2} = spikeTimes;
                        unitsFinal{count,1} = ch;
                        unitsFinal{count,3} = 'pyr';
                
                        count = count + 1;
                end
            end
            
             if ismember(ch,intChannels)
                
                ch_ind = find(ismember(intChannels,ch));
                
        %         spikeTimes = [];
                for unii = 1:length(int.unit{ch_ind})
                        un = int.unit{ch_ind}{unii};
                        times = matFile.spike.timestamp{un}-spin_offset; % need to decrease this by the spin_offset
                        times(times < 0 ) = [];
                        times = times / 1e3;
        %                 spikeTimes = [spikeTimes times];
                        spikeTimes = times;
        %         end
                
                        unitsFinal{count,2} = spikeTimes;
                        unitsFinal{count,1} = ch;
                        unitsFinal{count,3} = 'int';
                
                        count = count + 1;
                end
             end
            
             if ismember(ch,multChannels)
                
                ch_ind = find(ismember(multChannels,ch));
                
        %         spikeTimes = [];
                for unii = 1:length(mult.unit{ch_ind})
                        un = mult.unit{ch_ind}{unii};
                        times = matFile.spike.timestamp{un}-spin_offset; % need to decrease this by the spin_offset
                        times(times < 0 ) = [];
                        times = times / 1e3;
        %                 spikeTimes = [spikeTimes times];
                        spikeTimes = times;
        
                        unitsFinal{count,2} = spikeTimes;
                        unitsFinal{count,1} = ch;
                        unitsFinal{count,3} = 'mult';
                        
                        count = count + 1;
                end
        
        
                
             end
    

    end
end


















return