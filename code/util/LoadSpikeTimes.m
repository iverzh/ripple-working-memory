
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

    case 'RutishauserLab'
        spikeDirectory = sprintf('/space/seh10/6/halgdev/projects/iverzh/data/%s/preprocess/units', study);
        files = dir(sprintf('%s/*%s*_units.mat', spikeDirectory, subject));
        files = {files.name};
        filename = files{1};
        uObj = load(fullfile(spikeDirectory, filename));
       
        detectedChannels = unique(uObj.electrode_id);
        pyrChannels = detectedChannels;
        

        if strcmp(study,'Sternberg')
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

    end
end


















return