%% Single-Trial Co-Ripple Visualization
% This script generates single-trial visualizations of co-ripple events
% during encoding and probe phases of the Sternberg working memory task.
% It identifies synchronous ripple events across brain regions and compares
% neural activity during memory retrieval.

close all

%% Analysis parameters
fs = 1e3;  % Sampling frequency (Hz)
fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
ripWin = 200;  % Window size for visualization (ms)

% Bandpass filter for ripple band (70-100 Hz)
[b, a] = butter(3, [70 100]/(fs/2));

% Low-pass filter for LFP (250 Hz)
[bL, aL] = butter(3, 250/(fs/2), 'low');

%%  Co-ripple events during successful memory retrieval

for subj = 37  % Process specific subject
    subject = subj_list_full{subj};
    fprintf('Processing %s ... \n', subject)
    
    % Only process if subject has replay events
    if ~ismember(subj, replay_event{1})
        continue
    end
    
    %% Load spike data
    units = LoadSpikeTimes(subject, 'RutishauserLab', 'Sternberg');
    if isempty(units)
        continue
    end
    
    % Filter for valid unit types (pyramidal, interneuron, multi-unit)
    U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:, 3));
    uChan = cell2mat(units(:, 1));
    
    %% Load LFP data and channel labels
    f = contains(LFPfilesMicr, subject);
    micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
    chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');
    
    f = sprintf('%s/%s/%s_1kHz_unitsRemoved.mat', fldr, '../data_1kHz', subject);
    load(f, 'data', 'lfp_data')
    
    %% Parse channel locations and hemispheres
    splitPattern = '(_right|_left)';
    locations = regexp(chan_labels, splitPattern, 'split')';
    locations(cellfun(@(X) isempty(X), locations)) = [];
    
    % Convert to standardized region names
    locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
    locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
    locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
    locations = strrep(locations, 'amygdala', 'AMY');
    locations = strrep(locations, 'hippocampus', 'HIP');
    
    % Add hemisphere labels
    hem = regexp(chan_labels, splitPattern, 'match')';
    hemi = contains(hem, 'left');
    locations(hemi) = cellfun(@(X) ['L' X], locations(hemi), 'UniformOutput', false);
    hemi = contains(hem, 'right');
    locations(hemi) = cellfun(@(X) ['R' X], locations(hemi), 'UniformOutput', false);
    
    %% Load ripple statistics
    modifier = '1kHz_template_z25';
    tag = [recordingState, '_', location, '_', modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    rippleStats = microObj.rippleStats;
    chanNum = str2double(rippleStats.chanLabels);
    
    %% Create ripple mask
    rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
    
    for chRipp = 1:size(rippMask, 1)
        if rippleStats.density{chRipp} <= 1
            continue
        end
        
        iS = round(rippleStats.window{chRipp}(:, 1));
        iE = round(rippleStats.window{chRipp}(:, 2));
        
        for ii = 1:length(iE)
            if any([iS(ii) iE(ii)] <= 0)
                continue
            end
            rippMask(chRipp, iS(ii):iE(ii)) = 1;
        end
    end
    
    %% Preprocess LFP data
    data(isnan(data)) = 0;
    for ch = 1:size(data, 1)
        data(ch, :) = filtfilt(bL, aL, data(ch, :));  % Zero-phase low-pass filter
    end
    
    % Pad ending with NaN
    rippMask(:, end:end+3e3) = nan;
    data(:, end:end+3e3) = nan;
    
    %% Load task data
    f = contains(taskfiles, subject);
    trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
    
    %% Process replay events (co-ripple during successful retrieval)
    iiSubj = find(ismember(replay_event{1}, subj) & ismember(replay_event{4}, 89));
    
    for iS = 1:length(iiSubj)
        iT = replay_event{4}(iiSubj(iS));  % Trial number
        uA = 21;  % Unit A index
        chA = find(chanNum == uChan(uA));
        uB = 42;  % Unit B index
        chB = find(chanNum == uChan(uB));
        
        if isempty(chA) || isempty(chB)
            continue
        end
        
        % Get ripple amplitude thresholds for each channel
        threshA = min(rippleStats.rippleAmp{chA});
        threshB = min(rippleStats.rippleAmp{chB});
        
        % Get location labels
        locA = locations{chA};
        locB = locations{chB};
        
        % Create region-specific ripple masks
        datArip = sum(rippMask(contains(locations, locA), :), 'omitnan') > 0;
        datBrip = sum(rippMask(contains(locations, locB), :), 'omitnan') > 0;
        datCoR = datArip > 0 & datBrip > 0;  % Co-ripple events
        
        % Get spike times for each unit
        timesA = ceil(units{uA, 2} * fs);
        timesB = ceil(units{uB, 2} * fs);
        
        probeTime = round(trials.timestamps_Probe(iT) * 1e3);
        
        %% Analyze probe period
        % Bandpass filter data in ripple band
        bandpassA = filtfilt(b, a, data(chA, probeTime+1:probeTime+1e3));
        bandpassB = filtfilt(b, a, data(chB, probeTime+1:probeTime+1e3));
        
        % Extract spikes during co-ripple events
        probeA = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3 & datCoR(timesA)) - (probeTime+1);
        probeB = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3 & datCoR(timesB)) - (probeTime+1);
        
        % Extract all spikes in probe period
        probeAall = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3) - (probeTime+1);
        probeBall = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3) - (probeTime+1);
        
        % Find synchronous spikes (within 26 ms)
        s = probeA - probeB';
        [B, A] = find(abs(s) <= 26);
        
        for ii = 1:length(A)
            cntr = round(mean([probeA(A(ii)), probeB(B(ii))]));
            
            % Skip if too close to edges
            if cntr <= 50 || cntr > length(bandpassA) - 50
                continue
            end
            
            % Verify ripple presence in both channels
            if ~any(abs(hilbert(bandpassA(cntr-50:cntr+50))) > threshA) || ...
               ~any(abs(hilbert(bandpassB(cntr-50:cntr+50))) > threshB)
                continue
            end
            
            %% Create figure showing probe period
            figure('Position', [96 1198 1235 335]);
            
            % Channel A - Probe
            subplot(2, 2, 1)
            yyaxis left
            plot(data(chA, probeTime+1:probeTime+1e3), 'k');
            ylim([-70 70])
            
            yyaxis right
            pl = plot(bandpassA);
            hold on;
            pl.Color = [0 .7 0 0.5];
            plot(datArip(probeTime+1:probeTime+1e3));
            vline(probeAall)
            xlim([cntr-ripWin cntr+ripWin])
            ylim([-15 15])
            
            % Channel B - Probe
            subplot(2, 2, 3)
            yyaxis left
            plot(data(chB, probeTime+1:probeTime+1e3), 'k');
            ylim([-70 70])
            
            yyaxis right
            pl = plot(bandpassB);
            hold on;
            pl.Color = [0 .7 0 0.5];
            plot(datBrip(probeTime+1:probeTime+1e3));
            vline(probeBall)
            xlim([cntr-ripWin cntr+ripWin])
            ylim([-15 15])
            
            %% Analyze encoding period for matching item
            encIM = [trials.PicIDs_Encoding1(iT) trials.PicIDs_Encoding2(iT) trials.PicIDs_Encoding3(iT)];
            encTall = round([trials.timestamps_Encoding1(iT) trials.timestamps_Encoding2(iT) trials.timestamps_Encoding3(iT)] * 1e3);
            Em = find(encIM == trials.PicIDs_Probe(iT));
            encodeTime = [encTall(Em), encTall(Em)+2e3];
            
            % Bandpass filter encoding period
            bandpassA = filtfilt(b, a, data(chA, encodeTime(1):encodeTime(2)));
            bandpassB = filtfilt(b, a, data(chB, encodeTime(1):encodeTime(2)));
            
            % Extract spikes during co-ripple events in encoding
            encodeA = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2) & datCoR(timesA)) - encodeTime(1);
            encodeB = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2) & datCoR(timesB)) - encodeTime(1);
            
            % Extract all spikes in encoding period
            encodeAall = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2)) - encodeTime(1);
            encodeBall = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2)) - encodeTime(1);
            
            % Find synchronous spikes in encoding
            s = encodeA - encodeB';
            [B2, A2] = find(abs(s) <= 26);
            
            for jj = 1:length(B2)
                cntr = round(mean([encodeA(A2(jj)), encodeB(B2(jj))]));
                
                % Skip if too close to edges
                if cntr <= 50 || cntr > length(bandpassA) - 50
                    continue
                end
                
                % Verify ripple presence in both channels
                if ~any(abs(hilbert(bandpassA(cntr-50:cntr+50))) > threshA) || ...
                   ~any(abs(hilbert(bandpassB(cntr-50:cntr+50))) > threshB)
                    continue
                end
                
                % Channel A - Encoding
                subplot(2, 2, 2)
                yyaxis left
                plot(data(chA, encodeTime(1):encodeTime(2)), 'k');
                ylim([-70 70])
                
                yyaxis right
                pl = plot(bandpassA);
                hold on;
                pl.Color = [0 .7 0 0.5];
                plot(datArip(encodeTime(1):encodeTime(2)));
                vline(encodeAall)
                xlim([cntr-ripWin cntr+ripWin])
                ylim([-15 15])
                
                % Channel B - Encoding
                subplot(2, 2, 4)
                yyaxis left
                plot(data(chB, encodeTime(1):encodeTime(2)), 'k');
                ylim([-70 70])
                
                yyaxis right
                pl = plot(bandpassB);
                hold on;
                pl.Color = [0 .7 0 0.5];
                ylim([-15 15])
                plot(datBrip(encodeTime(1):encodeTime(2)));
                vline(encodeBall)
                xlim([cntr-ripWin cntr+ripWin])
                
                sgtitle(sprintf('A: %s B: %s', locA, locB))
                
                fig = gcf;
                fig.Color = 'w';
                savepdf(gcf, sprintf('%s_uA%i_uB%i_iT%i_%i_%i_rippMask.pdf', subject, uA, uB, iT, ii, jj))
                
                % Clear encoding subplots for next iteration
                subplot(2, 2, 2)
                cla
                subplot(2, 2, 4)
                cla
            end
            
            close all
        end
    end
end







































































