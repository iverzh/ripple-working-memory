

close all


fs = 1e3;
fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
[b,a] = butter(3,[70 100]/(fs/2));
[bL,aL] = butter(3,250/(fs/2), 'low');

ripWin = 200;




for subj = 37 %:length(subj_list_full)
        subject = subj_list_full{subj};

        fprintf('processing %s ... \n', subject)
        
        if ismember(subj, replay_event{1})
            

            units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
            if isempty(units); continue; end
            U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
            uChan = cell2mat(units(:,1));
            
            f = contains(LFPfilesMicr, subject);
            micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
            chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');

            f = sprintf('%s/%s/%s_1kHz_unitsRemoved.mat', fldr, '../data_1kHz', subject);
            load(f, 'data', 'lfp_data')

            splitPattern = '(_right|_left)';
            locations = regexp(chan_labels, splitPattern, 'split')';
            locations(cellfun(@(X) isempty(X), locations)) = [];
            locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
            locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
            locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
            locations = strrep(locations, 'amygdala', 'AMY');
            locations = strrep(locations, 'hippocampus', 'HIP');


            hem = regexp(chan_labels, splitPattern, 'match')';
            hemi = contains(hem, 'left');
            locations(hemi) = cellfun(@(X) ['L' X], locations(hemi),  'UniformOutput', false);
            hemi = contains(hem, 'right');
            locations(hemi) = cellfun(@(X) ['R' X], locations(hemi),  'UniformOutput', false);
            
            modifier = '1kHz_template_z25';
            tag = [recordingState,'_',location,'_',modifier];
            filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
            microObj = load(fullfile(matExportFolder, filename));
            rippleStats = microObj.rippleStats;
            chanNum = str2double(rippleStats.chanLabels);

            rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
            for chRipp = 1:size(rippMask,1) %sum(strcmp(rippleStats.recordingType, 'macro'))+1:size(rippMask,1) 
                if rippleStats.density{chRipp} > 1
        %   

                    iS = round(rippleStats.window{chRipp}(:,1));
                    iE = round(rippleStats.window{chRipp}(:,2));

                    for ii = 1:length(iE)
        %                 iS(ii) = find(times == iS(ii));
        %                 iE(ii) = find(times == iE(ii));
                        if any([iS(ii) iE(ii)] <= 0); continue; end
                        rippMask(chRipp,iS(ii):iE(ii)) = 1;
                    end
                end

            end
            data(isnan(data)) = 0;
            for ch = 1:size(data,1)
                                
                data(ch,:) = filtfilt(bL,aL,data(ch,:)); % zero phase
            end
            rippMask(:, end:end+3e3) = nan; %pad the ending
            data(:, end:end+3e3) = nan;



            
            
            f = contains(taskfiles, subject);
            trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
            
            iiSubj = find(ismember(replay_event{1}, subj) & ismember(replay_event{4}, 89));
            
            for iS = 1:length(iiSubj)
                iT = replay_event{4}(iiSubj(iS));
                uA = 21;replay_event{2}(iiSubj(iS));
                chA = find(chanNum == uChan(uA));
                uB = 42;replay_event{3}(iiSubj(iS));
                chB = find(chanNum == uChan(uB));
                
                
                if isempty(chA) || isempty(chB); continue; end
                
                threshA = min(rippleStats.rippleAmp{chA});
                threshB = min(rippleStats.rippleAmp{chB});
                
                locA = locations{chA}; locB = locations{chB};
                datArip = sum(rippMask(contains(locations, locA),:), 'omitnan') > 0;
                datBrip = sum(rippMask(contains(locations, locB),:), 'omitnan') > 0;
                datCoR = datArip > 0 & datBrip > 0;
                timesA = ceil(units{uA,2}*fs);
                timesB = ceil(units{uB,2}*fs);
                
                probeTime =  round(trials.timestamps_Probe(iT)*1e3);

                
                
                bandpassA = filtfilt(b,a,data(chA,probeTime+1:probeTime+1e3)); % zero phase
                bandpassB = filtfilt(b,a,data(chB,probeTime+1:probeTime+1e3)); % zero phase
                
                probeA = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3 & datCoR(timesA))-(probeTime+1);
                probeB = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3 & datCoR(timesB))-(probeTime+1);
                
                probeAall = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3)-(probeTime+1);
                probeBall = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3)-(probeTime+1);
                
                s = probeA - probeB';
                [B, A] = find(abs(s) <= 26);
                for ii = 1:length(A)
                    cntr = round(mean([probeA(A(ii)), probeB(B(ii))]));
                    if cntr <= 50 || cntr > length(bandpassA) - 50; continue; end
                    
                    if ~any(abs(hilbert(bandpassA(cntr- 50 : cntr + 50))) > threshA) || ~any(abs(hilbert(bandpassB(cntr- 50 : cntr + 50))) > threshB)
                        continue
                    end
                    
                    
                    
                    figure('Position', [96 1198 1235 335]); 
                    subplot(2,2,1)
                    yyaxis left
                    plot(data(chA, probeTime+1:probeTime+1e3), 'k');
                    ylim([-70 70])

                    yyaxis right
                    pl = plot(bandpassA); hold on;
                    pl.Color = [0 .7 0 0.5];
                    plot(datArip(probeTime+1:probeTime+1e3));
                    vline(probeAall)
                    xlim([cntr- ripWin cntr + ripWin])
                    ylim([-15 15])

                    subplot(2,2,3)
                    yyaxis left
                    plot(data(chB, probeTime+1:probeTime+1e3), 'k');
                    ylim([-70 70])

                    yyaxis right
                    pl = plot(bandpassB);hold on;
                    pl.Color = [0 .7 0 0.5];
                    plot(datBrip(probeTime+1:probeTime+1e3));
                    vline(probeBall)
                    xlim([cntr- ripWin cntr + ripWin])
                    ylim([-15 15])

                    encIM = [trials.PicIDs_Encoding1(iT) trials.PicIDs_Encoding2(iT) trials.PicIDs_Encoding3(iT)];
                    encTall = round([trials.timestamps_Encoding1(iT) trials.timestamps_Encoding2(iT) trials.timestamps_Encoding3(iT)]*1e3);
                    Em = find(encIM == trials.PicIDs_Probe(iT));
                    encodeTime = [encTall(Em), encTall(Em)+2e3];

                    bandpassA = filtfilt(b,a,data(chA,encodeTime(1):encodeTime(2))); % zero phase
                    bandpassB = filtfilt(b,a,data(chB,encodeTime(1):encodeTime(2))); % zero phase
                    encodeA = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2) & datCoR(timesA))-encodeTime(1);
                    encodeB = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2) & datCoR(timesB))-encodeTime(1);
                    
                    encodeAall = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2))-encodeTime(1);
                    encodeBall = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2))-encodeTime(1);
                    s = encodeA - encodeB';
                    [B2, A2] = find(abs(s) <= 26);
                    for jj = 1:length(B2)
                        cntr = round(mean([encodeA(A2(jj)), encodeB(B2(jj))]));
                        if cntr <= 50 || cntr > length(bandpassA) - 50; continue; end
                       
                        if ~any(abs(hilbert(bandpassA(cntr- 50 : cntr + 50))) > threshA) || ~any(abs(hilbert(bandpassB(cntr- 50 : cntr + 50))) > threshB)
                            continue
                        end
                        subplot(2,2,2)
                        yyaxis left
                        pl = plot(data(chA, encodeTime(1):encodeTime(2)), 'k');
                        ylim([-70 70])

                        yyaxis right
                        pl = plot(bandpassA); hold on;
                        pl.Color = [0 .7 0 0.5];
                        plot(datArip(encodeTime( 1):encodeTime(2)));
                        vline(encodeAall)
                        xlim([cntr- ripWin cntr + ripWin])
                        ylim([-15 15])


                        subplot(2,2,4)
                        yyaxis left
                        plot(data(chB, encodeTime(1):encodeTime(2)), 'k');
                        ylim([-70 70])

                        yyaxis right
                        pl = plot(bandpassB);hold on;
                        pl.Color = [0 .7 0 0.5];
                        ylim([-15 15])

                        plot(datBrip(encodeTime(1):encodeTime(2)));
                        vline(encodeBall)
                        xlim([cntr- ripWin cntr + ripWin])
                         sgtitle(sprintf('A: %s B: %s', locA,locB))
                        
                        fig = gcf;
                        fig.Color = 'w';
                        savepdf(gcf, sprintf('%s_uA%i_uB%i_iT%i_%i_%i_rippMask.pdf', subject, uA, uB,iT,ii,jj))

                        % Get the current subplot
                        ax = gca;
                        subplot(2,2,2)
                        % Clear left y-axis
                        yyaxis(ax, 'left');
                        cla;

                        % Clear right y-axis
                        yyaxis(ax, 'right');
                        cla;
                        subplot(2,2,4)
                        % Clear left y-axis
                        yyaxis(ax, 'left');
                        cla;

                        % Clear right y-axis
                        yyaxis(ax, 'right');
                        cla;
                    end
                    
                    close all
                
                end

                
            end
           

        end


end



%%



close all


fs = 1e3;
fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
[b,a] = butter(3,[70 100]/(fs/2));
[bL,aL] = butter(3,250/(fs/2), 'low');

ripWin = 200;




for subj = 37 %:length(subj_list_full)
        subject = subj_list_full{subj};

        fprintf('processing %s ... \n', subject)
        
        if ismember(subj, replay_event{1})
            

            units = LoadSpikeTimes(subject,'RutishauserLab', 'Sternberg');
            if isempty(units); continue; end
            U = cellfun(@(x) any(strcmp(x, {'pyr', 'int', 'mult'})), units(:,3));
            uChan = cell2mat(units(:,1));
            
            f = contains(LFPfilesMicr, subject);
            micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
            chan_labels = regexprep(micrObj.chan_locations, '[^a-zA-Z_]', '');

            f = sprintf('%s/%s/%s_1kHz_unitsRemoved.mat', fldr, '../data_1kHz', subject);
            load(f, 'data', 'lfp_data')

            splitPattern = '(_right|_left)';
            locations = regexp(chan_labels, splitPattern, 'split')';
            locations(cellfun(@(X) isempty(X), locations)) = [];
            locations = strrep(locations, 'ventral_medial_prefrontal_cortex', 'OFC');
            locations = strrep(locations, 'dorsal_anterior_cingulate_cortex', 'ACC');
            locations = strrep(locations, 'pre_supplementary_motor_area', 'SMA');
            locations = strrep(locations, 'amygdala', 'AMY');
            locations = strrep(locations, 'hippocampus', 'HIP');


            hem = regexp(chan_labels, splitPattern, 'match')';
            hemi = contains(hem, 'left');
            locations(hemi) = cellfun(@(X) ['L' X], locations(hemi),  'UniformOutput', false);
            hemi = contains(hem, 'right');
            locations(hemi) = cellfun(@(X) ['R' X], locations(hemi),  'UniformOutput', false);
            
            modifier = '1kHz_template_z25';
            tag = [recordingState,'_',location,'_',modifier];
            filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
            microObj = load(fullfile(matExportFolder, filename));
            rippleStats = microObj.rippleStats;
            chanNum = str2double(rippleStats.chanLabels);

            rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
            for chRipp = 1:size(rippMask,1) %sum(strcmp(rippleStats.recordingType, 'macro'))+1:size(rippMask,1) 
                if rippleStats.density{chRipp} > 1
        %   

                    iS = round(rippleStats.window{chRipp}(:,1));
                    iE = round(rippleStats.window{chRipp}(:,2));

                    for ii = 1:length(iE)
        %                 iS(ii) = find(times == iS(ii));
        %                 iE(ii) = find(times == iE(ii));
                        if any([iS(ii) iE(ii)] <= 0); continue; end
                        rippMask(chRipp,iS(ii):iE(ii)) = 1;
                    end
                end

            end
            data(isnan(data)) = 0;
            for ch = 1:size(data,1)
                                
                data(ch,:) = filtfilt(bL,aL,data(ch,:)); % zero phase
            end
            rippMask(:, end:end+3e3) = nan; %pad the ending
            data(:, end:end+3e3) = nan;



            
            
            f = contains(taskfiles, subject);
            trials = load(fullfile(dataDirectory, '../task', taskfiles{f}));
            
            iiSubj = find(~trials.probe_in_out);
            
            for iS = 1:length(iiSubj)
                iT = iiSubj(iS); %replay_event{4}(iiSubj(iS));
                uA = 21;
                chA = find(chanNum == uChan(uA));
                uB = 42;
                chB = find(chanNum == uChan(uB));
                
                
                if isempty(chA) || isempty(chB); continue; end
                
                threshA = min(rippleStats.rippleAmp{chA});
                threshB = min(rippleStats.rippleAmp{chB});
                
                locA = locations{chA}; locB = locations{chB};
                datArip = sum(rippMask(contains(locations, locA),:), 'omitnan') > 0;
                datBrip = sum(rippMask(contains(locations, locB),:), 'omitnan') > 0;
                datBnorip = sum(rippMask(contains(locations, locB),:), 'omitnan') == 0;
                datCoR = datArip > 0 & datBrip > 0;
                datOneR = (datArip > 0 & datBrip == 0) | (datArip == 0 & datBrip > 0);
                timesA = ceil(units{uA,2}*fs);
                timesB = ceil(units{uB,2}*fs);
                
                probeTime =  round(trials.timestamps_Probe(iT)*1e3);

                
                
                bandpassA = filtfilt(b,a,data(chA,probeTime+1:probeTime+1e3)); % zero phase
                bandpassB = filtfilt(b,a,data(chB,probeTime+1:probeTime+1e3)); % zero phase
                
                probeA = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3 & datOneR(timesA))-(probeTime+1);
                probeB = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3 & datOneR(timesB))-(probeTime+1);
                

                
                probeAall = timesA(timesA >= probeTime+1 & timesA <= probeTime+1e3)-(probeTime+1);
                probeBall = timesB(timesB >= probeTime+1 & timesB <= probeTime+1e3)-(probeTime+1);
                if isempty(probeBall); continue; end
                s = probeA - probeB';
                [B, A] = find(abs(s) <= 26);
                for ii = 1:length(probeA)
%                     cntr = round(mean([probeA(A(ii)), probeB(B(ii))]));
                    cntr = round(probeA(ii));
                    if cntr <= 50 || cntr > length(bandpassA) - 50; continue; end
                    
%                     if ~any(abs(hilbert(bandpassA(cntr- 50 : cntr + 50))) > threshA) || ~any(abs(hilbert(bandpassB(cntr- 50 : cntr + 50))) > threshB)
%                         continue
%                     end

                        if any(probeBall > cntr - 50 & probeBall < cntr + 50)
                            
                            continue
                        end
                    
                    
                    
                    figure('Position', [96 1198 1235 335]); 
                    subplot(2,2,1)
                    yyaxis left
                    plot(data(chA, probeTime+1:probeTime+1e3), 'k');
                    ylim([-70 70])

                    yyaxis right
                    pl = plot(bandpassA); hold on;
                    pl.Color = [0 .7 0 0.5];
                    plot(datArip(probeTime+1:probeTime+1e3));
                    vline(probeAall)
                    xlim([cntr- ripWin cntr + ripWin])
                    ylim([-15 15])

                    subplot(2,2,3)
                    yyaxis left
                    plot(data(chB, probeTime+1:probeTime+1e3), 'k');
                    ylim([-70 70])

                    yyaxis right
                    pl = plot(bandpassB);hold on;
                    pl.Color = [0 .7 0 0.5];
                    plot(datBrip(probeTime+1:probeTime+1e3));
                    vline(probeBall)
                    xlim([cntr- ripWin cntr + ripWin])
                    ylim([-15 15])

                    encIM = [trials.PicIDs_Encoding1(iT) trials.PicIDs_Encoding2(iT) trials.PicIDs_Encoding3(iT)];
                    encTall = round([trials.timestamps_Encoding1(iT) trials.timestamps_Encoding2(iT) trials.timestamps_Encoding3(iT)]*1e3);
                    Em = find(encIM == trials.PicIDs_Probe(iT));
                    encodeTime = [encTall(1), encTall(1)+2e3];

                    bandpassA = filtfilt(b,a,data(chA,encodeTime(1):encodeTime(2))); % zero phase
                    bandpassB = filtfilt(b,a,data(chB,encodeTime(1):encodeTime(2))); % zero phase
                    encodeA = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2) & datOneR(timesA))-encodeTime(1);
                    encodeB = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2) & datOneR(timesB))-encodeTime(1);
                    
                    encodeAall = timesA(timesA >= encodeTime(1) & timesA <= encodeTime(2))-encodeTime(1);
                    encodeBall = timesB(timesB >= encodeTime(1) & timesB <= encodeTime(2))-encodeTime(1);
                    
                    if isempty(encodeBall); continue; end
                    s = encodeA - encodeBall';
                    [B2, A2] = find(abs(s) >= 100 & abs(s) <= 500);
                    for jj = 1:length(encodeA)
                        cntr = encodeA(jj);
                        if cntr <= 50 || cntr > length(bandpassA) - 50; continue; end
                       
%                         if ~any(abs(hilbert(bandpassA(cntr- 50 : cntr + 50))) > threshA) 
%                             continue
%                         end
                        
                        if any(encodeBall > cntr - 50 & encodeBall < cntr + 50)
                            
                            continue
                        end
                        subplot(2,2,2)
                        yyaxis left
                        pl = plot(data(chA, encodeTime(1):encodeTime(2)), 'k');
                        ylim([-70 70])

                        yyaxis right
                        pl = plot(bandpassA); hold on;
                        pl.Color = [0 .7 0 0.5];
                        plot(datArip(encodeTime( 1):encodeTime(2)));
                        vline(encodeAall)
                        xlim([cntr- ripWin cntr + ripWin])
                        ylim([-15 15])


                        subplot(2,2,4)
                        yyaxis left
                        plot(data(chB, encodeTime(1):encodeTime(2)), 'k');
                        ylim([-70 70])

                        yyaxis right
                        pl = plot(bandpassB);hold on;
                        pl.Color = [0 .7 0 0.5];
                        ylim([-15 15])

                        plot(datBrip(encodeTime(1):encodeTime(2)));
                        vline(encodeBall)
                        xlim([cntr- ripWin cntr + ripWin])
                         sgtitle(sprintf('A: %s B: %s', locA,locB))
                        
                        fig = gcf;
                        fig.Color = 'w';
                        savepdf(gcf, sprintf('%s_uA%i_uB%i_iT%i_%i_%i_mismatch.pdf', subject, uA, uB,iT,ii,jj))

                        % Get the current subplot
                        ax = gca;
                        subplot(2,2,2)
                        % Clear left y-axis
                        yyaxis(ax, 'left');
                        cla;

                        % Clear right y-axis
                        yyaxis(ax, 'right');
                        cla;
                        subplot(2,2,4)
                        % Clear left y-axis
                        yyaxis(ax, 'left');
                        cla;

                        % Clear right y-axis
                        yyaxis(ax, 'right');
                        cla;
                    end
                    
                    close all
                
                end

                
            end
            

        end


end





























































































