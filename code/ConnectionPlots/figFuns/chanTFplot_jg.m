
function hand = chanTFplot_jg(cleanedEvents, ch, win, fs, subject, chan_labels, mode, maxfreq, minfreq)

% 10/12/21 pulled from Ilya to add max plot frequency argument
if ~exist('minfreq','var')
    minfreq=1;
end
switch mode
    case 'TF'
        h = figure(13); clf;
        data = cleanedEvents.raw{ch};

        if size(data,1) > 5000
            ind  = randperm(size(data,1));
            data = data(ind(1:5000),:)';
        else
            data = data';
        end

        [ersp, ~, ~, times, freqs, ~] = newtimef(data, 2*win+1, [-win win], fs,...
                    'cycles', [0], ...
                    'alpha', 0.05, ...
                    'plotitc', 'off',...
                    'winsize', fs, ...
                    'maxfreq', maxfreq, 'freqs', [minfreq maxfreq], 'freqscale', 'log',...
                    'plotphase', 'off', ...
                    'padratio', 4, 'naccu', 100,'baseline', [-win -win+0.5*fs], 'baseboot', [-win -win+0.5*fs], 'mcorrect', 'fdr');

       

%         xticks([2 4 10 16 70 190])
        % extract plot data
          h3 = get(h,'Children'); 

%         h3(4).Children(2)
%         xd(1) = h3(4).Children(2).XData(1);   
%         xd(2) = h3(4).Children(2).XData(end);  
        close 
        
       


       

        hand = figure('Position', [316 30 815 564], 'Visible', 'off');
        
        %fix power spikes in the notch filtered frequencies 
        
        for n = 60:60:fs/2
            inds = find(abs(freqs-n) < 1);
            for i = inds
                delta = 1;
                while ismember(i+delta,inds) || ismember(i-delta,inds)
                    delta = delta+1;
                end
                
                ersp(i,:) = mean(ersp([i-delta, i+delta],:));
            end
        end
        
     
        
%         yyaxis left
        %im = imagesc([-fs*1.5 fs*1.5], [1 fs/2], ersp);hold on
        im = imagesc([-fs*1.5 fs*1.5], [minfreq maxfreq], ersp);hold on
        ax = gca;
        ax.YDir = 'normal';
        ax=imgca;
        set(ax,'YScale', 'log')
        % xticks([xd(1) (2/3)*xd(1) (1/3)*xd(1) 0 (1/3)*xd(2) (2/3)*xd(2) xd(2)])
        yTickVals=[1 2 4 10 16 30 40 60 90 120 200 400];
%         yTickVals(yTickVals>0.9*maxfreq)=[];
        yticks(yTickVals)
        caxis([-max(abs([prctile(ersp(:),0.01) prctile(ersp(:),99.99)])) max(abs([prctile(ersp(:),0.01) prctile(ersp(:),99.99)]))])
        colormap(jet)
        %xlim([-fs*1.5 fs*1.5])
        ylim([minfreq maxfreq])
        v  = vline(0, 'k:');
        v.LineWidth = 1.5;
        xlabel('time to ripple center (ms)')
        ylabel('frequency (Hz)')
        c=colorbar;
        c.Label.String = 'ERSP (dB)';
        set(gca, 'FontSize', 10)
        set(gcf, 'Color', [1 1 1])
        
        ax.XTickLabel = ax.XTick * (1/fs) * 1000; %xaxis in ms 
        xlabel('time to ripple center [ms]')
        
        
        rippleTimeWindow = find(times >= -50 & times <= 50);
        rippleFreqWindow = find(freqs >= 60 & freqs <= 120);
        rippleFreqs = freqs(rippleFreqWindow);
        rippleBand = ersp(rippleFreqWindow, rippleTimeWindow);
        meanRipple = mean(rippleBand,2);
        rippleOscFreq = rippleFreqs(meanRipple == max(meanRipple));

        nEvents = size(cleanedEvents.raw{ch},1);
        density = nEvents/sum(cleanedEvents.recordingLength);
        density = density * fs * 60;
        title(sprintf('%s, %s, center freq = %.1f Hz (%i events)', ...
                      subject, chan_labels{ch}, rippleOscFreq, nEvents))
                  
%         yyaxis right
%         meanLFP(ch,:) = mean(cleanedEvents.raw{ch});
%         center = round(length(meanLFP(ch,:))/2);
%         ln = plot(-fs*1.5:fs*1.5, meanLFP(ch,center-fs*1.5:center+fs*1.5), 'k-'); hold on
%         ln.LineWidth = 1.5;
        
        ax = gca;
        ax.YAxis.Color = 'k';
        
    case 'LFP'
        meanLFP(ch,:) = mean(cleanedEvents.raw{ch},1);
    
%         yyaxis right
        hand = figure('Position', [316 30 815 564], 'Visible', 'off');

        center = round(length(meanLFP(ch,:))/2);
        ln = plot(-win:win, meanLFP(ch,center-win:center+win), 'k-'); hold on
        ln.LineWidth = 3;
% 
%         ax = gca;
%         ax.YAxis(1).Color = 'b';
%         ax.YAxis(2).Color = 'k';

end
fig = gcf;
fig.Color = [1 1 1];


return


