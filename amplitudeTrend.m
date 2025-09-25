

% Edited subjects
% sub-4_ses-1
%         49-52
% sub-29_ses-1
%         4








dataDirectory = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/OrigUpload';
flst = dir(fullfile(dataDirectory, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '_LFP_micro.mat', '');
subj_list_full = flst;


matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults';
if ~isfolder(exportDir); mkdir(exportDir); end
exportDirFigs = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/Tasks';
if ~isfolder(exportDirFigs); mkdir(exportDirFigs); end

unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));

LFPfilesMicr = {LFPfilesMicr.name}';
taskfiles = dir(fullfile(dataDirectory, '../task/*task*'));
taskfiles = {taskfiles.name}';
bpFiles = dir(fullfile(dataDirectory, '../../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../../out', '*micro*'));
micrFiles = {micrFiles.name}';

recordingState = 'wake';
location = 'NC';
% tag = 'vHG';
sfreq = 1e3;
copmuteTF = false; 
% regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', 'LSPE', ...
%            'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP', 'RSPE'};
regions = {'OFC', 'ACC', 'SMA', 'AMY', 'HIP'} ;  
%%

amplDev = [];
subjName = {};
channelName = [];
tot = 0;
thresh = 0.7;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj}

    modifier = '1kHz_template_z25';

    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    
    fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz';
    load(fullfile(fldr, sprintf('%s_1kHz_unitsRemoved.mat', subject)))
    B = [];
    for ch = 1:length(microObj.rippleStats.locs)
        A = false(1, microObj.rippleStats.recordingLength); 
        A(microObj.rippleStats.locs{ch}) = true; % Example array
        
        num_bins = 5;
        bin_size = ceil(microObj.rippleStats.recordingLength / num_bins);

        % Compute padding size
        pad_size = num_bins * bin_size - length(A);

        % Pad with zeros
        A_padded = [A, zeros(1, pad_size)];
        dataCh =  data(ch,:);
        dataCh = [dataCh, zeros(1, pad_size)];
        dataCh = reshape(dataCh, bin_size, []);

        % Reshape and sum
        Btmp = sum(reshape(A_padded, bin_size, []), 1);
        
        Btmp = [Btmp - mean(Btmp)] / mean(Btmp);


        % Reshape and compute the mean, ignoring NaNs
        B(ch,:) = Btmp;
        subjName = [subjName, subject];
        channelName = [channelName, ch];
        
        if any(Btmp > thresh)
            fs = 1000; % Sampling frequency in Hz
            figure('Position', [1219 803 589 770]);
            for iB = 1:num_bins
                x = dataCh(:,iB)'; 

                % Compute the Power Spectral Density (PSD) using Welch's method
                [pxx, f] = pwelch(x, [], [], [], fs);
                pxx = smoothdata(pxx, 'gaussian', 500);

                % Plot the PSD

                pl = plot(f, 10*log10(pxx)); hold on; % Convert power to dB/Hz
                pl.LineWidth = 1.5;
                pl.Color = [0 0 iB/5];
%                 plot(f, pxx); hold on; % Convert power to dB/Hz
                
            end
            xlabel('Frequency (Hz)');
            ylabel('Power/Frequency (dB/Hz)');
            title('Power Spectral Density (PSD)');
            grid on;
            savepdf(gcf, sprintf('%s_chan%i-PSD.pdf', subject, ch))
            
            figure('Position', [271 1160 1518 413]);
            plot(data(ch,:)); 
            vline(microObj.rippleStats.locs{ch});
            savepdf(gcf, sprintf('%s_chan%i-LFP.pdf', subject, ch))
            
            close all

            
        end

        

    end

    amplDev = [amplDev; B];
    tot = tot + length(dev);
end

maxDev = abs(max(amplDev,[],2));

%%

fs = 1000; % Sampling frequency in Hz
t = 0:1/fs:1-1/fs; % Time vector (1 second duration)
x = data(35,:); % Example signal (50Hz + 120Hz + noise)

% Compute the Power Spectral Density (PSD) using Welch's method
[pxx, f] = pwelch(x, [], [], [], fs);

% Plot the PSD
figure;
plot(f, 10*log10(pxx)); % Convert power to dB/Hz
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density (PSD)');
grid on;


%%

amplDev = [];
subjName = {};
channelName = [];
tot = 0;
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj}

    modifier = '1kHz_template_detrend_z25';

    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    microObj = load(fullfile(matExportFolder, filename));
    B = [];
    for ch = 1:size(microObj.detrendDat,1)
        A = microObj.detrendDat(ch,:); % Example array
        bin_size =  round(length(A)/5);
        % Compute how many NaNs need to be added
        pad_size = bin_size - mod(length(A), bin_size);
        if pad_size == bin_size
            pad_size = 0;
        end

        % Pad with NaNs
        A_padded = [A, NaN(1, pad_size)];


        % Reshape and compute the mean, ignoring NaNs
        B(ch,:) = mean(reshape(A_padded, bin_size, []), 1, 'omitnan') - 1;
        subjName = [subjName, subject];
        channelName = [channelName, ch];

    end

    amplDev = [amplDev; B];
    tot = tot + length(dev);
end
maxDev = abs(max(B,[],2));

%%
thresh = 0.5;
sum(amplDev >= thresh) / tot

close all
figure; 
histogram(amplDev, 0:0.1:5)
vline(thresh)


%%
badChan = find(amplDev >= 1);
for iB = 1:length(badChan)
    ii = badChan(iB);
    subject = subjName{ii};
    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    rip = load(fullfile(matExportFolder, filename));
    
    lfp = load(fullfile('/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess/data_1kHz', sprintf('%s_1kHz_unitsRemoved.mat', subject)));
    figure; 
    yyaxis left
    plot(data(channelName(ii),:));
    yyaxis right
    plot(detrendDat(channelName(ii),:));
    
    title(sprintf('%s %i', subject, channelName(ii)))
    waitforbuttonpress; close;

    
end















