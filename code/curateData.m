clear 
close all

addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/RippleDetection'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))

addpath(genpath('/space/seh10/6/halgdev/projects/cdickey/packages'))
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/PNAS/code/LFP_Tools/'))

% addpath(genpath('/home/jgarret/mat')) %mask2bounds
addpath(genpath('/home/jgarret/ArtifactDetection'))
addpath('/home/jgarret/CortRipple/')
% addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/'))
% rmpath(genpath('/space/seh8/1/halgdev/projects/cdickey/eeglab14_1_2b/functions/octavefunc'))

%%

fldr = '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/out';
% fldr = '/space/seh10/6/halgdev/projects/iverzh/data/Sternberg/preprocess';
flst = dir(fullfile(fldr, '*LFP_micro*'));
flst = {flst.name}';
flst = strrep(flst, '.mat', '');
subj_list_full = flst;

runList = 29:length(subj_list_full); %[1,2,4,14,15];%[1,2,4,5,6,7,8,9,11,15]; %[1,2,4,5,6,7,8,9,11,15];%[6,8,10,13,14];%[1,3,4,19,20,21];%[1,2,5,6,7,8,9,11,14];%[1,8,11];
recordingState = 'wake'; % for cc wake, use sleep established polarity check
locations = {'NC'}; %{'NC', 'TH'};%,'HC','TH'}; %{'NC','HC','TH','AMY'};
% modifier = 'ECoG_WN_taskPeriod_v3';
modifier = 'IISremv_singleWire_v3';

masterFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/Figures';
if ~isfolder(masterFolder); mkdir(masterFolder); end

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
% matExportFolder = '/space/seh8/1/halgdev/projects/cdickey/shared/matFiles';
if ~isfolder(matExportFolder); mkdir(matExportFolder); end

preProcExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/preProc';
% preProcExportFolder = '/space/seh8/1/halgdev/projects/cdickey/shared/preProc';
if ~isfolder(preProcExportFolder); mkdir(preProcExportFolder); end
parcelDir = '/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/CC/PPT_labels/';
badChanDir = '/space/seh8/1/halgdev/projects/cdickey/ripple/CortRipple/bad_chan';

scale = 1e6; %1e6; 
IISflag = true; %1 - check for spikes. 0 - do not check for spikes
detrendFlag = false;
artifactFlag = true;

winSec = 2; % pre and post event window in seconds
spikeCheck = 500; %+/- for spike check (10 Hz) abs(LFP) > 1000
rippleWindowSec = 0.100; %+/- window used to compute ripple ocsillation freq. in seconds
IIScheckSec = 0.010; %+/- in sec

fs = 1000; %400; %Hz sample rate
% fs = 2500; %Hz sample rate
% fs = 500;
% rippleband = [120 200];
% rippleband = [65 120];
% rippleband = [55 65];
rippleband = [70 110];

RejectParams.RBHFzscore = [0, 0]; %Reject zscore(Ripple) - zscore(HF) < RBHFzscore
RejectParams.sharpzscore = [7, 1]; %Reject zscore(100+hZ) > sharpzscore
RejectParams.LFPdiff = [50, 0]; %Reject zscore(LFPdiff) > LFPdiffzscore 
RejectParams.LFPdiffzscore = [4, 0]; %Reject stepwise jumps in UAB data
% RejectParams.RBAmp = 2; %Reject RBAmp < 2
RejectParams.RBzscore = [3, 1]; %Reject RBAmp < 2
RejectParams.LFPthresh = [1e3, 1];
RejectParams.minDuration = [0.025, 1]; % in seconds

RejectParams.bndParams.smthKrnl=100;
RejectParams.bndParams.srchWinSz=1/2; %of kernel size
RejectParams.bndParams.scoreThresh=0.75;


% 
sleepfiles_set_CC = {[2,3,7],[2,5,9,10],[3,5,6,7],[1],...
        [2,3,5],[4,10,12],[5,6,7,9,11],[2,3,6,8],[1,2,4,5,7],...
        [3,7,8],[3,7,8,11],[1,2,5],[2,4,5,6,7,8,10,11],[3,8,9,12],...
        [2,6,8],[2,6,7],[1,2,3,5,7,10],[1,3,4,5],[1,3,4]};

sleepfiles_set = {1,1,1,1,1,1,1,1,1,1,1,1,...
                        1,1,1,1,1,1,1,1,1,1,1,1,...
                        1,1,1,1,1,1,1,1,1,1,1,1,...
                        1,1,1,1,1,1,1,1,1,1,1,1,...
                        1,1,1,1,1,1,1,1,1,1,1,1,...
                        1,1,1,1,1,1,1,1,1,1,1,1};

% data = readtable('Ueli Movie Datasat Curation - Bad Channels.csv');
% 
% subj = 2;
% bundles = table2cell(data(2,3:end-1));   
% ii = find(strcmp(bundles, 'R-SMA'));
% badChan = table2cell(data(subj+2, ii+2));
% badChan = str2double(split(badChan,','));
    
    %%
fprintf('Running Ripple Selection ...\n')

for subj = runList% [2,4,7] %[14,15,20] %[1,2,3,5,6,7,8,9,11] %[11,16,24] %[5,8,9,11,16,24] %[19,21,22] %25:27 %  11:numel(subj_list_full) %[19,21,22] %[3,14:17] %[17:20] %[4,5,8,15] %1:numel(subj_list_full) %[4,5,7] %20:22 %[1,3,6,7,10,11,13,16,17] %%  %[14, 15] %[11,13,14,15] %21 %15:17 %[3,4,5,8,10,11,13,14,15]  %loop through subjects
    subject = subj_list_full{subj};
    
    
        

    for l = 1:length(locations)

        location = locations{l};
        
        fprintf('Loading %s Data...\n', location)
        
        tag = [recordingState,'_',location,'_',modifier];
        
        for s = 1:length(sleepfiles_set{subj})
            sleep_ID = sleepfiles_set{subj}(s);
            
            fprintf('Loading segment %i ...\n', sleep_ID)

            f = sprintf('%s/%s.mat', fldr, subject);
            load(f, 'data', 'lfp_data', 'chan_labels')

            if exist('lfp_data','var')
                BroadbandData = lfp_data;
            elseif exist('data','var')
                BroadbandData = data;


            end
                
               
            if ~isa(BroadbandData, 'double')
                BroadbandData = double(BroadbandData);
            end
    
            if size(BroadbandData,1) > size(BroadbandData,2); BroadbandData = BroadbandData'; end
                
            BroadbandData = scale * BroadbandData;
            bundleSize = 8;
            if mod(size(BroadbandData,1),bundleSize) ~= 0; error('%s microwire lfp data missing', subject); end
            
            bbDataSpike = BroadbandData;
            uSubj = strrep(subject,'_LFP_micro', '');
            units = LoadSpikeTimes(uSubj,'RutishauserLab', 'bmovie');
            uCh = cell2mat(units(:,1));
            unitFree = [];
            for iB = 25:bundleSize:size(BroadbandData,1)
                bundleChanAll = iB:(iB+(bundleSize-1));
                unitsBundle = units(ismember(uCh, bundleChanAll),:);
                unitsBundle(:,1) = num2cell(cellfun(@(X) X - (iB - 1), unitsBundle(:,1)));

                chan_labels(bundleChanAll)
                
                figure(1);
                fig = gcf;
                fig.Position = [1 -213 1504 1534];

                if isempty(unitsBundle)
                    lamplot(-BroadbandData(bundleChanAll,:), fs, true, false, 500)
                else
                    lamplot_units(-BroadbandData(bundleChanAll,:), fs, true, true, 500, unitsBundle, chan_labels(bundleChanAll))
                end
                title(chan_labels{iB})


                if IISflag
                    
                    microFlag = false;
                    bbDataSpike = BroadbandData(bundleChanAll,:);
    
                    [spikeMask,spikeBounds,IISscore] = detectSpikes_IAV(bbDataSpike', fs, microFlag);

                    figure(2);
                    fig = gcf;
                    fig.Position =[1505 -213 1504 1534];
    
                    lamplot(-IISscore, fs, true, true, 110)

                    title(subject, 'Interpreter','latex')

                        
                end

                disp('Press "n" to continue to next bundle');

                while true
                    % Wait for a key press
                    key = get(gcf, 'CurrentCharacter');
                    fprintf(key)
                    
                    % If 'q' is pressed, exit the loop
                    if key == 'n'
                        disp('key "n" was pressed. loading next bundle...');
                        break;
                    elseif key == 'q'
                        close all
                        error('key "q" was pressed. quiting curation...');
                    end
                    
                    % Wait for the next key press
                    waitforbuttonpress;
                end
 
                
                
                
                close all;
            end
        end

    end
end