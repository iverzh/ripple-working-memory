
close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))
%%

subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P48CS_R1', 'P48CS_R2','P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
dataDirectory =     '/space/seh10/6/halgdev/projects/iverzh/data/bmovie/preprocess';
exportDir = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coRipple';
if ~isfolder(exportDir); mkdir(exportDir); end
exportDirPDF = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/corip';


unitfiles = dir(fullfile(dataDirectory, '*unit*'));
unitfiles = {unitfiles.name}';
LFPfilesMacr = dir(fullfile(dataDirectory, '*macro*'));
LFPfilesMacr = {LFPfilesMacr.name}';
LFPfilesMicr = dir(fullfile(dataDirectory, '*micro*'));
LFPfilesMicr = {LFPfilesMicr.name}';
bpFiles = dir(fullfile(dataDirectory, '../out', '*macro*'));
bpFiles = {bpFiles.name}';
micrFiles = dir(fullfile(dataDirectory, '../out', '*micro*'));
micrFiles = {micrFiles.name}';

recordingState = 'wake';
location = 'NC';

regions = {'LOFC', 'LACC', 'LSMA', 'LAMY', 'LHIP', ...
           'ROFC', 'RACC', 'RSMA', 'RAMY', 'RHIP'};
        
computeCoRipple = false;
rippleFR = nan(length(subj_list_full), 100);
CoRippleRatesRegionAll = nan(length(regions), length(regions), length(subj_list_full));
figure('Position', [912 128 1462 868], 'Units', 'pixels');
figure('Position', [58.4600 138.0700 342.5500 274.6550], 'Units', 'pixels');
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};

     switch computeCoRipple
        case true
            modifier = 'IISremv';
            tag = [recordingState,'_',location,'_',modifier];
            filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'macro_bp', tag);
            load(fullfile(matExportFolder, filename))
            
            modifier = '';
            tag = [recordingState,'_',location,'_',modifier];
            filename = sprintf('%s_LFP_%s_ripple_stats_%s.mat', subject, 'micro', tag);
            microObj = load(fullfile(matExportFolder, filename));
                
            
            f = contains(unitfiles, subject);
            uObj = load(fullfile(dataDirectory, unitfiles{f}));
            
            f = contains(LFPfilesMacr, subject);
            macrObj = load(fullfile(dataDirectory, LFPfilesMacr{f}));
            LFPtimeMacr = macrObj.times;
        
            f = contains(LFPfilesMicr, subject);
            micrObj = load(fullfile(dataDirectory, LFPfilesMicr{f}));
            LFPtimeMicr = micrObj.times;
        
            f = contains(bpFiles, subject);
            bObj = load(fullfile(dataDirectory, '../out', bpFiles{f}));
        
            f = contains(micrFiles, subject);
            micrObj = load(fullfile(dataDirectory, '../out', micrFiles{f}));
            
            rippleStats.recordingType = [repmat({'macro'}, [1 length(rippleStats.locs)]) repmat({'micro'}, [1 length(microObj.rippleStats.locs)])];
            rippleStats.locs = [cellfun(@(X) LFPtimeMacr(X), rippleStats.locs, 'UniformOutput', false), ...
                                cellfun(@(X) LFPtimeMicr(X), microObj.rippleStats.locs, 'UniformOutput', false)];
            rippleStats.window(cellfun(@(X) isempty(X), rippleStats.window)) = {[1 1]};
            microObj.rippleStats.window(cellfun(@(X) isempty(X),  microObj.rippleStats.window)) = {[1 1]};
            rippleStats.window = [cellfun(@(X) [LFPtimeMacr(X(:,1)); LFPtimeMacr(X(:,2))]', rippleStats.window, 'UniformOutput', false), ...
                                  cellfun(@(X) [LFPtimeMicr(X(:,1)); LFPtimeMicr(X(:,2))]', microObj.rippleStats.window, 'UniformOutput', false)];
        
            rippleStats.chanLabels = [rippleStats.chanLabels; microObj.rippleStats.chanLabels];
            rippleStats.macroTimes = LFPtimeMacr;
            rippleStats.microTimes = LFPtimeMicr;
            
        
            chA_all = find(strcmp(rippleStats.recordingType, 'micro')); %1:length(rippleStats.locs);
            chB_all = find(strcmp(rippleStats.recordingType, 'micro')); %1:length(rippleStats.locs);
        
            rippMask = zeros(length(rippleStats.chanLabels), rippleStats.recordingLength);
            for chRipp = 1:size(rippMask,1) 
                if length(rippleStats.window{chRipp}) > 2
%                     if strcmp(rippleStats.recordingType{chRipp}, 'macro'); times = rippleStats.macroTimes;
                    if strcmp(rippleStats.recordingType{chRipp}, 'macro'); continue;
                    elseif strcmp(rippleStats.recordingType{chRipp}, 'micro'); times = rippleStats.microTimes; end
                    
                    iS = round(rippleStats.window{chRipp}(:,1) * 1e3);
                    iE = round(rippleStats.window{chRipp}(:,2) * 1e3);
                
                    for ii = 1:length(iE)
%                         iS(ii) = find(times == iS(ii));
%                         iE(ii) = find(times == iE(ii));
                        if any([iS(ii) iE(ii)] <= 0); continue; end
                        rippMask(chRipp,iS(ii):iE(ii)) = 1;
                    end
                end
                
            end
        
            totalRipples = sum(rippMask);
            
            unitTimes = cell2mat(uObj.spike_times);
            nUnits = length(uObj.spike_times);
            baselineFR = length(unitTimes) / nUnits / (rippleStats.recordingLength / rippleStats.fs);
            for nR = 0:max(totalRipples)
                b = mask2bounds(totalRipples == nR) / 1e3;
                nSpike = 0;
                if size(b,1) > 1
                    for iB = 1:length(b)
                        nSpike = nSpike + sum(unitTimes >= b(iB,1) & unitTimes <= b(iB,2));
                    end
                elseif ~isempty(b)
                    nSpike = nSpike + sum(unitTimes >= b(1) & unitTimes <= b(2));
                end
                dur = sum(b(:,2) - b(:,1));
        
                if dur >  1
                    FR = nSpike / dur / nUnits;
                    rippleFR(subj, nR+1) = FR / baselineFR;
                end
        
            end
            
        %     clf
        %     boundedline(0:99, mean(rippleFR, 'omitnan'), std(rippleFR, 'omitnan')/sqrt(length(subj_list_full)));
        %     xlim([0 15])
        
        
        
            coords = [bObj.bpCoords; micrObj.chan_coords];    
   
            CoRippleRates = nan(length(chA_all), length(chB_all));
            CoRippleN = nan(length(chA_all), length(chB_all));
            euDistances = nan(length(chA_all), length(chB_all));
            for iiA = 1:length(chA_all)
        
                for iiB = 1:length(chB_all)
                    chA = chA_all(iiA);
                    chB = chB_all(iiB);
                    if isempty(rippleStats.locs{chA}) || isempty(rippleStats.locs{chB}) || ...
                       (chA == chB)

                        continue
                    end
                    s = (rippleStats.locs{chA} - rippleStats.locs{chB}');
                    r = sum(abs(s(:)) < (0.05)) / length(rippleStats.locs{chA});
                    
%                     bnd = mask2bounds(rippMask(chA,:) & rippMask(chB,:));
%                     bnd = bnd(:,2) - bnd(:,1);
%                     r = sum(bnd > 25) / length(rippleStats.locs{chA});
                    CoRippleRates(iiA,iiB) = r; 
                    CoRippleN(iiA,iiB) = sum(abs(s(:)) < (0.05)); 
            
                    labelA = rippleStats.chanLabels{chA};
                    labelB = rippleStats.chanLabels{chB};
    
                    euDistances(iiA,iiB) = pdist([coords(chA,:); coords(chB,:)]);
    
            
                
                %             bnd = mask2bounds(rippMask(chA,:) & rippMask(chB,:));
                %             bnd = bnd(:,2) - bnd(:,1);
                %             coRippRate = sum(bnd > (25*(rippleStats.fs/1e3))) / length(rippleStats.locs{chA});
        %             CoRippleRates(chA,chB) = coRippRate;
                    
                    if mod(chB, 5) == 1
                        clc
                        fprintf('computing coRipple Rates for session %s\n', subject)
                        fprintf('channel %i / %i\n', chA_all(iiA), length(chA_all))
                        
                        prog = ceil(iiB/length(chB_all)*10);
                        Bar = repmat('*', [1 prog]);
                        Remain = repmat('-', [1 10-prog]);
                        progBar = [Bar, Remain];
                        fprintf('%s\n', progBar)
                    end
    
            
              
                
                
                end 
                
            end
            
            CoRippleRatesRegion = nan(length(regions));

            chAregions = rippleStats.chanLabels(chA_all);
            chBregions = rippleStats.chanLabels(chB_all);
            for iRa = 1:length(regions)
                for iRb = 1:length(regions)
                    iiA = contains(chAregions, regions{iRa});
                    iiB = contains(chBregions, regions{iRb});

                    coRAB = CoRippleRates(iiA, iiB);

                    if iiA == iiB; coRAB(logical(eye(sum(iiA)))) = []; end

                    CoRippleRatesRegion(iRa, iRb) = mean(coRAB(:), 'omitnan');

                end
            end
            filename = sprintf('%s_coRipple-microOnly-fast_Orig.mat', subject);
            save(fullfile(exportDir, filename), 'CoRippleRates', 'CoRippleN', 'euDistances', 'CoRippleRatesRegion', '-v7.3')
         otherwise
            f = contains(LFPfilesMacr, subject);
            macrObj = load(fullfile(dataDirectory, LFPfilesMacr{f}), 'chan_coords');

            f = contains(bpFiles, subject);
            bObj = load(fullfile(dataDirectory, '../out', bpFiles{f}), 'bpCoords');

            figure(1);
            filename = sprintf('%s_coRipple-microOnly-fast_Orig.mat', subject);
            load(fullfile(exportDir, filename));
            subplot(1,2,1)
            pl = scatter(euDistances(euDistances > 0), CoRippleRates(euDistances > 0), 'o'); hold on; ylim([0 0.45]);
            pl.MarkerFaceColor = pl.CData;

            pl.SizeData = 20;
            alpha(0.02);


            subplot(1,2,2)
            pl = plot3(macrObj.chan_coords(:,1), macrObj.chan_coords(:,2), macrObj.chan_coords(:,3), '.'); hold on;
            pl.Color = [0.7 0.7 0.7];

            pl = plot3(bObj.bpCoords(:,1), bObj.bpCoords(:,2), bObj.bpCoords(:,3), 'r.'); hold on;

            CoRippleRatesRegionAll(:,:,subj) = CoRippleRatesRegion;

            savepdf(gcf, fullfile(exportDirPDF, 'CoRipDistance_micro_orig.pdf'))

            figure(2);
            imagesc(mean(CoRippleRatesRegionAll, 3, 'omitnan'), [0 0.2])
            ax = gca;
            ax.XTick = 1:length(regions);
            ax.YTick = 1:length(regions);
            ax.XTickLabel = regions;
            ax.YTickLabel = regions;
            savepdf(gcf, fullfile(exportDirPDF, 'CoRipRegion_micro_orig.pdf'))


    end

    

    pause(0.5)






end


figure; 
boundedline(0:99, mean(rippleFR, 'omitnan'), std(rippleFR, 'omitnan')/sqrt(length(subj_list_full)));
xlim([0 20])

figure; 

imagesc(CoRippleRatesRegion, [0 0.1])
ax = gca;
ax.XTick = 1:length(regions);
ax.YTick = 1:length(regions);
ax.XTickLabel = regions;
ax.YTickLabel = regions;

















