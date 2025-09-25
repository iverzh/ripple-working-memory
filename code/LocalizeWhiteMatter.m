


close all 
% clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))

%%

subj_list_full = {'P41CS_R1', 'P41CS_R2', 'P42CS_R1', 'P42CS_R2', 'P43CS_R1', 'P43CS_R2', 'P44CS_R1', 'P47CS_R1', 'P47CS_R2', ...
                  'P49CS_R1', 'P49CS_R2', 'P51CS_R1', 'P51CS_R2', 'P53CS_R1', 'P53CS_R2', 'P54CS_R1', 'P54CS_R2', ...
                  'P55CS_R1', 'P55CS_R2', 'P56CS_R1', 'P56CS_R2', 'P57CS_R1', 'P57CS_R2', 'P58CS_R1', 'P60CS_R1', ...
                  'P62CS_R1', 'P62CS_R2'};

matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';
exportDir =     '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures/WM_localization';
if ~isfolder(exportDir); mkdir(exportDir); end
recordingState = 'wake';
location = 'NC';
figure('Position', [1219 728 1734 845])
for subj = 1:length(subj_list_full)
    subject = subj_list_full{subj};
    modifier = 'macro';
    tag = [recordingState,'_',location,'_',modifier];
    filename = sprintf('%s_ripple_stats_%s.mat', subject, tag);
    load(fullfile(matExportFolder, filename))
    
    lab = rippleStats.chanLabels;
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
        idx = find(nMsk,1,'first');
        alph{iL} = lab{iL}(1:idx-1);
        numr{iL} = str2double(lab{iL}(idx:end));
      end
    end
    
    shankLocations = unique(alph);
    for iL = 1:length(shankLocations)
        iiL = contains(alph, shankLocations{iL});
        XX = cellfun(@(X) mean(X), rippleStats.rippleAmp(iiL));
        subplot(3,4,iL);
        plot(XX, 'o-')
        ax = gca;
        ax.XTick = 1:sum(iiL);
        ax.XTickLabel = rippleStats.chanLabels(iiL);
        hline(2)
    end

    sgtitle(strrep(subject, '_',' '))

    fig = gcf;
    fig.Color = 'w';

%     WMii = 

    savepdf(gcf, fullfile(exportDir, sprintf('%s_rippleAmplitudes.pdf', subject)))






end




























































