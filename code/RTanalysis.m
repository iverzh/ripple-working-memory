


close all
clc

figure('Position',[669 805 210 round(729/5)]);
% plot_condition = true(1, length(correct_trials_all));
plot_condition_2 = load_trials_all == 3 & respLatency > median(respLatency(load_trials_all == 3));
plot_condition = load_trials_all == 3 & respLatency <= median(respLatency(load_trials_all == 3));
condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');

baseln = sum(nChRegions,2); % .* (nChRegions(:,iR)-1)/2;
    
baseln_1 = repmat(baseln(plot_condition), [1 size(whole_trial_all,2)]);
plot_data    = smoothdata(whole_trial_all(plot_condition,:)./baseln_1, 2, 'gaussian',100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition));
plot_data_mu(isnan(whole_trial_all(condTrial1,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial1,:))) = nan;

[bl1, bf] = boundedline([1:length(plot_data_mu)]-taskMarkers(end), plot_data_mu, plot_data_sem, 'r', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
bl1.Color = [1 0 0];
bl1.LineWidth = 1;
bf.FaceColor = [1 0 0];
bl1.MarkerFaceColor = bl1.Color;
bf.FaceAlpha = 0.3;

baseln_3 = repmat(baseln(plot_condition_2), [1 size(whole_trial_all,2)]);
plot_data = smoothdata(whole_trial_all(plot_condition_2,:)./baseln_3, 2, 'gaussian',100);
plot_data_mu = mean(plot_data, 'omitnan');
plot_data_sem = std(plot_data, 'omitnan')/sqrt(sum(plot_condition_2));
plot_data_mu(isnan(whole_trial_all(condTrial2,:)))  = nan;
plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
[bl2, bf] = boundedline([1:length(plot_data_mu)]-taskMarkers(end), plot_data_mu, plot_data_sem, 'k--', 'nan', 'gap'); hold on;
bf.FaceAlpha = 0.7;
bl2.Color = 0.75*[1 0 0];
bl2.LineWidth = 1;
bf.FaceColor = 0.75*[1 0 0];
bl2.MarkerFaceColor = bl2.Color;
bf.FaceAlpha = 0.3;


baseln = repmat(baseln, [1 1e3]);
hl = hline(mean(whole_trial_all(:,1:1e3)./baseln, 'all', 'omitnan'));
hl.Color = 'k';
hl.LineStyle = '-';
hl.LineWidth = 0.5;
xlim([-200 1.5e3])

ax = gca;
ax.FontSize = 8;

vl = vline([0]); 
for iV = 1:length(vl); vl(iV).LineWidth = 1.5; end

fig = gcf;
fig.Color = 'w';
% ylabel('co-ripples')
% xlabel('time from probe [ms]')
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTask_allChannels_RT_%s.pdf', tag)))


dat = mean(whole_trial_all(plot_condition | plot_condition_2,taskMarkers(end):taskMarkers(end)+(1e3)), 2) ./ baseln(plot_condition | plot_condition_2);


% tab = table(subjID(~isnan(dat)), cond(~isnan(dat)), dat(~isnan(dat)), ...
%       'VariableNames', {'subject',   'condition',                 'rates'});
% lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
% pLME = lme.Coefficients(2,6); c = c +1;


regionColors =  brewermap(12, 'Dark2');


% plot_condition_2 = load_trials_all == 3;
% plot_condition = load_trials_all == 1;

condTrial1 = find(plot_condition, 1, 'first');
condTrial2 = find(plot_condition_2, 1, 'first');
% plot_condition = logical(novel_trials_all(~isnan(correct_trials_all)));


figure('Position',[976 401 210 729]);

blall = [];
binSz = 1000;
bins = 1:binSz:size(whole_trials_region_all,2);
times = 1:size(whole_trials_region_all,2);
for iR = 4:5 %1:5 %1:length(regions)
    subplot(5,1,iR)

    nTrial = arrayfun(@(X) sum(subjID == X), 1:max(subjID));
    baseln = nChRegions(:,iR); % .* (nChRegions(:,iR)-1)/2;
    
    plot_data_1 = whole_trials_region_all(plot_condition_2,:,iR);    
    baseln_1 = repmat(baseln(plot_condition_2), [1 size(plot_data_1,2)]);
%     plot_data_1 = (plot_data_1 - baseln_1) ./ baseln_1;
    plot_data_1 = plot_data_1 ./ baseln_1;
    
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));
    plot_data_mu_1 = smoothdata(plot_data_mu_1, 2, 'gaussian',100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',100);
    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline([1:length(plot_data_1)]-taskMarkers(end), plot_data_mu_1, plot_data_sem, 'k--', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl2.Color = 0.75*regionColors(iR,:);
    bl2.LineWidth = 1;
    bf.FaceColor = 0.75*regionColors(iR,:);
    bl2.MarkerFaceColor = bl2.Color;
    bf.FaceAlpha = 0.3;    
    
    plot_data_3    = whole_trials_region_all(plot_condition,:,iR);
    baseln_3 = repmat(baseln(plot_condition), [1 size(plot_data_3,2)]);
    
%     plot_data_3 = (plot_data_3 - baseln_3) ./ baseln_3;
    plot_data_3 = plot_data_3 ./ baseln_3;
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));
    plot_data_mu_3 = smoothdata(plot_data_mu_3, 2, 'gaussian',100);
    plot_data_sem = smoothdata(plot_data_sem, 2, 'gaussian',100);
    
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline([1:length(plot_data_3)]-taskMarkers(end), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bl1.LineWidth = 1;
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
   
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end)+250:taskMarkers(end)+1.0e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end)+250:taskMarkers(end)+1.0e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    
    
    
    subjCond = subjID(plot_condition | plot_condition_2);
    cond = zeros(length(plot_condition),1); cond(plot_condition) = 2; cond(plot_condition_2) = 1;
    cond(cond == 0) = [];
    
    c = 1;
    dat = mean(whole_trials_region_all(plot_condition | plot_condition_2,taskMarkers(end):taskMarkers(end)+(1e3),iR), 2);% ./ baseln(plot_condition | plot_condition_2);


    tab = table(subjID(~isnan(dat)), cond(~isnan(dat)), dat(~isnan(dat)), ...
          'VariableNames', {'subject',   'condition',                 'rates'});
    lme = fitlme(tab, 'rates ~ condition + (1|subject)');   
    pLME = lme.Coefficients(2,6); c = c +1;
%         [~,~,p_perm] = statcond({dat(~isnan(dat) & cond == 2)',dat(~isnan(dat) & cond == 1)'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    

    fprintf('%s ripples maintenance p = %.2f\n', regions{iR}, pLME)



   
    xlim([-200 1.5e3])
    vl = vline([500 2600 4700 6800 9400]-500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
%     ylabel('rip. [bsl. norm.]')
%     xlabel('time from probe stim [ms]')
    
    ax = gca;
    ax.FontSize = 8;
end
fig = gcf; fig.Color = 'w';
savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegions_RT_%s.pdf', tag)))
%%

figure('Position',[669 402 140 729]);

blall = [];
binSz = 1000;
bins = 1:binSz:size(whole_trials_region_all,2);
times = 1:size(whole_trials_region_all,2);
for iR = 1:length(regions)
    subplot(5,1,iR)

    nTrial = arrayfun(@(X) sum(subjID == X), 1:max(subjID));
    baseln = arrayfun(@(X) mean(whole_trials_units_sm_region_all(subjID == X,1:1e3,iR), 'all', 'omitnan'), 1:max(subjID));
    baseln(baseln < 1e-4) = nan;
    baseln = repelem(baseln, nTrial)';
    
    plot_data_1 = whole_trials_units_sm_region_all(plot_condition_2,:,iR);    
    baseln_1 = repmat(baseln(plot_condition_2), [1 size(plot_data_1,2)]);
%     plot_data_1 = (plot_data_1 - baseln_1) ./ baseln_1;
%     plot_data_1 = plot_data_1 ./ baseln_1;
    
%     plot_data(isnan(whole_trials_region_load1_all(:,:,iR))) = nan;
    plot_data_mu_1 = mean(plot_data_1, 'omitnan');
    plot_data_sem = std(plot_data_1, 'omitnan')/sqrt(sum(~isnan(plot_data_1(:,1))));

    plot_data_mu_1(isnan(whole_trial_all(condTrial2,:)))  = nan;
    plot_data_sem(isnan(whole_trial_all(condTrial2,:))) = nan;
    [bl2, bf] = boundedline([1:length(plot_data_1)]-taskMarkers(end), plot_data_mu_1, plot_data_sem, 'k', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl2.Color = 0.5*regionColors(iR,:);
    bf.FaceColor = 0.5*regionColors(iR,:);
    bl2.MarkerFaceColor = bl2.Color;
    bf.FaceAlpha = 0.3;    
    
    plot_data_3    = whole_trials_units_sm_region_all(plot_condition,:,iR);
    baseln_3 = repmat(baseln(plot_condition), [1 size(plot_data_3,2)]);
    
%     plot_data_3 = (plot_data_3 - baseln_3) ./ baseln_3;
%     plot_data_3 = plot_data_3 ./ baseln_3;
%     plot_data(isnan(whole_trials_region_load3_all(:,:,iR))) = nan;
    plot_data_mu_3 = mean(plot_data_3, 'omitnan');
    plot_data_sem = std(plot_data_3, 'omitnan')/sqrt(sum(~isnan(plot_data_3(:,1))));

    
    plot_data_mu_3(isnan((whole_trial_all(condTrial1,:))))  = nan;
    plot_data_sem(isnan((whole_trial_all(condTrial1,:)))) = nan;

    [bl1, bf] = boundedline([1:length(plot_data_3)]-taskMarkers(end), plot_data_mu_3, plot_data_sem, 'r', 'nan', 'gap'); hold on;
    bf.FaceAlpha = 0.7;
    bl1.Color = regionColors(iR,:);
    bf.FaceColor = regionColors(iR,:);
    bl1.MarkerFaceColor = bl1.Color;
    bf.FaceAlpha = 0.3;
    
   
    
    hl = hline(mean([plot_data_mu_1(1:1e3), plot_data_mu_3(1:1e3)]));
    hl.Color = 'k';
    hl.LineStyle = '-';
    hl.LineWidth = 0.5;
    
    probeResp_1 = trapz(plot_data_1(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(plot_data_3(:,taskMarkers(end):taskMarkers(end)+1.0e3), 2); probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s ripples response p = %.2f\n', regions{iR}, p_perm)
    
    probeResp_1 = trapz(fillmissing(plot_data_1(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
    probeResp_1(isnan(probeResp_1)) = [];
    probeResp_3 = trapz(fillmissing(plot_data_3(:,taskMarkers(end-1):taskMarkers(end)), 'linear', 2), 2); 
    probeResp_3(isnan(probeResp_3)) = [];
    [~,~,p_perm] = statcond({probeResp_3',probeResp_1'},'paired','off','method', 'perm', 'naccu', 10000, 'verbose','off');
    
    fprintf('%s ripples maintenance p = %.2f\n', regions{iR}, p_perm)



    for iB = 1:length(bins)-1
        
        
%         data = recog_trials_all_region(plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(plot_condition)), 'o'); hold on;
%         bl1.Color = regionColors(iR,:);
        bl1.MarkerFaceColor = bl1.Color;
%     
%         bf.FaceColor = regionColors(iR,:);
% 
%         data = recog_trials_all_region(~plot_condition,times >= bins(iB) & times < bins(iB+1), iR);
%         bl1 = errorbar(bins(iB) + (binSz/2), mean(data(:), 'omitnan'), std(data(:), 'omitnan')/sqrt(sum(~plot_condition)), 'ko');
%         hold on;
%         bl1.MarkerFaceColor = bl1.Color;

    

    
    end
    


%     title(regions{iR})

%     legend([bl1, bl2], 'load3', 'load1', 'location', 'best')
%     legend([bl1, bl2], 'probe in', 'probe out', 'location', 'best')


%     blall(iR) =bl;
%     xlim(Lim)
%     ylim([0.3 1])
%     xlim([max(taskMarkers)-1e3 max(taskMarkers) + 1e3])
    xlim([-200 1.5e3])
    vl = vline([500 2600 4700 6800 9400]-500); 
    for iV = 1:length(vl); vl(iV).LineWidth = 1.0; end
    
%     ylabel('rip. [bsl. norm.]')
%     xlabel('time from probe stim [ms]')
    
    ax = gca;
    ax.FontSize = 10;
end
fig = gcf; fig.Color = 'w';
% savepdf(gcf, fullfile(exportDirFigs, sprintf('coRipTaskWithinRegions_%s.pdf', tag)))