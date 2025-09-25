



close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))

%%




subj_list_full = {'T11'};
study = 'BG';
state = 'NREM'; %NREM, wake or Simon
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';


exportDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles/%s', study);
if ~isfolder(exportDir); mkdir(exportDir); end
fs = 1e3;
cd(exportDir); pause(0.5);


%% Preprocess co ripple and control spike times

useRb = false; %flag whether you want to analyze the ripple band phase
cUnit = 0; 
saveLFP = false;
padLength = 0; %sec
for coRipPadMS = 250 %[1000, 200, 100] %co ripple padding time for the CCGs
    coRipPad = round(coRipPadMS * fs / 1e3);
    unitCCGstep = fs / 1e3;
    histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
    for subj = 1 %loop through subjects (only 1 in brain gate data)
        subject = subj_list_full{subj};
        unitsM = LoadSpikeTimes(subject,'CellExplorer','medial');
        unitsL = LoadSpikeTimes(subject,'CellExplorer','lateral');

        for ii = 1:size(unitsL,1); unitsL{ii,1} = unitsL{ii,1} + 96; end

        rippMedial  = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat'])); %uncomment 
        rippLateral  = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));

%             rippMedial  = load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_medial.mat']));
%             rippLateral  = load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_lateral.mat']));

        units = [unitsM; unitsL];


        siteA = 1:96;
        siteB = 97:192;        

        chan_labels = [siteA, siteB];

        if strcmp(study, 'Sternberg')
            chan_labels = cellfun(@(X) str2double(X), chan_labels);
        end


        timesCCG = 1:rippMedial.rippleStats.recordingLength;


        load(sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/%s_stageMask_%s.mat', subject, subject, state))
        iiStage = find(StageMask_1kHz);
        iiStage = iiStage(1:1800e3); %update the indices here to select a subset of the recording stage so that the code runs faster.
        StageMask_1kHz(~ismember(1:length(StageMask_1kHz),iiStage)) = false;

        rippMaskM = zeros(length(rippMedial.rippleStats.chanLabels), rippMedial.rippleStats.recordingLength);
        rippMaskL = zeros(length(rippLateral.rippleStats.chanLabels), rippLateral.rippleStats.recordingLength);
        for chRipp = 1:size(rippMaskM,1)  %loop through channels across both UAs
            if ~isempty(rippMedial.rippleStats.window{chRipp})
                iS = rippMedial.rippleStats.window{chRipp}(:,1);
                iE = rippMedial.rippleStats.window{chRipp}(:,2);

                for ii = 1:length(iE)
                    rippMaskM(chRipp,iS(ii):iE(ii)) = 1;
                end
            end

            if ~isempty(rippLateral.rippleStats.window{chRipp})
                iS = rippLateral.rippleStats.window{chRipp}(:,1);
                iE = rippLateral.rippleStats.window{chRipp}(:,2);

                for ii = 1:length(iE)
                    rippMaskL(chRipp,iS(ii):iE(ii)) = 1;
                end
            end

        end

        rippMask = [rippMaskM; rippMaskL];
        rippMask(:,~StageMask_1kHz) = 0;
        clear rippMaskL rippMaskM;





        for uAi = 1:size(units, 1) %loop through A units
            coRipSpike = cell(length(units),length(units));
            coRipDur = cell(length(units),length(units));
            controlDur = cell(length(units),length(units));
            noRipDur = cell(length(units),length(units));
            noRipAspike = cell(length(units),length(units));
            noRipSpike = cell(length(units),length(units));
            interactionType = cell(length(units),length(units));
            ccgNoR = cell(length(units),length(units),2);
            ccgCoR = cell(length(units),length(units));
            ccgNoRcontrol = cell(length(units),length(units));
            spikeData = cell(length(units),length(units));

 
            for uBi = 1:size(units, 1) % loop through B units 

                uA = uAi; uB = uBi;

                chA = units{uA, 1};
                chB = units{uB, 1};
                typeA = units{uA,3}; typeB = units{uB,3};
                if isempty(chA) || isempty(chB); continue; end





                if chA ~= chB && any(strcmp(typeA,{'pyr', 'int'})) && any(strcmp(typeB,{'pyr', 'int'}))

                    uTa = (units{uA,2}) * rippLateral.rippleStats.fs; uTa = uTa(uTa >= 0.5);
                    uTb = (units{uB,2}) * rippLateral.rippleStats.fs; uTb = uTb(uTb >= 0.5);


                    coRipMask = rippMask(chA, :) & rippMask(chB, :) & StageMask_1kHz;
                    noRipMask = ~rippMask(chA, :) & ~rippMask(chB, :) & StageMask_1kHz;



                    uTcoRa = uTa(coRipMask(round(uTa)));
                    uTcoRb = uTb(coRipMask(round(uTb)));

                    uTnoRa = uTa(noRipMask(round(uTa)));
                    uTnoRb = uTb(noRipMask(round(uTb)));

                    controlMask = zeros(size(coRipMask));

                    if ~isempty(uTcoRa) && ~isempty(uTcoRb)


                        bnds = mask2bounds(coRipMask);
                        dur = [bnds(:,2) - bnds(:,1)] / rippLateral.rippleStats.fs * 1e3;
                        bnds(dur < 25, :) = [];
                        dur(dur < 25) = [];

                        if size(bnds,1) < 2; continue; end

                        coSpikeCount = nan(2,length(bnds));
                        spikeDataAB = cell(7,length(bnds));
                        coSpikePhaseAB = cell(1,length(bnds));
                        coSpikeLFPAB = cell(1,length(uTcoRa)); cnt = 1;
                        noSpikeCount = nan(2,length(bnds));

                        rbA = nan(length(bnds), 201);
                        rbB = nan(length(bnds), 201);
                        bndsBaselineAll = nan(size(bnds)); 
                        tme = 0;
                        for ibnd  = 1:length(bnds) % loop through co ripple bounds
                            aSpikesCo = sum(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                            bSpikesCo = sum(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));

                            % find duration match control periods with
                            % no ripple activity
                            shift = round(2e3 * rippLateral.rippleStats.fs / 1e3);
                            bndsBaseline = bnds(ibnd,:) - shift;
                            if any(bndsBaseline <= 0); bndsBaselineAll(ibnd, :) = nan; continue; end % co-ripple is too close to beginning of recording

                            loopCount = 0;
                            while ~all(noRipMask(bndsBaseline(1):bndsBaseline(2)))  
                                bndsBaseline = bndsBaseline - 100;
                                loopCount = loopCount + 1;
                                if any(bndsBaseline < 1) || loopCount > 300  
                                    break
                                end
                            end


%                                 
                            % store relevant outputs 
                            aSpikesNo = sum(uTnoRa >= bndsBaseline(1)   &  uTnoRa <= bndsBaseline(2)); % # spikes in A for each no ripple control 
                            bSpikesNo = sum(uTnoRb >= bndsBaseline(1)   &  uTnoRb <= bndsBaseline(2)); % # spikes in B for each no ripple control 

                            uTcoRaBnd = uTa(uTa >= bnds(ibnd,1)  &  uTa <= bnds(ibnd,2));
                            coSpikeCount(:, ibnd)     = [aSpikesCo, bSpikesCo];  % # spikes in A and # spikes in B for each co ripple  
                            spikeDataAB{1, ibnd}   = uTa(uTa >= (bnds(ibnd,1)-coRipPad)  &  uTa <= (bnds(ibnd,2)+coRipPad)); % A spikes +- padding window
                            spikeDataAB{2, ibnd}   = uTcoRaBnd; % A spikes only during co ripple
                            spikeDataAB{3, ibnd}   = uTb(uTb >= (bnds(ibnd,1)-coRipPad)  &  uTb <= (bnds(ibnd,2)+coRipPad)); % B spikes +- padding window

                            uTcontrolRaBnd = uTa(uTa >= bndsBaseline(1)  &  uTa <=  bndsBaseli.ne(2) );
                            spikeDataAB{4, ibnd}   = uTcontrolRaBnd; % A spikes during duration matched no ripple control period
                            spikeDataAB{5, ibnd}   = uTb(uTb >= ( bndsBaseline(1) -coRipPad)  &  uTb <= (bndsBaseline(2)+coRipPad)); % B spikes during no ripple control period +- padding window

                            if useRB

                                spikeDataAB{6, ibnd}   = RBphaseAll(chA,round(spikeDataAB{2, ibnd})); %local ripple phase during A spikes
                                spikeDataAB{7, ibnd}   = RBphaseAll(chB,round(spikeDataAB{3, ibnd})); %local ripple phase during B spikes

                                coSpikePhaseAB{ibnd}(1,:) = RBphaseAll(chA,bnds(ibnd,1) :bnds(ibnd,2) ); % A phase during co-ripples
                                coSpikePhaseAB{ibnd}(2,:) = RBphaseAll(chB,bnds(ibnd,1) :bnds(ibnd,2) ); % B phase during co-ripples
                            end        

                            noSpikeCount(:, ibnd) = [aSpikesNo bSpikesNo]; % # spikes in A and # spikes in B for each no ripple control 
                            bndsBaselineAll(ibnd, :) = bndsBaseline;



                        end
%                             
%                             




                        coRipSpike{uA,uB} = coSpikeCount;
                        dur(isnan(bndsBaselineAll(:,1))) = nan;
                        coRipDur{uA,uB} = dur;
                        controlDur{uA,uB} = bndsBaselineAll(:,2) - bndsBaselineAll(:,1);
                        noRipDur{uA,uB} = sum(noRipMask);
                        noRipAspike{uA,uB} = length(uTnoRa);
                        noRipSpike{uA,uB} = noSpikeCount;

                        coRipMask = bounds2mask(bnds, length(coRipMask));
                        uTcoRa = uTa(coRipMask(round(uTa)));
                        coRipMask = bounds2mask(bnds, length(coRipMask), coRipPad);
                        uTcoRb = uTb(coRipMask(round(uTb)));


                        spikeData{uA, uB} = spikeDataAB;

                        ccgNoR{uA,uB,1} = uTnoRa;
                        ccgNoR{uA,uB,2} = uTnoRb;




                    end
                end

                clc
            fprintf('... %s co firing unit [%i %i] / %i ... \n co ripple padding %03i\n', subject, uAi, uBi, length(units), coRipPad)
            end


            uCh = units(:,1);
            uLoc = units(:,end); %cellfun(@(X,Y) [X '_' n
            filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_u%i_%s-cut.mat', subject, coRipPad, uAi, state);
            save(fullfile(exportDir, filename), 'coRipSpike', 'coRipDur', 'controlDur', 'noRipSpike', ...
                                            'noRipDur', 'noRipAspike', ...
                                            'ccgCoR','ccgNoR', 'ccgNoRcontrol', 'uCh', 'uLoc',...
                                            'spikeData', '-v7.3') % save data for each A unit


             cUnit=cUnit+1;
        end
%             

        uCh = units(:,1);
        uLoc = units(:,end); 




    end
end


%% Compute CCGs and prepare for plotting 



fs = 1e3; %1e3, 400;

for ang = 2 % keep this for now, but it does not do anything in this version of the script 
    tic





    phLagsAll = nan(1, 1e6); cLags = 1;
    cAP =1;
    coRipPad = 0.25 * fs;
    durTrough = 0; durPeak = 0;
    unitCCGstep = fs / 1e3;
    histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
    coFd = 20;
    minDur = 25;


    lfpA =  nan(1e6, length(-100:100)); lfpB =  nan(1e6, length(-100:100)); 
    rbA =  nan(1e6, length(-100:100)); rbB = nan(1e6, length(-100:100));
    phLag =  nan(1,1e6);
    ccgAllcoR = nan(1e6, length(histEdges)-1);
    ccgAllcoRpeak =  nan(1e6, length(histEdges)-1);
    ccgAllcoRtrough =  nan(1e6, length(histEdges)-1);
    ccgAllnoRAll = nan(1e6, length(histEdges)-1);
    ccgAllnoRepoch =  nan(1e6, length(histEdges)-1);
    durAllcoR =  nan(1,1e6);
    durAllcoRpeak =  nan(1,1e6);
    durAllcoRtrough =  nan(1,1e6);
    durAllnoRAll =  nan(1,1e6);
    durAllnoRepoch =  nan(1,1e6);
    nAspikeAllcoR =  nan(1,1e6);
    nAspikeAllcoRpeak =  nan(1,1e6);
    nAspikeAllcoRtrough =  nan(1,1e6);
    nAspikeAllnoRAll =  nan(1,1e6);
    nAspikeAllnoRepoch =  nan(1,1e6);
    compStringAll_LFP = cell(1,1e6);
    compStringAll_phase = cell(1,1e6);
    pairType = cell(1e6,2);
    ccgCrossArray = nan(1,1e6);
    aSpPhase = nan(1,1e6);
    bSpPhase = nan(1,1e6);
    Nphase = zeros(1, length(histEdges)-1);
    countLFP = 0;countPhase = 0;countCCG = 1;
    for subj = 1:length(subj_list_full) %[1:17 20:length(subj_list_full)] % [1:15, 18:length(subj_list_full)]
        subject = subj_list_full{subj};
        tractDistance =  []; compString = [];
        phLagsSubj = []; 


        flist = dir(fullfile(exportDir, '*.mat'));
        flist = {flist.name}';



        unitsM = LoadSpikeTimes(subject,'CellExplorer','medial');
        unitsL = LoadSpikeTimes(subject,'CellExplorer','lateral');

        for ii = 1:size(unitsL,1); unitsL{ii,1} = unitsL{ii,1} + 96; end
        units = [unitsM; unitsL];
        uCh = units(:,1);

        SUii = find(strcmp(units(:,3),'pyr') | strcmp(units(:,3),'int'));

        cSubj = 0;
         for uAi = 1:length(SUii)
            uA = SUii(uAi);
            filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_IISremv_SU_u%i_%s.mat', subject,fs/4, uA, state); %update this to include output of previous section
            load(fullfile(exportDir, filename)); 
            for uBi = 1:length(SUii)
                uB = SUii(uBi);

                typeA = units(uA,3);
                typeB = units(uB,3);



                coRip = coRipSpike{uA,uB};
                noRip = noRipSpike{uA,uB};
                dur = (sum(coRipDur{uA,uB})); %+ 2*0*length(coRipDur{uA,uB})) /fs;
                durNo = noRipDur{uA,uB} / fs;

                if dur > 0 && uCh{uA} ~= uCh{uB} 
                    sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                    nCoSpikeCo = size(sp,2); 

                    sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                    nCoSpikeNo = size(sp,2);


                    if (uCh{uA} <= 96 && uCh{uB} > 96) || (uCh{uA} > 96 && uCh{uB} <= 96)
                    d = 1; else
                    d = 0;
                    end


                    if d >= 0 

                        ripDur = coRipDur{uA,uB}' >= minDur;



                        spMat = coSpikeTimes{uA, uB};
    %                     ;
                        uTa = cell2mat(spMat(2,ripDur)'); uTb = cell2mat(spMat(3,ripDur)');
                        if ~isempty(uTa) && ~isempty(uTb)
                            s = uTa - uTb';
                            cAP = cAP + length(uTa);
                            Nall = histcounts(s, histEdges);
                        else
                            Nall = zeros(1,length(histEdges)-1);
                        end


                        durAll = sum(coRipDur{uA,uB}(ripDur)) / fs;
                        nAspikeAllcoR(countCCG) =  length(uTa);


                        uTa = cell2mat(spMat(4,:)'); uTb = cell2mat(spMat(5,:)');
                        if ~isempty(uTa) && ~isempty(uTb)
                            s = uTa - uTb';
                            Ncontrol = histcounts(s, histEdges);
                        else
                            Ncontrol = zeros(1,length(histEdges)-1);
                        end
                        nAspikeAllnoRepoch(countCCG) = length(uTa);


                        ccgAllcoR(countCCG, :) = Nall;
                        durAllcoR(countCCG) = (durAll);



                        durAllnoRAll(countCCG) = durNo;
%                             nAspikeAllnoRAll(countCCG) = length(uTnoRa);

                        ccgAllnoRepoch(countCCG, :) = Ncontrol;
                        durAllnoRepoch(countCCG) = durAll;
                        ccgCrossArray(countCCG) = d;
                        pairType(countCCG,1) = typeA;
                        pairType(countCCG,2) = typeB;


                        chA = uCh{uA}; chB = uCh{uB};

                        countCCG = countCCG + 1;


                    end



                end
            end
            clc
            fprintf('unit %i done...\n', uA)
         end

        fprintf('%s   angle %i coFire ver %i\n', subject, ang, v)



    end
    disp('computing CCGs complete.')
    toc


end


%% Plot CCGs

% close all
clr = brewermap(10,'Paired');
smWin = 9; %*fs/1e3;
smoothMethod = 'gaussian';
tag = 'bmovie_perSpike';
% tag = 'Sternberg_perSpike';
% tag = 'Sternberg';
% tag = 'bothStudies';

xtimes = -fs/4:fs/4;
xl = [-100 100];

ii = 3;
figure('Position',  [846 419 1158 902]);

typeA = {'pyr' , 'pyr', 'int', 'int'};
typeB = {'pyr' , 'int', 'pyr', 'int'};

ccgCrossArray(isnan(ccgCrossArray)) = [];
IIdel = cell2mat(cellfun(@(X) isempty(X), pairType,  'UniformOutput', false));
pairType(IIdel(:,1), :) = [];
for iT = 1:4
    subplot(2,2,iT)
    keep = logical(ccgCrossArray); % ~logical(ccgCrossArray) for within array cell pairs


    xT = xtimes >= xl(1) & xtimes <= xl(2);
     
    keep = keep & strcmp(pairType(:,1),typeA{iT})' &  strcmp(pairType(:,2),typeB{iT})';

    D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
    D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
    D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
    D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';

%     D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)' * 1e3;
%     D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)' * 1e3;
%     D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)' * 1e3;
%     D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)' * 1e3;   




    % 
    % 
    eval(sprintf('D = D%i;', ii))
    D = smoothdata(D, 2, smoothMethod, smWin);
    Dnull = smoothdata(D4, 2, smoothMethod, smWin);


    [pl1, pa1] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(D, 'omitnan'), std(D, 'omitnan')/sqrt(size(D,1))); hold on; 
    pl1.Color = clr(8,:);
    pa1.FaceColor = clr(8,:);
    pl1.LineWidth = 2;
    pa1.FaceAlpha = 0.3;

    [pl2, pa2] = boundedline(-coRipPad:unitCCGstep:coRipPad, mean(Dnull, 'omitnan'), std(Dnull, 'omitnan')/sqrt(size(Dnull,1))); hold on; 
    pl2.Color = clr(1,:);
    pa2.FaceColor = clr(1,:);
    pl2.LineWidth = 2;
    pa2.FaceAlpha = 0.3;
    vline([-99:11:99]*fs/1e3);
    vl = vline(0); vl.Color = 'k'; vl.LineStyle = '-'; 

    xlim((xl*fs)/1e3)
    %             ylim(yl)
    ylabel('firing rate')
    xlabel('time [ms]')
    %             % grid on;
    fig = gcf;
    fig.Color = 'w';

    title(sprintf('%s <--> %s',typeA{iT}, typeB{iT}))
end





