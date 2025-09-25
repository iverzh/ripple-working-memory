



close all 
clc
clear


addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/UtahArray'))
addpath(genpath('/space/seh10/6/halgdev/projects/iverzh/ripples/code/util'))
addpath(genpath('/space/seh8/1/halgdev/projects/cdickey/packages'))

%%




subj_list_full = {'T11'};
study = 'BG';
state = 'wake'; %NREM %Simon
matExportFolder = '/space/seh10/6/halgdev/projects/iverzh/ripples/matFiles';


exportDir = sprintf('/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/processedResults/coFire/subjectFiles/%s', study);
if ~isfolder(exportDir); mkdir(exportDir); end
fs = 1e3;
cd(exportDir); pause(0.5);


%%
nUnit = 0;
nUnitNoLFP = 0;
for v = 3
    cUnit = 0; 
    saveLFP = false;
    padLength = 0; %sec
    for coRipPadMS = 250 %[1000, 200, 100] 
        coRipPad = round(coRipPadMS * fs / 1e3);
        unitCCGstep = fs / 1e3;
        histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
        for subj = 1
            subject = subj_list_full{subj};
            unitsM = LoadSpikeTimes(subject,'CellExplorer','medial');
            unitsL = LoadSpikeTimes(subject,'CellExplorer','lateral');

            for ii = 1:size(unitsL,1); unitsL{ii,1} = unitsL{ii,1} + 96; end
            
%             rippMedial  = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_medial.mat']));
%             rippLateral  = load(fullfile(matExportFolder,[subject,'_ripple_stats_sleep_NC_NREMonly_lateral.mat']));
            
            rippMedial  = load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_medial.mat']));
            rippLateral  = load(fullfile(matExportFolder,[subject,'_ripple_stats_wake_NC_wakeonly_lateral.mat']));

            units = [unitsM; unitsL];


            siteA = 1:96;
            siteB = 97:192;        

            chan_labels = [siteA, siteB];
    
            if strcmp(study, 'Sternberg')
                chan_labels = cellfun(@(X) str2double(X), chan_labels);
            end

            
            timesCCG = 1:rippMedial.rippleStats.recordingLength;

            
            load(sprintf('/space/seh10/6/halgdev/projects/iverzh/data/UtahArrayData/%s/data_1kHz/%s_stageMask_%s.mat', subject, subject, state))



            rippMaskM = zeros(length(rippMedial.rippleStats.chanLabels), rippMedial.rippleStats.recordingLength);
            rippMaskL = zeros(length(rippLateral.rippleStats.chanLabels), rippLateral.rippleStats.recordingLength);
            for chRipp = 1:size(rippMaskM,1) 
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

   
        
        
            
            for uAi = 1:size(units, 1)
                coRipSpike = cell(length(units),length(units));
                coRipDur = cell(length(units),length(units));
                controlDur = cell(length(units),length(units));
                noRipDur = cell(length(units),length(units));
                noRipAspike = cell(length(units),length(units));
                noRipSpike = cell(length(units),length(units));
                interactionType = cell(length(units),length(units));
                coFireProbNoR = cell(length(units),length(units));
                ccgNoR = cell(length(units),length(units),2);
                ccgCoR = cell(length(units),length(units));
                ccgNoRcontrol = cell(length(units),length(units));
                coSpikeTimes = cell(length(units),length(units));

                
                for uBi = 1:size(units, 1)

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
                            coSpikeTimesAB = cell(7,length(bnds));
                            coSpikePhaseAB = cell(1,length(bnds));
                            coSpikeLFPAB = cell(1,length(uTcoRa)); cnt = 1;
                            noSpikeCount = nan(2,length(bnds));
                            
                            rbA = nan(length(bnds), 201);
                            rbB = nan(length(bnds), 201);
                            bndsBaselineAll = nan(size(bnds)); 
                            tme = 0;
                            for ibnd  = 1:length(bnds)
                                aSpikesCo = sum(uTcoRa >= bnds(ibnd,1)   &  uTcoRa <= bnds(ibnd,2));
                                bSpikesCo = sum(uTcoRb >= bnds(ibnd,1)   &  uTcoRb <= bnds(ibnd,2));
            
                                shift = round(1e3 * rippLateral.rippleStats.fs / 1e3);
                                bndsBaseline = bnds(ibnd,:) - shift;
                                if any(bndsBaseline <= 0); continue; end
            
                                loopCount = 0;
                                while ~all(noRipMask(bndsBaseline(1):bndsBaseline(2))) && loopCount < 300  
                                    bndsBaseline = bndsBaseline - 100;
                                    loopCount = loopCount + 1;
                                    if any(bndsBaseline < 1)
                                        loopCount = 300;
                                        break
                                    end
                                end
            

%                                 
            
                                aSpikesNo = sum(uTnoRa >= bndsBaseline(1)   &  uTnoRa <= bndsBaseline(2));
                                bSpikesNo = sum(uTnoRb >= bndsBaseline(1)   &  uTnoRb <= bndsBaseline(2));
                                
                                uTcoRaBnd = uTa(uTa >= bnds(ibnd,1)  &  uTa <= bnds(ibnd,2));
                                coSpikeCount(:, ibnd)     = [aSpikesCo, bSpikesCo];
                                coSpikeTimesAB{1, ibnd}   = uTa(uTa >= (bnds(ibnd,1)-coRipPad)  &  uTa <= (bnds(ibnd,2)+coRipPad));
                                coSpikeTimesAB{2, ibnd}   = uTcoRaBnd;
                                coSpikeTimesAB{3, ibnd}   = uTb(uTb >= (bnds(ibnd,1)-coRipPad)  &  uTb <= (bnds(ibnd,2)+coRipPad));
                                
                                uTcontrolRaBnd = uTa(uTa >= bndsBaseline(1)  &  uTa <=  bndsBaseline(2) );
                                coSpikeTimesAB{4, ibnd}   = uTcontrolRaBnd;
                                coSpikeTimesAB{5, ibnd}   = uTb(uTb >= ( bndsBaseline(1) -coRipPad)  &  uTb <= (bndsBaseline(2)+coRipPad));
                                
                                coSpikeTimesAB{6, ibnd}   = RBphaseAll(chA,round(coSpikeTimesAB{2, ibnd}));
                                coSpikeTimesAB{7, ibnd}   = RBphaseAll(chB,round(coSpikeTimesAB{3, ibnd}));
    
%                                 coSpikePhaseAB{ibnd}(1,:) = RBphaseAll(chA,bnds(ibnd,1) :bnds(ibnd,2) );
%                                 coSpikePhaseAB{ibnd}(2,:) = RBphaseAll(chB,bnds(ibnd,1) :bnds(ibnd,2) );
        
                                noSpikeCount(:, ibnd) = [aSpikesNo bSpikesNo];
                                bndsBaselineAll(ibnd, :) = bndsBaseline;
    
                                center = round(mean(bnds(ibnd, :)));

                                
            
                            end
%                             
%                             
                            
            
%                             controlMask = bounds2mask(bndsBaselineAll, length(coRipMask));
%                             uTcontrola = uTa(controlMask(round(uTa)));
%                             controlMask = bounds2mask(bndsBaselineAll, length(coRipMask), coRipPad);
%                             uTcontrolb = uTb(controlMask(round(uTb)));
%                 
%                             s = uTcontrola - uTcontrolb';
%                             Ncontrol = histcounts(s, histEdges);
            
                       
                            coRipSpike{uA,uB} = coSpikeCount;
                            coRipDur{uA,uB} = dur;
                            controlDur{uA,uB} = bndsBaselineAll(:,2) - bndsBaselineAll(:,1);
                            noRipDur{uA,uB} = sum(noRipMask);
                            noRipAspike{uA,uB} = length(uTnoRa);
                            noRipSpike{uA,uB} = noSpikeCount;
    
                            coRipMask = bounds2mask(bnds, length(coRipMask));
                            uTcoRa = uTa(coRipMask(round(uTa)));
                            coRipMask = bounds2mask(bnds, length(coRipMask), coRipPad);
                            uTcoRb = uTb(coRipMask(round(uTb)));

        
                            coSpikeTimes{uA, uB} = coSpikeTimesAB;

                            ccgNoR{uA,uB,1} = uTnoRa;
                            ccgNoR{uA,uB,2} = uTnoRb;
                            
                           
                            
                            
                        end
                    end
                    
                    clc
                fprintf('... %s co firing unit [%i %i] / %i ... \n co ripple padding %03i\n', subject, uAi, uBi, length(units), coRipPad)
                end
                
                
                uCh = units(:,1);
                uLoc = units(:,end); %cellfun(@(X,Y) [X '_' n
                filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_IISremv_SU_u%i_%s.mat', subject, coRipPad, uAi, state);
                save(fullfile(exportDir, filename), 'coRipSpike', 'coRipDur', 'controlDur', 'noRipSpike', ...
                                                'noRipDur', 'noRipAspike', 'coFireProbNoR', ...
                                                'ccgCoR','ccgNoR', 'ccgNoRcontrol', 'uCh', 'uLoc',...
                                                'coSpikeTimes', '-v7.3')
                
                
                 cUnit=cUnit+1;
            end
%             
            nCoSpikeCoAll = [];
            nCoSpikeNoAll = [];
            
            distanceAll = [];
            c = 1;
            
            nCoSpikeNo = [];
            nCoSpikeCo = [];
            %     dur = 0;
            for uAi = 1:size(units,1)
                for uBi = 1:size(units,1)
            %         uA = units(uAi); uB = units(uBi);
                    uA = uAi; uB = uBi;
    %                 chA = find(chan_labels == units{uA,1}); chB = find(chan_labels == units{uB,1});
                    chA = find(strcmp(chan_labels, units{uA,end})); chB = find(strcmp(chan_labels, units{uB,end}));
    
                    if isempty(chA) || isempty(chB); continue; end
                    d = pdist([micrObj.chan_coords(chA,:); micrObj.chan_coords(chB,:)]);
        
                    
                    coRip = coRipSpike{uA,uB};
                    noRip = noRipSpike{uA,uB};
                    dur = sum(coRipDur{uA,uB}) / rippLateral.rippleStats.fs;
        
                    if ~isempty(coRip) 
                        sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                        nCoSpikeCo = size(sp,2); %/dur; % / (dur);
                        
                        sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0); 
                        nCoSpikeNo = size(sp,2); %/dur; % / (dur);
        %                     dur = dur + (sum(coRipDur{uA,uB}) / rippleStats.fs);
            
                        nCoSpikeCoAll = [nCoSpikeCoAll, nCoSpikeCo];
                        nCoSpikeNoAll = [nCoSpikeNoAll, nCoSpikeNo];
                        distanceAll = [distanceAll, d];
                    end
                    
                        
            
                end
            end
            uCh = units(:,1);
            uLoc = units(:,end); 
            
            filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_IISremv_SU_v%i.mat', subject, coRipPad, 6);
            save(fullfile(exportDir, filename), 'nCoSpikeCoAll', 'nCoSpikeNoAll', 'distanceAll', ...
                                                'coRipSpike', 'coRipDur', 'controlDur', 'noRipSpike', ...
                                                'noRipDur', 'noRipAspike', 'coFireProbNoR', ...
                                                'ccgNoR', 'ccgCoR','ccgNoRcontrol', 'uCh', 'uLoc',...
                                                'coSpikeTimes', 'coSpikePhase', 'coSpikeLFP', 'LFPall', 'rbAll', '-v7.3')
            
        
            
        end
    end
end

%%




pdfExport = '/space/seh10/6/halgdev/projects/iverzh/ripples/RutishauserLab/figures';
% tractLenStruct = load(fullfile(exportDir, '../../..', 'averageConnectivity_tractLengths.mat'));
fs = 1e3; %1e3, 400;
ctxParc = {'ACC', 'SMA', 'OFC' ,'HIP', 'AMY'};
broadman = {'p24', '8BM', 's32', 'H', 'TGd'};
for v = 6
    for ang = 2 %2:5 %2:6
        tic
        
        
        
        
               
        phLagsAll = nan(1, 1e6); cLags = 1;
        cAP =1;
        coRipPad = 0.25 * fs;
        durTrough = 0; durPeak = 0;
        unitCCGstep = fs / 1e3;
        histEdges = -(coRipPad+(unitCCGstep/2)):unitCCGstep:(coRipPad+(unitCCGstep/2));
        coFd = 20;
        minDur = 50;
        
        
        lfpA =  nan(1e6, length(-100:100)); lfpB =  nan(1e6, length(-100:100)); 
        rbA =  nan(1e6, length(-100:100)); rbB = nan(1e6, length(-100:100));
        compStringAllSubj = []; phLag =  nan(1,1e6);
        nCoSpikeCoAllSubj = []; nCoSpikeNoAllSubj =  [];
        distanceAllSubj =[]; 
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
        c = []; countLFP = 0;countPhase = 0;countCCG = 1;
        for subj = 1:length(subj_list_full) %[1:17 20:length(subj_list_full)] % [1:15, 18:length(subj_list_full)]
            subject = subj_list_full{subj};
            nCo = []; nNo = []; durationAll = [];
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
                filename = sprintf('%s_coRipple_coFire_coRipPad-%01i_IISremv_SU_u%i_%s.mat', subject,fs/4, uA, state);
                load(fullfile(exportDir, filename)); 
                for uBi = 1:length(SUii)
                    uB = SUii(uBi);
                    
                    typeA = units(uA,3);
                    typeB = units(uB,3);



                    coRip = coRipSpike{uA,uB};
                    noRip = noRipSpike{uA,uB};
                    dur = (sum(coRipDur{uA,uB})); %+ 2*0*length(coRipDur{uA,uB})) /fs;
                    durNo = noRipDur{uA,uB} / fs;

                    if dur > 0 && uCh{uA} ~= uCh{uB} %&& strcmp(interactionType{uA, uB}, interact) % || strcmp(interactionType{uA, uB}, 'int-pyr'))
                        sp = coRip(:,coRip(1,:) > 0 & coRip(2,:) > 0);
                        nCoSpikeCo = size(sp,2); 

                        sp = noRip(:,noRip(1,:) > 0 & noRip(2,:) > 0);
                        nCoSpikeNo = size(sp,2);

                        nCo = [nCo, nCoSpikeCo];
                        nNo = [nNo, nCoSpikeNo];
                        durationAll = [durationAll, dur];
                        
                        if (uCh{uA} <= 96 && uCh{uB} > 96) || (uCh{uA} > 96 && uCh{uB} <= 96)
                        d = 1; else
                        d = 0;
                        end


        %                 phLag =


                        if d >= 0 %&& ~contains(compString(end), {'AMY', 'HIP'}) 


%                             zeroLag = (phLags >= - pi/ang & phLags <= pi/ang) & coRipDur{uA,uB}' >= minDur;
%                             nzeroLag = (phLags < - (ang-1)*pi/ang | phLags > (ang-1)*pi/ang)  & coRipDur{uA,uB}' >= minDur;
                            allSpike = coRipDur{uA,uB}' >= minDur;

%                             allLag = phLags((zeroLag | nzeroLag) & (coRip(1,:) > 0 & coRip(2,:) > 0));
%                             allLag = [phLags(coRip(1,:) > 0 & coRip(2,:) > 0)];
%                             phLagsSubj = [phLagsSubj allLag];
%                             phLagsAll(cLags:cLags+length(allLag)-1) = allLag;
%                             cLags = cLags + length(allLag);

                            spMat = coSpikeTimes{uA, uB};
        %                     ;
                            uTa = cell2mat(spMat(2,allSpike)'); uTb = cell2mat(spMat(3,allSpike)');
                            if ~isempty(uTa) && ~isempty(uTb)
                                s = uTa - uTb';
                                cAP = cAP + length(uTa);
                                Nz = histcounts(s, histEdges);
                            else
                                Nz = zeros(1,length(histEdges)-1);
                            end
                            nAspikeAllcoRpeak(countCCG) = length(uTa);


        %                     durPeak = [durPeak + sum(coRipDur{uA,uB}(zeroLag)) / fs];
                            durAll = sum(coRipDur{uA,uB}) / fs;
                            durPeak = sum(coRipDur{uA,uB}(allSpike)) / fs;
    %                         Nz = Nz/durPeak;
                            Nnorm = Nz; %/sum(Nz); durPeak; %sum(Nz);
%                             Nnorm = Nz/durPeak; %sum(Nz);
                            ccgAllcoRpeak(countCCG, :) = Nnorm;
                            durAllcoRpeak(countCCG) = durPeak;

%                             uTa = cell2mat(spMat(2,nzeroLag)); uTb = cell2mat(spMat(3,nzeroLag));
%                             if ~isempty(uTa) && ~isempty(uTb)
%                                 s = uTa - uTb';
%                                 Nnz = histcounts(s, histEdges);
%                             else
                                Nnz = zeros(1,length(histEdges)-1);
%                             end
                            nAspikeAllcoRtrough(countCCG) = 0;
                            nAspikeAllcoR(countCCG) =  nAspikeAllcoRtrough(countCCG) + nAspikeAllcoRpeak(countCCG);


                            uTa = cell2mat(spMat(4,:)'); uTb = cell2mat(spMat(5,:)');
                            if ~isempty(uTa) && ~isempty(uTb)
                                s = uTa - uTb';
                                Ncontrol = histcounts(s, histEdges);
                            else
                                Ncontrol = zeros(1,length(histEdges)-1);
                            end
                            nAspikeAllnoRepoch(countCCG) = length(uTa);



        %                     durTrough = [durTrough + sum(coRipDur{uA,uB}(nzeroLag)) / fs];
%                             durTrough = sum(coRipDur{uA,uB}(nzeroLag)) / fs;
%                             Nnorm = Nnz; %/sum(Nnz); %durTrough; %sum(Nnz);
% %                             Nnorm = Nnz/durTrough; %durTrough; %sum(Nnz);
%                             ccgAllcoRtrough(countCCG, :) = Nnorm;
%                             durAllcoRtrough(countCCG) = durTrough;

                            Nall = Nz + Nnz; %ccgCoR{uA, uB};
        %                     N = N/dur;
        %                     Nnorm = Nall/sum(Nall); %(durTrough+durPeak); %sum(Nall);
%                             Nnorm = Nall/(durTrough+durPeak); %sum(Nall);
        %                     Nnorm = Nall/dur; %sum(Nall);
                            Nnorm = Nall; %/sum(Nall);
                            ccgAllcoR(countCCG, :) = Nnorm;
                            durAllcoR(countCCG) = (durTrough+durPeak);

%                             uTnoRa = ccgNoR{uA,uB,1};
%                             uTnoRb = ccgNoR{uA,uB,2};
%                             s = uTnoRa - uTnoRb';
                            Nab = histcounts(s, histEdges);
                            

                            ccgAllnoRAll(countCCG, :) = Nab;
                            durAllnoRAll(countCCG) = durNo;
%                             nAspikeAllnoRAll(countCCG) = length(uTnoRa);
                            
                            ccgAllnoRepoch(countCCG, :) = Ncontrol;
                            durAllnoRepoch(countCCG) = durAll;
                            ccgCrossArray(countCCG) = d;
                            pairType(countCCG,1) = typeA;
                            pairType(countCCG,2) = typeB;
%                             iiZ = ismember(iCoRip, find(zeroLag));
%                             iiZ = true(1, length(iCoRip));

%                             lfpA(countLFP+1:countLFP + sum(iiZ),:) = lfpMatA(:,iiZ)';
%                             lfpB(countLFP+1:countLFP + sum(iiZ),:) = lfpMatB(:,iiZ)';
%                             rbA(countLFP+1:countLFP + sum(iiZ),:)  = rbMatA(:,iiZ)';
%                             rbB(countLFP+1:countLFP + sum(iiZ),:)  = rbMatB(:,iiZ)';

                     
                            chA = uCh{uA}; chB = uCh{uB};
%                             if ~isempty(LFPall{chA,chB,1})
%                                 iiZ = size(LFPall{chA,chB,1}, 1);
%                                 lfpA(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,1};
%                                 lfpB(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,2};
%                                 rbA(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,1};
%                                 rbB(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,2};
%                                 
%                                 compStringAll_LFP(countLFP+1:countLFP + iiZ) = repmat(compString(end), [1 iiZ]);
% 
% 
%                                 countLFP = countLFP + iiZ;
% 
%                             end
                            
                            if d > 0
                                uTa = cell2mat(spMat(2,allSpike)'); uTb = cell2mat(spMat(3,allSpike)');
                                uPha = cell2mat(spMat(6,allSpike)); uPhb = cell2mat(spMat(7,allSpike));
                                s = (uTa - uTb');
%                                 [b, a] = find(abs(s) <= coFd & abs(s) >= 5);
                                [a, b] = find(s <= coFd & s >= 5);
                                if ~isempty(a)
                                    aii = find(isnan(aSpPhase), 1, 'first');
                                    aSpPhase(aii:aii+(length(a)-1)) = uPha(a);
                                    bii = find(isnan(bSpPhase), 1, 'first');
                                    bSpPhase(bii:bii+(length(b)-1)) = uPhb(b);
                                    
                                    Nphase = Nphase + histcounts(s(a,b), histEdges);

                                end
                                
                                
                            end
                           
                            
                            
%                             compStringAll_phase(countPhase+1:countPhase + length(allLag)) = repmat(compString(end), [1 length(allLag)]);

%                             countLFP = countLFP + sum(iiZ);
%                             countPhase = countPhase + length(allLag);
                            countCCG = countCCG + 1;

    %                         figure; 
    %                         subplot(2,1,1)
    %                         plot(rbMatA(:,iiZ));
    %                         subplot(2,1,2)
    %                         plot(rbMatB(:,iiZ));

    %                         waitforbuttonpress; close;

                        end



                    end
                end
                clc
                fprintf('unit %i done...\n', uA)
             end
             
%              uCh = cell2mat(uCh);
%              for iA = 1:length(uCh)
%                  for iB = chA+1:length(uCh)
%                      chA = uCh(iA); chB = uCh(iB); 
%                     if ~isempty(LFPall{chA,chB,1})
%                         iiZ = size(LFPall{chA,chB,1}, 1);
%                         lfpA(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,1};
%                         lfpB(countLFP+1:countLFP + iiZ,:) = LFPall{chA,chB,2};
%                         rbA(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,1};
%                         rbB(countLFP+1:countLFP + iiZ,:)  = rbAll{chA,chB,2};
% 
%                         countLFP = countLFP + iiZ;
% 
%                     end
%                  end
%              end
             
             

    %          figure;
    % %     plot(-100:100, mean(lfpA([1:1.5e4],:)), 'r'); hold on;
    % %     plot(-100:100, mean(lfpB([1:1.5e4],:)), 'b'); hold on;
    % %     plot(-100:100, mean(lfpA(2.5e4:end,:)), 'r'); hold on;
        % %     plot(-100:100, mean(lfpB(2.5e4:end,:)), 'b'); hold on;
        %         plot(-100:100, mean(lfpA), 'r'); hold on;
        %         plot(-100:100, mean(lfpB), 'b'); hold on;
        %         vline(-99:11:99)
        %         xlim([-50 50])
        %         title(subject)
        % 
        %         waitforbuttonpress; close;

                 c = [c cSubj];

            %     zAll = zscore([nCoSpikeCoAll(distanceAll > 0) nCoSpikeNoAll(distanceAll > 0)]);
            %     nCoSpikeCoAll = zAll(1:sum(distanceAll > 0));
            %     nCoSpikeNoAll = zAll(sum(distanceAll > 0)+1:end);

%                 zAll = zscore([nCoSpikeCoAll(tractDistance >= 0) nCoSpikeNoAll(tractDistance >= 0)]);
%                 nCoSpikeCoAll = zAll(1:sum(tractDistance >= 0));
%                 nCoSpikeNoAll = zAll(sum(tractDistance >= 0)+1:end);
                compStringAllSubj = [compStringAllSubj compString(tractDistance >= 0)];

                nCoSpikeCoAllSubj = [nCoSpikeCoAllSubj nCoSpikeCoAll];
                nCoSpikeNoAllSubj = [nCoSpikeNoAllSubj nCoSpikeNoAll];
                distanceAllSubj = [distanceAllSubj tractDistance(tractDistance >= 0)];
                fprintf('%s   angle %i coFire ver %i\n', subject, ang, v)
                try
            %         figure; polarhistogram(phLagsSubj);%waitforbuttonpress; 
            %         savepdf(gcf, fullfile(pdfExport, sprintf('%s_coRipPL.pdf', subject)));
            %         close all;
                end
                
                
        end
        disp('computing CCGs complete.')
        toc
        
        filename = sprintf('coRipple_coFire_ang%i_v%i_fireProp.mat', ang, v);
%         save(fullfile(exportDir, '..', filename), '-v7.3')
        
    end
end

%%

% close all
clr = brewermap(10,'Paired');
smWin = 7; %*fs/1e3;
smoothMethod = 'gaussian';
tag = 'bmovie_perSpike';
% tag = 'Sternberg_perSpike';
% tag = 'Sternberg';
% tag = 'bothStudies';

xtimes = -fs/4:fs/4;
xl = [-250 250];

ii = 3;
figure('Position',  [846 419 1158 902]);

typeA = {'pyr' , 'pyr', 'int', 'int'};
typeB = {'pyr' , 'int', 'pyr', 'int'};

ccgCrossArray(isnan(ccgCrossArray)) = [];
IIdel = cell2mat(cellfun(@(X) isempty(X), pairType,  'UniformOutput', false));
pairType(IIdel(:,1), :) = [];
for iT = 1:4
    subplot(2,2,iT)
    keep = ~logical(ccgCrossArray);


    %     keep = true(1, size(ccgAllcoR, 1));  
    % keep = keep & sum(ccgAllnoRAll, 2) > 0;
    % keep = sum(ccgAllcoRtrough, 2) > 0;
    %     keep = sum(ccgAllcoR, 2) > 0;
    xT = xtimes >= xl(1) & xtimes <= xl(2);
     
    keep = keep & strcmp(pairType(:,1),typeA{iT})' &  strcmp(pairType(:,2),typeB{iT})';% & sum(ccgAllcoR(1:length(keep),xT), 2) > 0;
    % times 

    D1 = ccgAllcoRpeak(keep, :) ./ durAllcoRpeak(keep)';
    D2 = ccgAllcoRtrough(keep, :) ./ durAllcoRtrough(keep)';
    D3 = ccgAllcoR(keep, :) ./ durAllcoR(keep)';
    D4 = ccgAllnoRepoch(keep, :) ./ durAllnoRepoch(keep)';
    % D4 = ccgAllnoRAll(keep, :) ./ durAllnoRAll(keep)';

    D1 = ccgAllcoRpeak(keep, :) ./ nAspikeAllcoRpeak(keep)' * 1e3;
    D2 = ccgAllcoRtrough(keep, :) ./ nAspikeAllcoRtrough(keep)' * 1e3;
    D3 = ccgAllcoR(keep, :) ./ nAspikeAllcoR(keep)' * 1e3;
    D4 = ccgAllnoRepoch(keep, :) ./ nAspikeAllnoRepoch(keep)' * 1e3;   
    
%     D1 = ccgAllcoRpeak(keep, :);
%     D2 = ccgAllcoRtrough(keep, :) ;
%     D3 = ccgAllcoR(keep, :) ;
%     D4 = ccgAllnoRepoch(keep, :);
    %         D4 = ccgAllnoRAll(keep, :) ./ nAspikeAllnoRAll(keep)';

    %         D1 = ccgAllcoRpeak(keep, :);
    %         D2 = ccgAllcoRtrough(keep, :);
    %         D3 = ccgAllcoR(keep, :);
    %         D4 = ccgAllnoRepoch(keep, :);


    titles = {'zero lag', 'phase offset', 'combined'};
    % 
    % 
    eval(sprintf('D = D%i;', ii))
    D = smoothdata(D, 2, smoothMethod, smWin);
    Dnull = smoothdata(D4, 2, smoothMethod, smWin);
    %             D = D./Dnull;
    %             D(~isfinite(D)) = nan;

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

    %             xlim([-50 50]*fs/1e3)
    %             yl = [min(mean(D(:,xT), 'omitnan'))-0.01 max(mean(D(:,xT), 'omitnan'))+0.01];
    xlim((xl*fs)/1e3)
    %             ylim(yl)
    ylabel('firing rate')
    xlabel('time [ms]')
    %             % grid on;
    fig = gcf;
    fig.Color = 'w';
    % 
    % %             title(sprintf('%s <-> %s %s', ctxParc{iRa}, ctxParc{iRb}, titles{ii}))

%     title(sprintf('%i', sum(ccgAllcoR(keep, xT), 'all')))
    title(sprintf('%s <--> %s',typeA{iT}, typeB{iT}))
end





