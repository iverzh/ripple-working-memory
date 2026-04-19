function out = aggregateRegionDataShuff(in)
% Computes the co-ripple replay rate and a within-subject shuffle null
% distribution across nIter iterations, for a given region pairing.
% Input/output follow the same struct convention as aggregateRegionData.
%
% The shuffle independently permutes the trial indices where channel 1
% fires and where channel 2 fires (within subject groups), then counts
% how often the shuffled indices overlap ? producing a null distribution
% of the co-ripple replay rate.

% Defaults for optional fields
defaults = struct('sameHemisphere', [], ...
                  'loads',          [], ...
                  'loadParam',      3,  ...
                  'nIter',          1e3);
optFields = fieldnames(defaults);
for k = 1:numel(optFields)
    if ~isfield(in, optFields{k})
        in.(optFields{k}) = defaults.(optFields{k});
    end
end

% Unpack inputs
regions        = in.regions;
regionsPlot    = in.regionsPlot;
targetRegionsA = in.targetRegionsA;
targetRegionsB = in.targetRegionsB;
coRcoF         = in.coRcoF;
coRap          = in.coRap;
coRap_mm       = in.coRap_mm;
noRap          = in.noRap;
noRap_mm       = in.noRap_mm;
coRdur         = in.coRdur;
coRdur_mm      = in.coRdur_mm;
noRdur         = in.noRdur;
noRdur_mm      = in.noRdur_mm;
subjID         = in.subjID;
sameHemisphere = in.sameHemisphere;
loads          = in.loads;
loadParam      = in.loadParam;
nIter          = in.nIter;

% Get indices for target regions
iRapl = find(contains(regions, regionsPlot(targetRegionsA)));
iRbpl = find(contains(regions, regionsPlot(targetRegionsB)));

fprintf('aggregateRegionDataShuff: %s  ?  %s  (loadParam = %s, nIter = %d)\n', ...
    strjoin(regionsPlot(targetRegionsA), '/'), ...
    strjoin(regionsPlot(targetRegionsB), '/'), ...
    num2str(loadParam), nIter)

% Compute actual co-ripple replay rate (once, outside the shuffle loop)
coR_replay = zeros(2, 1);
coRap_total = 0;  coRap_mm_total = 0;
noRap_total = 0;  noRap_mm_total = 0;
coRdur_total = 0; coRdur_mm_total = 0;
noRdur_total = 0; noRdur_mm_total = 0;

for Aloop = iRapl
    for Bloop = iRbpl
        [iiA, iiB] = resolveIndices(Aloop, Bloop, iRapl, iRbpl);

        if ~isempty(sameHemisphere)
            if strcmp(regions{iiA}(1), regions{iiB}(1)) ~= sameHemisphere
                continue
            end
        end

        loadMask = ismember(loads(:, iiA, iiB), loadParam);

        ii1 = (coRcoF(:,1,iiA,iiB) > 0  & coRcoF(:,2,iiA,iiB) > 0)  & loadMask;
        ii2 = (coRcoF(:,1,iiA,iiB) == 0 | coRcoF(:,2,iiA,iiB) == 0) & loadMask;
        coR_replay = coR_replay + [sum(ii1); sum(ii2)];

        coRap_total    = coRap_total    + sum(coRap(:,:,iiA,iiB),    'all');
        coRap_mm_total = coRap_mm_total + sum(coRap_mm(:,:,iiA,iiB), 'all');
        noRap_total    = noRap_total    + sum(noRap(:,:,iiA,iiB),    'all');
        noRap_mm_total = noRap_mm_total + sum(noRap_mm(:,:,iiA,iiB), 'all');

        coRdur_total    = coRdur_total    + sum(coRdur(:,:,iiA,iiB),    'all');
        coRdur_mm_total = coRdur_mm_total + sum(coRdur_mm(:,:,iiA,iiB), 'all');
        noRdur_total    = noRdur_total    + sum(noRdur(:,:,iiA,iiB),    'all');
        noRdur_mm_total = noRdur_mm_total + sum(noRdur_mm(:,:,iiA,iiB), 'all');
    end
end

% Shuffle null distribution
shuff = nan(1, nIter);

for iter = 1:nIter
    shuff_replay = zeros(2, 1);

    for Aloop = iRapl
        for Bloop = iRbpl
            [iiA, iiB] = resolveIndices(Aloop, Bloop, iRapl, iRbpl);

            if ~isempty(sameHemisphere)
                if strcmp(regions{iiA}(1), regions{iiB}(1)) ~= sameHemisphere
                    continue
                end
            end

            loadMask = ismember(loads(:, iiA, iiB), loadParam);

            % Shuffle channel-1 firing indices within subject groups
            ind1  = find(coRcoF(:,1,iiA,iiB) > 0 & loadMask);
            group = subjID(coRcoF(:,1,iiA,iiB) > 0 & loadMask, iiA, iiB);
            nPairTrialGroup = arrayfun(@(X) sum(subjID(loadMask, iiA, iiB) == X), unique(group));
            ind1  = shuffle_within_group(ind1, group, nPairTrialGroup);

            % Shuffle channel-2 firing indices within subject groups
            ind2  = find(coRcoF(:,2,iiA,iiB) > 0 & loadMask);
            group = subjID(coRcoF(:,2,iiA,iiB) > 0 & loadMask, iiA, iiB);
            nPairTrialGroup = arrayfun(@(X) sum(subjID(loadMask, iiA, iiB) == X), unique(group));
            ind2  = shuffle_within_group(ind2, group, nPairTrialGroup);

            % Count shuffled overlap and total trials per subject group
            nOverlap = sum(ismember(ind1, ind2));
            nPairTrialGroup = arrayfun(@(X) sum(subjID(loadMask, iiA, iiB) == X), ...
                                       unique(subjID(loadMask, iiA, iiB)));

            shuff_replay = shuff_replay + [nOverlap; sum(nPairTrialGroup) - nOverlap];
        end
    end

    shuff(iter) = shuff_replay(1) / sum(shuff_replay) * 100;

    if mod(iter, 100) == 0
        fprintf('.\n')
    elseif mod(iter, 10) == 0
        fprintf('.')
    end
end
fprintf('\n')

% Pack outputs
out.coR_replay     = coR_replay;                                  % [2x1] actual counts (cofire; not-cofire)
out.coR_rate       = coR_replay(1) / sum(coR_replay) * 100;      % actual replay rate (%)
out.shuff          = shuff;                                        % [1 x nIter] shuffle null distribution (%)
out.coRap_total    = coRap_total;
out.coRap_mm_total = coRap_mm_total;
out.noRap_total    = noRap_total;
out.noRap_mm_total = noRap_mm_total;
out.coRdur_total   = coRdur_total;
out.coRdur_mm_total = coRdur_mm_total;
out.noRdur_total   = noRdur_total;
out.noRdur_mm_total = noRdur_mm_total;

end

% -------------------------------------------------------------------------
function [iiA, iiB] = resolveIndices(Aloop, Bloop, iRapl, iRbpl)
% Return the canonical (min, max) index pair, or (A, B) when sets overlap.
    if all(ismember(iRapl, iRbpl))
        iiA = Aloop;
        iiB = Bloop;
    else
        iiA = min([Aloop Bloop]);
        iiB = max([Aloop Bloop]);
    end
end
