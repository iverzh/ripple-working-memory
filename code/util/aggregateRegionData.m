function out = aggregateRegionData(in)

% Apply defaults for optional fields
defaults = struct('sameHemisphere', [], ...
                  'respTimes',      [], ...
                  'loads',          [], ...
                  'loadParam',      3);
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
coRcoF_mm      = in.coRcoF_mm;
noRcoF         = in.noRcoF;
noRcoF_mm      = in.noRcoF_mm;
coRap          = in.coRap;
coRap_mm       = in.coRap_mm;
noRap          = in.noRap;
noRap_mm       = in.noRap_mm;
coRdur         = in.coRdur;
coRdur_mm      = in.coRdur_mm;
noRdur         = in.noRdur;
noRdur_mm      = in.noRdur_mm;
coR_replay     = in.coR_replay;
sameHemisphere = in.sameHemisphere;
respTimes      = in.respTimes;
loads          = in.loads;
loadParam      = in.loadParam;

% Get indices for target regions
iRapl = find(contains(regions, regionsPlot(targetRegionsA)));
iRbpl = find(contains(regions, regionsPlot(targetRegionsB)));


% Initialize accumulators
replay_match = zeros(2, 1);
replay_mismatch = zeros(2, 1);
replay_norip_match = zeros(2, 1);
replay_norip_mismatch = zeros(2, 1);

% Fast/slow RT split accumulators (matched trials only)
splitByRT = ~isempty(respTimes) && ~isempty(loads) && ~isempty(loadParam);
replay_fast_match = zeros(2, 1);
replay_slow_match = zeros(2, 1);
if splitByRT
    rtMedian = median(respTimes(ismember(loads, loadParam)), 'all', 'omitnan');
end

coRap_total = 0;
coRap_mm_total = 0;
noRap_total = 0;
noRap_mm_total = 0;

coRdur_total = 0;
coRdur_mm_total = 0;
noRdur_total = 0;
noRdur_mm_total = 0;

% Loop through region pairs
for Aloop = iRapl
    for Bloop = iRbpl
        if all(ismember(iRapl, iRbpl))
            iiA = Aloop;
            iiB = Bloop;
        else
            iiA = min([Aloop Bloop]);
            iiB = max([Aloop Bloop]);
        end


        if sum(coR_replay(:, iiA, iiB)) == 0; continue; end

        % Check hemisphere constraint if specified
        if ~isempty(sameHemisphere)
            if strcmp(regions{iiA}(1), regions{iiB}(1)) ~= sameHemisphere
                continue
            end
        end

        % Aggregate matched trials
        dat = [sum(coRcoF(:, 1, iiA, iiB) > 0 & coRcoF(:, 2, iiA, iiB) > 0); ...
               sum(coRcoF(:, 1, iiA, iiB) == 0 | coRcoF(:, 2, iiA, iiB) == 0)];
        replay_match = replay_match + dat;

        % Fast/slow RT split for matched trials
        if splitByRT
            loadMask = ismember(loads(:, iiA, iiB), loadParam);
            fastMask = respTimes(:, iiA, iiB) <= rtMedian & loadMask;
            slowMask = respTimes(:, iiA, iiB) >  rtMedian & loadMask;

            ii1 = (coRcoF(:, 1, iiA, iiB) > 0  & coRcoF(:, 2, iiA, iiB) > 0)  & fastMask;
            ii2 = (coRcoF(:, 1, iiA, iiB) == 0 | coRcoF(:, 2, iiA, iiB) == 0) & fastMask;
            replay_fast_match = replay_fast_match + [sum(ii1); sum(ii2)];

            ii1 = (coRcoF(:, 1, iiA, iiB) > 0  & coRcoF(:, 2, iiA, iiB) > 0)  & slowMask;
            ii2 = (coRcoF(:, 1, iiA, iiB) == 0 | coRcoF(:, 2, iiA, iiB) == 0) & slowMask;
            replay_slow_match = replay_slow_match + [sum(ii1); sum(ii2)];
        end

        % Aggregate mismatched trials
        dat = [sum(coRcoF_mm(:, 1, iiA, iiB) > 0 & coRcoF_mm(:, 2, iiA, iiB) > 0); ...
               sum(coRcoF_mm(:, 1, iiA, iiB) == 0 | coRcoF_mm(:, 2, iiA, iiB) == 0)];
        replay_mismatch = replay_mismatch + dat;

        % Aggregate no-ripple controls (average across iterations)
        dat = [];
        for iter = 1:25
            dat = [dat, [sum(noRcoF(:, 1, iiA, iiB, iter) > 0 & noRcoF(:, 2, iiA, iiB, iter) > 0); ...
                        sum(noRcoF(:, 1, iiA, iiB, iter) == 0 | noRcoF(:, 2, iiA, iiB, iter) == 0)]];
        end
        replay_norip_match = replay_norip_match + round(mean(dat, 2));

        dat = [];
        for iter = 1:25
            dat = [dat, [sum(noRcoF_mm(:, 1, iiA, iiB, iter) > 0 & noRcoF_mm(:, 2, iiA, iiB, iter) > 0); ...
                        sum(noRcoF_mm(:, 1, iiA, iiB, iter) == 0 | noRcoF_mm(:, 2, iiA, iiB, iter) == 0)]];
        end
        replay_norip_mismatch = replay_norip_mismatch + round(mean(dat, 2));

        % Aggregate co-firing rates
        coRap_total = coRap_total + sum(coRap(:, :, iiA, iiB), 'all');
        coRap_mm_total = coRap_mm_total + sum(coRap_mm(:, :, iiA, iiB), 'all');
        noRap_total = noRap_total + sum(noRap(:, :, iiA, iiB), 'all');
        noRap_mm_total = noRap_mm_total + sum(noRap_mm(:, :, iiA, iiB), 'all');

        coRdur_total = coRdur_total + sum(coRdur(:, :, iiA, iiB), 'all');
        coRdur_mm_total = coRdur_mm_total + sum(coRdur_mm(:, :, iiA, iiB), 'all');
        noRdur_total = noRdur_total + sum(noRdur(:, :, iiA, iiB), 'all');
        noRdur_mm_total = noRdur_mm_total + sum(noRdur_mm(:, :, iiA, iiB), 'all');
    end
end

% Pack outputs
out.replay_match          = replay_match;
out.replay_mismatch       = replay_mismatch;
out.replay_norip_match    = replay_norip_match;
out.replay_norip_mismatch = replay_norip_mismatch;
out.replay_fast_match     = replay_fast_match;
out.replay_slow_match     = replay_slow_match;
out.coRap_total           = coRap_total;
out.coRap_mm_total        = coRap_mm_total;
out.noRap_total           = noRap_total;
out.noRap_mm_total        = noRap_mm_total;
out.coRdur_total          = coRdur_total;
out.coRdur_mm_total       = coRdur_mm_total;
out.noRdur_total          = noRdur_total;
out.noRdur_mm_total       = noRdur_mm_total;

end
