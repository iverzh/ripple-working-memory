function out = aggregateRegionDataSubj(in)
% Per-subject replay rates for co-ripple, no-ripple, fast RT, and slow RT.
% Input/output use the same struct convention as aggregateRegionData.

% Defaults for optional fields
defaults = struct('respTimes', [], ...
                  'loads',     [], ...
                  'loadParam', 3,  ...
                  'nSubj',     []);
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
noRcoF         = in.noRcoF;
subjID         = in.subjID;
respTimes      = in.respTimes;
loads          = in.loads;
loadParam      = in.loadParam;

if isempty(in.nSubj)
    nSubj = max(subjID(:));
else
    nSubj = in.nSubj;
end

splitByRT = ~isempty(respTimes) && ~isempty(loads);
if splitByRT
    rtMedian = median(respTimes(ismember(loads, loadParam)), 'all', 'omitnan');
end

% Get indices for target regions
iRapl = find(contains(regions, regionsPlot(targetRegionsA)));
iRbpl = find(contains(regions, regionsPlot(targetRegionsB)));

% Per-subject accumulators [2 x nSubj]: row 1 = co-fire, row 2 = not co-fire
coR_replay = zeros(2, nSubj);
noR_replay = zeros(2, nSubj);
coR_fast   = zeros(2, nSubj);
coR_slow   = zeros(2, nSubj);

for iS = 1:nSubj
    for Aloop = iRapl
        for Bloop = iRbpl
            if all(ismember(iRapl, iRbpl))
                iiA = Aloop;
                iiB = Bloop;
            else
                iiA = min([Aloop Bloop]);
                iiB = max([Aloop Bloop]);
            end

            % Skip same-electrode pairs
            if strcmp(regions{iiA}, regions{iiB}); continue; end

            subjMask = subjID(:, iiA, iiB) == iS;
            loadMask = ismember(loads(:, iiA, iiB), loadParam);
            baseMask = subjMask & loadMask;

            % Co-ripple replay (all trials for this subject)
            ii1 = (coRcoF(:,1,iiA,iiB) > 0  & coRcoF(:,2,iiA,iiB) > 0)  & baseMask;
            ii2 = (coRcoF(:,1,iiA,iiB) == 0 | coRcoF(:,2,iiA,iiB) == 0) & baseMask;
            coR_replay(:, iS) = coR_replay(:, iS) + [sum(ii1,'omitnan'); sum(ii2,'omitnan')];

            % No-ripple control
            ii1 = (noRcoF(:,1,iiA,iiB) > 0  & noRcoF(:,2,iiA,iiB) > 0)  & baseMask;
            ii2 = (noRcoF(:,1,iiA,iiB) == 0 | noRcoF(:,2,iiA,iiB) == 0) & baseMask;
            noR_replay(:, iS) = noR_replay(:, iS) + [sum(ii1,'omitnan'); sum(ii2,'omitnan')];

            % Fast / slow RT split
            if splitByRT
                fastMask = respTimes(:, iiA, iiB) <= rtMedian & baseMask;
                slowMask = respTimes(:, iiA, iiB) >  rtMedian & baseMask;

                ii1 = (coRcoF(:,1,iiA,iiB) > 0  & coRcoF(:,2,iiA,iiB) > 0)  & fastMask;
                ii2 = (coRcoF(:,1,iiA,iiB) == 0 | coRcoF(:,2,iiA,iiB) == 0) & fastMask;
                coR_fast(:, iS) = coR_fast(:, iS) + [sum(ii1,'omitnan'); sum(ii2,'omitnan')];

                ii1 = (coRcoF(:,1,iiA,iiB) > 0  & coRcoF(:,2,iiA,iiB) > 0)  & slowMask;
                ii2 = (coRcoF(:,1,iiA,iiB) == 0 | coRcoF(:,2,iiA,iiB) == 0) & slowMask;
                coR_slow(:, iS) = coR_slow(:, iS) + [sum(ii1,'omitnan'); sum(ii2,'omitnan')];
            end
        end
    end
end

% Compute per-subject replay rates (fraction of trials with co-fire)
coSubj     = coR_replay(1,:) ./ sum(coR_replay, 1) * 100;
noSubj     = noR_replay(1,:) ./ sum(noR_replay, 1) * 100;
coFastSubj = coR_fast(1,:)   ./ sum(coR_fast,   1) * 100;
coSlowSubj = coR_slow(1,:)   ./ sum(coR_slow,   1) * 100;

% Pack outputs
out.coR_replay  = coR_replay;   % raw counts [2 x nSubj]: co-ripple match
out.noR_replay  = noR_replay;   % raw counts [2 x nSubj]: no-ripple control
out.coR_fast    = coR_fast;     % raw counts [2 x nSubj]: fast-RT trials
out.coR_slow    = coR_slow;     % raw counts [2 x nSubj]: slow-RT trials
out.coSubj      = coSubj;       % co-ripple replay rate per subject (%)
out.noSubj      = noSubj;       % no-ripple replay rate per subject (%)
out.coFastSubj  = coFastSubj;   % fast-RT replay rate per subject (%)
out.coSlowSubj  = coSlowSubj;   % slow-RT replay rate per subject (%)

end
