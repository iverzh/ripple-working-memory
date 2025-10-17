

function [bl, bf] = errobarPatch(x, y, c, w, varargin)

 % Parse optional parameters
p = inputParser;
addParameter(p, 'DevType', 'MAD', @isstr);
parse(p, varargin{:});
DevType = p.Results.DevType;
    

switch DevType
    case 'SEM'
        YY = mean(y, 'omitnan');
        YYdev(1) = std(y, 'omitnan')./sqrt(sum(~isnan(y), 1));
        YYdev(2) = std(y, 'omitnan')./sqrt(sum(~isnan(y), 1));
    case 'MAD'
        YY = median(y, 'omitnan');

        % Compute MAD
        mad_Y = median(abs(y - median(y, 'omitnan')), 'omitnan');

        % Compute robust standard error
        YYdev(1) = mad_Y / sqrt(sum(~isnan(y), 1));
        YYdev(2) = mad_Y / sqrt(sum(~isnan(y), 1));
    case 'IQR'
        YY = median(y, 'omitnan');
        YYdev = abs(quantile(y, [0.25 0.75]) - median(y));
end
if size(x,2) == 1
    for ii = 1:length(x)

        bf = patch([x(ii)-w x(ii)+w x(ii)+w x(ii)-w], [YY(ii)-YYdev(1) YY(ii)-YYdev(1) YY(ii)+YYdev(2) YY(ii)+YYdev(2)], c); hold on;
        bf.FaceAlpha = 0.4;
        bf.EdgeAlpha = 0;

        bl = plot([x(ii)-w x(ii)+w], [YY(ii) YY(ii)], '-'); hold on;
        bl.LineWidth = 2;
        bl.Color = c;

    end
    
else
    [bl, bf] = boundedline(x, YY,  YYdev, 'rs-', 'nan', 'gap'); hold on;
    bl.Color = c;
    bl.LineWidth = 1.5;
    bl.MarkerSize = 4;
    bl.MarkerFaceColor = bl.Color;

    bf.FaceColor = c;
    bf.FaceAlpha = 0.2;
    
end





