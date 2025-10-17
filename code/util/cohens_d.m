function d = cohens_d(group1, group2, varargin)
% COHENS_D Compute Cohen's d effect size
%
% Syntax:
%   d = cohens_d(group1, group2)
%   d = cohens_d(group1, group2, 'pooled')
%   d = cohens_d(group1, group2, 'unpooled')
%
% Inputs:
%   group1 - Vector of values for group 1
%   group2 - Vector of values for group 2
%   method - (optional) 'pooled' (default) or 'unpooled'
%            'pooled' uses pooled standard deviation
%            'unpooled' uses group 1 standard deviation as denominator
%
% Output:
%   d - Cohen's d effect size
%
% Example:
%   group1 = [1, 2, 3, 4, 5];
%   group2 = [3, 4, 5, 6, 7];
%   d = cohens_d(group1, group2);
%   fprintf('Cohen''s d = %.3f\n', d);

    % Input validation
    if nargin < 2
        error('At least two input groups are required');
    end
    
    % Convert to column vectors
    group1 = group1(:);
    group2 = group2(:);
    
    % Remove NaN values
    group1 = group1(~isnan(group1));
    group2 = group2(~isnan(group2));
    
    % Check for empty groups
    if isempty(group1) || isempty(group2)
        error('Groups cannot be empty after removing NaN values');
    end
    
    % Parse optional arguments
    method = 'pooled'; % default
    if nargin > 2
        method = varargin{1};
    end
    
    % Calculate means
    mean1 = mean(group1);
    mean2 = mean(group2);
    
    % Calculate standard deviations
    std1 = std(group1, 0); % using N-1 denominator
    std2 = std(group2, 0);
    
    % Calculate sample sizes
    n1 = length(group1);
    n2 = length(group2);
    
    % Calculate Cohen's d based on method
    switch lower(method)
        case 'pooled'
            % Pooled standard deviation
            pooled_std = sqrt(((n1-1)*std1^2 + (n2-1)*std2^2) / (n1+n2-2));
            d = (mean1 - mean2) / pooled_std;
            
        case 'unpooled'
            % Use group 1 standard deviation
            d = (mean1 - mean2) / std1;
            
        otherwise
            error('Method must be either ''pooled'' or ''unpooled''');
    end
end



