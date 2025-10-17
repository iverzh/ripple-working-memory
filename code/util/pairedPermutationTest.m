function [pval, perm_stats, obs_stat] = pairedPermutationTest(data1, data2, numPermutations)
% pairedPermutationTest performs a paired permutation test on two paired samples.
%
%   p = pairedPermutationTest(data1, data2, numPermutations) computes the
%   p-value for the null hypothesis that the mean difference between paired
%   observations in data1 and data2 is zero.
%
%   INPUTS:
%     data1, data2        - Vectors of paired observations (must be of equal length)
%     numPermutations     - Number of permutations to run (default is 10000)
%
%   OUTPUT:
%     p                   - Two-tailed p-value from the permutation test

    % Set default number of permutations if not provided
    if nargin < 3
        numPermutations = 10000;
    end

    % Ensure the vectors are of equal length
    if length(data1) ~= length(data2)
        error('data1 and data2 must have the same number of elements.');
    end

    % Calculate the observed paired differences and the test statistic
    d = data1 - data2;
    obs_stat = median(d);

    % Initialize an array to store permutation statistics
    perm_stats = zeros(numPermutations, 1);
    n = length(d);

    % Loop over the number of permutations
    for i = 1:numPermutations
        % Generate random sign flips (+1 or -1) for each paired difference
        signs = 2 * (rand(n, 1) > 0.5) - 1;
        % Apply sign flips to the differences
        perm_d = d .* signs;
        % Compute the mean difference for this permutation
        perm_stats(i) = median(perm_d);
    end

    % Calculate the two-tailed p-value: proportion of permuted test statistics
    % that are as or more extreme than the observed statistic.
    pval = sum(abs(perm_stats) >= abs(obs_stat)) / numPermutations;
end
