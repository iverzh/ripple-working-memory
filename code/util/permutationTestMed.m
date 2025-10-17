function [pval, permStats, obsStat] = permutationTestMed(x, y, numPermutations)
% permutationTest performs a permutation test comparing two independent groups.
%
% This function calculates the observed difference in means between groups x and y,
% then creates a distribution of the test statistic by randomly permuting the data.
%
% Inputs:
%   x - vector of data from group 1
%   y - vector of data from group 2
%   numPermutations - number of permutations to perform (e.g., 10000)
%
% Outputs:
%   pval - the p-value from the permutation test (two-tailed)
%   permStats - the distribution of permuted test statistics
%   obsStat - the observed difference in means (x - y)
%
% Example usage:
%   x = [0, 0, 0, 1.2, 2.3, 3.4]; % Group 1 (with many zeros)
%   y = [0, 0, 1.0, 1.8, 2.7, 3.6]; % Group 2
%   [pval, permStats, obsStat] = permutationTest(x, y, 10000);
%   fprintf('Observed difference: %.4f\n', obsStat);
%   fprintf('Permutation p-value: %.4f\n', pval);

    % Compute the observed test statistic: difference in means
    obsStat = median(x) - median(y);
%     obsStat = mean(x) - mean(y);
    
    % Combine the two groups into one dataset
    combinedData = [x(:); y(:)];  % ensure column vector
    n1 = length(x);
    totalN = length(combinedData);
    
    % Initialize array to store the test statistic for each permutation
    permStats = zeros(numPermutations, 1);
    
    % Perform the permutation test
    for i = 1:numPermutations
        % Shuffle the combined data indices
        shuffledIndices = randperm(totalN);
        % Split the data according to the size of group x
        permGroup1 = combinedData(shuffledIndices(1:n1));
        permGroup2 = combinedData(shuffledIndices(n1+1:end));
        % Compute the difference in means for the permuted groups
        permStats(i) = median(permGroup1) - median(permGroup2);
%         permStats(i) = mean(permGroup1) - mean(permGroup2);
    end
    
    % Calculate the two-tailed p-value:
    % Count how many permuted statistics are as or more extreme than the observed statistic.
    pval = mean(abs(permStats) >= abs(obsStat));
end
