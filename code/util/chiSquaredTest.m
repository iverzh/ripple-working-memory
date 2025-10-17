







function [chi2stat, pValue, expected] = chiSquaredTest(matrix)
    total = sum(matrix(:));
    rowSums = sum(matrix,2);
    colSums = sum(matrix,1);
    expected = (rowSums * colSums) / total; % Compute expected frequencies
    
    % Compute Chi-Squared test
    chi2stat = sum((matrix(:) - expected(:)).^2 ./ expected(:));
    pValue = 1 - chi2cdf(chi2stat, 1); % df = (2-1)*(2-1) = 1
end


