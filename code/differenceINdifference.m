%% Difference-in-Differences Test for Eyes Open/Closed Analysis
% This script tests whether the performance drop from eyes-open to 
% eyes-closed differs between conditions A and B using:
% 1. Standard error approach with z-test
% 2. Bootstrap confidence intervals

clc; close all;

%% Step 1: Enter your data
% Replace these with your actual counts
% Format: number of participants who got the task RIGHT in each condition
% Condition A
n_A_open_correct = coR_replay_pl(1);    % Number correct with eyes open in condition A
n_A_open_total = sum(coR_replay_pl);       % Total participants with eyes open in condition A
n_A_closed_correct = coR_replay_mm_pl(1);   % Number correct with eyes closed in condition A
n_A_closed_total = sum(coR_replay_mm_pl);     % Total participants with eyes closed in condition A

% Condition B
n_B_open_correct = noR_replay_pl(1);    % Number correct with eyes open in condition B
n_B_open_total = sum(noR_replay_pl);       % Total participants with eyes open in condition B
n_B_closed_correct = noR_replay_mm_pl(1);    % Number correct with eyes closed in condition B
n_B_closed_total = sum(noR_replay_mm_pl);     % Total participants with eyes closed in condition B

%% Step 2: Calculate observed proportions and differences
% Calculate proportions
p_A_open = n_A_open_correct / n_A_open_total;
p_A_closed = n_A_closed_correct / n_A_closed_total;
p_B_open = n_B_open_correct / n_B_open_total;
p_B_closed = n_B_closed_correct / n_B_closed_total;

% Calculate differences (performance drop when closing eyes)
diff_A = p_A_open - p_A_closed;  % Drop in condition A
diff_B = p_B_open - p_B_closed;  % Drop in condition B

% Calculate difference-in-differences
diff_in_diff = diff_A - diff_B;

%% Step 3: Method 1 - Standard Error Approach with Z-test
fprintf('\n========================================\n');
fprintf('METHOD 1: STANDARD ERROR APPROACH\n');
fprintf('========================================\n\n');

% Calculate standard errors for each proportion
se_A_open = sqrt(p_A_open * (1 - p_A_open) / n_A_open_total);
se_A_closed = sqrt(p_A_closed * (1 - p_A_closed) / n_A_closed_total);
se_B_open = sqrt(p_B_open * (1 - p_B_open) / n_B_open_total);
se_B_closed = sqrt(p_B_closed * (1 - p_B_closed) / n_B_closed_total);

% Calculate standard errors for the differences
% Assuming independence between eyes open and closed conditions
se_diff_A = sqrt(se_A_open^2 + se_A_closed^2);
se_diff_B = sqrt(se_B_open^2 + se_B_closed^2);

% Calculate standard error for difference-in-differences
se_diff_in_diff = sqrt(se_diff_A^2 + se_diff_B^2);

% Calculate z-statistic
z_stat = diff_in_diff / se_diff_in_diff;

% Calculate p-value (two-tailed)
p_value = 2 * (1 - normcdf(abs(z_stat)));

% Calculate 95% confidence interval
ci_lower = diff_in_diff - 1.96 * se_diff_in_diff;
ci_upper = diff_in_diff + 1.96 * se_diff_in_diff;

% Display results
fprintf('Observed Proportions:\n');
fprintf('--------------------\n');
fprintf('Condition A - Eyes Open:    %.3f (SE = %.3f)\n', p_A_open, se_A_open);
fprintf('Condition A - Eyes Closed:  %.3f (SE = %.3f)\n', p_A_closed, se_A_closed);
fprintf('Condition B - Eyes Open:    %.3f (SE = %.3f)\n', p_B_open, se_B_open);
fprintf('Condition B - Eyes Closed:  %.3f (SE = %.3f)\n\n', p_B_closed, se_B_closed);

fprintf('Performance Drops (open - closed):\n');
fprintf('----------------------------------\n');
fprintf('Condition A: %.3f (SE = %.3f)\n', diff_A, se_diff_A);
fprintf('Condition B: %.3f (SE = %.3f)\n\n', diff_B, se_diff_B);

fprintf('Difference-in-Differences Test:\n');
fprintf('-------------------------------\n');
fprintf('Difference-in-differences: %.3f\n', diff_in_diff);
fprintf('Standard error:            %.3f\n', se_diff_in_diff);
fprintf('Z-statistic:               %.3f\n', z_stat);
fprintf('P-value (two-tailed):      %.4f\n', p_value);
fprintf('95%% CI:                    [%.3f, %.3f]\n\n', ci_lower, ci_upper);

% Interpretation
if p_value < 0.05
    fprintf('Result: SIGNIFICANT (p = %.4f)\n', p_value);
    if diff_in_diff > 0
        fprintf('The performance drop is SMALLER in Condition A than in Condition B.\n');
        fprintf('Condition A shows %.1f percentage points less drop than Condition B.\n\n', ...
                abs(diff_in_diff)*100);
    else
        fprintf('The performance drop is LARGER in Condition A than in Condition B.\n');
        fprintf('Condition A shows %.1f percentage points more drop than Condition B.\n\n', ...
                abs(diff_in_diff)*100);
    end
else
    fprintf('Result: NOT SIGNIFICANT (p = %.4f)\n', p_value);
    fprintf('No statistical evidence that performance drops differ between conditions.\n\n');
end

%% Step 4: Method 2 - Bootstrap Confidence Intervals
fprintf('========================================\n');
fprintf('METHOD 2: BOOTSTRAP CONFIDENCE INTERVALS\n');
fprintf('========================================\n\n');

% Set bootstrap parameters
n_bootstrap = 10000;
rng(123); % Set seed for reproducibility

% Initialize storage for bootstrap samples
boot_diff_in_diff = zeros(n_bootstrap, 1);

% Run bootstrap
fprintf('Running %d bootstrap iterations...\n', n_bootstrap);
for i = 1:n_bootstrap
    % Resample for Condition A - Eyes Open
    boot_A_open = binornd(1, p_A_open, n_A_open_total, 1);
    boot_p_A_open = mean(boot_A_open);
    
    % Resample for Condition A - Eyes Closed
    boot_A_closed = binornd(1, p_A_closed, n_A_closed_total, 1);
    boot_p_A_closed = mean(boot_A_closed);
    
    % Resample for Condition B - Eyes Open
    boot_B_open = binornd(1, p_B_open, n_B_open_total, 1);
    boot_p_B_open = mean(boot_B_open);
    
    % Resample for Condition B - Eyes Closed
    boot_B_closed = binornd(1, p_B_closed, n_B_closed_total, 1);
    boot_p_B_closed = mean(boot_B_closed);
    
    % Calculate bootstrap difference-in-differences
    boot_diff_A = boot_p_A_open - boot_p_A_closed;
    boot_diff_B = boot_p_B_open - boot_p_B_closed;
    boot_diff_in_diff(i) = boot_diff_A - boot_diff_B;
end

% Calculate bootstrap statistics
boot_mean = mean(boot_diff_in_diff);
boot_se = std(boot_diff_in_diff);
boot_ci = prctile(boot_diff_in_diff, [2.5, 97.5]);
boot_p_value = 2 * min(sum(boot_diff_in_diff <= 0), sum(boot_diff_in_diff >= 0)) / n_bootstrap;

% Display bootstrap results
fprintf('\nBootstrap Results (%d iterations):\n', n_bootstrap);
fprintf('----------------------------------\n');
fprintf('Mean difference-in-differences:    %.3f\n', boot_mean);
fprintf('Bootstrap standard error:          %.3f\n', boot_se);
fprintf('95%% Bootstrap CI (percentile):     [%.3f, %.3f]\n', boot_ci(1), boot_ci(2));
fprintf('Bootstrap p-value:                 %.4f\n\n', boot_p_value);

% Check if CI includes zero
if boot_ci(1) > 0 || boot_ci(2) < 0
    fprintf('Result: SIGNIFICANT (CI does not include 0)\n');
    if boot_mean > 0
        fprintf('The performance drop is SMALLER in Condition A than in Condition B.\n\n');
    else
        fprintf('The performance drop is LARGER in Condition A than in Condition B.\n\n');
    end
else
    fprintf('Result: NOT SIGNIFICANT (CI includes 0)\n');
    fprintf('No statistical evidence that performance drops differ between conditions.\n\n');
end

%% Step 5: Create visualizations
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Bar plot with differences
subplot(1,3,1);
conditions = categorical({'Condition A', 'Condition B'});
differences = [diff_A, diff_B];
b = bar(conditions, differences);
b.FaceColor = [0.3 0.5 0.7];
hold on;

% Add error bars
errorbar(1:2, differences, [se_diff_A, se_diff_B], 'k.', 'LineWidth', 1.5);

% Add horizontal line at 0
yline(0, 'k--', 'LineWidth', 1);

% Formatting
ylabel('Performance Drop (Open - Closed)');
title('Performance Drops by Condition');
% ylim([min(differences)-0.1, max(differences)+0.1]);
grid on;

% Add text showing difference-in-differences
text(1.5, max(differences)*0.8, sprintf('Diff-in-Diff = %.3f', diff_in_diff), ...
     'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);

% Subplot 2: Interaction plot
subplot(1,3,2);
eyes_states = [0 1]; % 0 = closed, 1 = open
perf_A = [p_A_closed, p_A_open];
perf_B = [p_B_closed, p_B_open];

plot(eyes_states, perf_A, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(eyes_states, perf_B, 'bs-', 'LineWidth', 2, 'MarkerSize', 8);

% Add annotations showing drops
mid_A = (perf_A(1) + perf_A(2)) / 2;
mid_B = (perf_B(1) + perf_B(2)) / 2;
text(0.5, mid_A, sprintf('Drop = %.3f', diff_A), ...
     'HorizontalAlignment', 'center', 'Color', 'r', 'FontWeight', 'bold');
text(0.5, mid_B, sprintf('Drop = %.3f', diff_B), ...
     'HorizontalAlignment', 'center', 'Color', 'b', 'FontWeight', 'bold');

% Formatting
set(gca, 'XTick', [0 1]);
set(gca, 'XTickLabel', {'Eyes Closed', 'Eyes Open'});
ylabel('Proportion Correct');
title('Interaction Plot');
legend({'Condition A', 'Condition B'}, 'Location', 'best');
% ylim([0 1]);
grid on;

% Subplot 3: Bootstrap distribution
subplot(1,3,3);
histogram(boot_diff_in_diff, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
hold on;

% Add vertical lines for observed value and CI
xline(diff_in_diff, 'r-', 'LineWidth', 2);
xline(boot_ci(1), 'b--', 'LineWidth', 1.5);
xline(boot_ci(2), 'b--', 'LineWidth', 1.5);
xline(0, 'k--', 'LineWidth', 1);

% Formatting
xlabel('Difference-in-Differences');
ylabel('Frequency');
title(sprintf('Bootstrap Distribution (n=%d)', n_bootstrap));
legend({'Bootstrap samples', 'Observed', '95% CI', '', 'Null (0)'}, ...
       'Location', 'best');
grid on;

sgtitle('Difference-in-Differences Analysis');

%% Step 6: Effect size calculation
fprintf('========================================\n');
fprintf('EFFECT SIZE MEASURES\n');
fprintf('========================================\n\n');

% Cohen's h for each comparison
h_A = 2 * asin(sqrt(p_A_open)) - 2 * asin(sqrt(p_A_closed));
h_B = 2 * asin(sqrt(p_B_open)) - 2 * asin(sqrt(p_B_closed));
h_diff = h_A - h_B;

fprintf('Cohen''s h (effect size for proportions):\n');
fprintf('----------------------------------------\n');
fprintf('Condition A effect size: %.3f\n', h_A);
fprintf('Condition B effect size: %.3f\n', h_B);
fprintf('Difference in effect sizes: %.3f\n\n', h_diff);

% Interpretation of Cohen's h
if abs(h_diff) < 0.2
    effect_interp = 'small';
elseif abs(h_diff) < 0.5
    effect_interp = 'small to medium';
elseif abs(h_diff) < 0.8
    effect_interp = 'medium to large';
else
    effect_interp = 'large';
end

fprintf('Effect size interpretation: %s\n\n', effect_interp);

%% Step 7: Summary table
fprintf('========================================\n');
fprintf('SUMMARY TABLE\n');
fprintf('========================================\n\n');

% Create summary data
summary_data = {
    'Condition A - Eyes Open', p_A_open, se_A_open;
    'Condition A - Eyes Closed', p_A_closed, se_A_closed;
    'Condition B - Eyes Open', p_B_open, se_B_open;
    'Condition B - Eyes Closed', p_B_closed, se_B_closed;
    '', NaN, NaN;
    'Drop in Condition A', diff_A, se_diff_A;
    'Drop in Condition B', diff_B, se_diff_B;
    'Difference-in-Differences', diff_in_diff, se_diff_in_diff;
};

% Display table
fprintf('%-30s %10s %10s\n', 'Group', 'Proportion', 'SE');
fprintf('%-30s %10s %10s\n', repmat('-', 1, 30), repmat('-', 1, 10), repmat('-', 1, 10));
for i = 1:size(summary_data, 1)
    if ~isnan(summary_data{i,2})
        fprintf('%-30s %10.3f %10.3f\n', summary_data{i,:});
    else
        fprintf('\n');
    end
end

fprintf('\n========================================\n');
fprintf('Analysis complete. See figure for visualizations.\n');
fprintf('========================================\n');