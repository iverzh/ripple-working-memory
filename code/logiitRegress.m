%% Logistic Regression with Interaction Term for Eyes Open/Closed Analysis
% This script tests whether the performance drop from eyes-open to 
% eyes-closed differs between conditions A and B

clear; clc;

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

%% Step 2: Create the data matrix for logistic regression
% We need to expand the data into individual observations

% Calculate number of incorrect responses
n_A_open_incorrect = n_A_open_total - n_A_open_correct;
n_A_closed_incorrect = n_A_closed_total - n_A_closed_correct;
n_B_open_incorrect = n_B_open_total - n_B_open_correct;
n_B_closed_incorrect = n_B_closed_total - n_B_closed_correct;

% Create response variable (1 = correct, 0 = incorrect)
y = [ones(n_A_open_correct, 1); zeros(n_A_open_incorrect, 1);     % A, eyes open
     ones(n_A_closed_correct, 1); zeros(n_A_closed_incorrect, 1);  % A, eyes closed
     ones(n_B_open_correct, 1); zeros(n_B_open_incorrect, 1);     % B, eyes open
     ones(n_B_closed_correct, 1); zeros(n_B_closed_incorrect, 1)]; % B, eyes closed

% Create predictor variables
% Condition: 0 = A, 1 = B
% Eyes: 0 = closed, 1 = open

% Condition A, eyes open
condition_A_open = zeros(n_A_open_total, 1);
eyes_A_open = ones(n_A_open_total, 1);

% Condition A, eyes closed
condition_A_closed = zeros(n_A_closed_total, 1);
eyes_A_closed = zeros(n_A_closed_total, 1);

% Condition B, eyes open
condition_B_open = ones(n_B_open_total, 1);
eyes_B_open = ones(n_B_open_total, 1);

% Condition B, eyes closed
condition_B_closed = ones(n_B_closed_total, 1);
eyes_B_closed = zeros(n_B_closed_total, 1);

% Combine all predictors
condition = [condition_A_open; condition_A_closed; 
             condition_B_open; condition_B_closed];
eyes = [eyes_A_open; eyes_A_closed; 
        eyes_B_open; eyes_B_closed];

% Create interaction term
interaction = condition .* eyes;

% Combine predictors into design matrix (including intercept)
X = [ones(size(y)), condition, eyes, interaction];

%% Step 3: Fit the logistic regression model
% Using glmfit for logistic regression
[b, dev, stats] = glmfit([condition, eyes, interaction], y, ...
                         'binomial', 'link', 'logit');

%% Step 4: Display results
fprintf('\n========================================\n');
fprintf('LOGISTIC REGRESSION RESULTS\n');
fprintf('========================================\n\n');

% Parameter estimates
fprintf('Parameter Estimates:\n');
fprintf('-------------------\n');
fprintf('Intercept (??):                     %8.4f (SE = %.4f)\n', b(1), stats.se(1));
fprintf('Condition (??):                     %8.4f (SE = %.4f)\n', b(2), stats.se(2));
fprintf('Eyes (??):                          %8.4f (SE = %.4f)\n', b(3), stats.se(3));
fprintf('Condition ? Eyes Interaction (??):  %8.4f (SE = %.4f)\n\n', b(4), stats.se(4));

% Calculate z-statistics and p-values
z_scores = b ./ stats.se;
p_values = 2 * (1 - normcdf(abs(z_scores)));

fprintf('Statistical Tests:\n');
fprintf('-----------------\n');
fprintf('Parameter                           z-score    p-value\n');
fprintf('Intercept:                          %7.3f    %.4f\n', z_scores(1), p_values(1));
fprintf('Condition:                          %7.3f    %.4f\n', z_scores(2), p_values(2));
fprintf('Eyes:                               %7.3f    %.4f\n', z_scores(3), p_values(3));
fprintf('Condition ? Eyes Interaction:       %7.3f    %.4f\n\n', z_scores(4), p_values(4));

%% Step 5: Interpret the interaction
fprintf('========================================\n');
fprintf('INTERPRETATION OF INTERACTION\n');
fprintf('========================================\n\n');

if p_values(4) < 0.05
    fprintf('The interaction is SIGNIFICANT (p = %.4f)\n', p_values(4));
    if b(4) > 0
        fprintf('Positive interaction coefficient (?? = %.4f) suggests:\n', b(4));
        fprintf('The performance drop from eyes-open to eyes-closed is SMALLER\n');
        fprintf('in Condition B compared to Condition A.\n\n');
    else
        fprintf('Negative interaction coefficient (?? = %.4f) suggests:\n', b(4));
        fprintf('The performance drop from eyes-open to eyes-closed is LARGER\n');
        fprintf('in Condition B compared to Condition A.\n\n');
    end
else
    fprintf('The interaction is NOT SIGNIFICANT (p = %.4f)\n', p_values(4));
    fprintf('There is no statistical evidence that the performance drop\n');
    fprintf('differs between conditions A and B.\n\n');
end

%% Step 6: Calculate and display odds ratios
fprintf('========================================\n');
fprintf('ODDS RATIOS\n');
fprintf('========================================\n\n');

OR = exp(b);
OR_CI_lower = exp(b - 1.96 * stats.se);
OR_CI_upper = exp(b + 1.96 * stats.se);

fprintf('Odds Ratios (95%% CI):\n');
fprintf('--------------------\n');
fprintf('Condition effect:              %.3f (%.3f - %.3f)\n', ...
        OR(2), OR_CI_lower(2), OR_CI_upper(2));
fprintf('Eyes effect:                   %.3f (%.3f - %.3f)\n', ...
        OR(3), OR_CI_lower(3), OR_CI_upper(3));
fprintf('Interaction effect:             %.3f (%.3f - %.3f)\n\n', ...
        OR(4), OR_CI_lower(4), OR_CI_upper(4));

%% Step 7: Calculate predicted probabilities for each group
fprintf('========================================\n');
fprintf('PREDICTED PROBABILITIES\n');
fprintf('========================================\n\n');

% Calculate predicted log-odds for each combination
logodds_A_open = b(1) + b(2)*0 + b(3)*1 + b(4)*0*1;
logodds_A_closed = b(1) + b(2)*0 + b(3)*0 + b(4)*0*0;
logodds_B_open = b(1) + b(2)*1 + b(3)*1 + b(4)*1*1;
logodds_B_closed = b(1) + b(2)*1 + b(3)*0 + b(4)*1*0;

% Convert to probabilities
prob_A_open = exp(logodds_A_open) / (1 + exp(logodds_A_open));
prob_A_closed = exp(logodds_A_closed) / (1 + exp(logodds_A_closed));
prob_B_open = exp(logodds_B_open) / (1 + exp(logodds_B_open));
prob_B_closed = exp(logodds_B_closed) / (1 + exp(logodds_B_closed));

fprintf('Model-predicted probabilities of success:\n');
fprintf('-----------------------------------------\n');
fprintf('Condition A - Eyes Open:    %.3f\n', prob_A_open);
fprintf('Condition A - Eyes Closed:  %.3f\n', prob_A_closed);
fprintf('Condition B - Eyes Open:    %.3f\n', prob_B_open);
fprintf('Condition B - Eyes Closed:  %.3f\n\n', prob_B_closed);

fprintf('Performance drops (open - closed):\n');
fprintf('----------------------------------\n');
fprintf('Condition A: %.3f\n', prob_A_open - prob_A_closed);
fprintf('Condition B: %.3f\n', prob_B_open - prob_B_closed);
fprintf('Difference in drops (A - B): %.3f\n\n', ...
        (prob_A_open - prob_A_closed) - (prob_B_open - prob_B_closed));

%% Step 8: Observed proportions for comparison
fprintf('========================================\n');
fprintf('OBSERVED PROPORTIONS\n');
fprintf('========================================\n\n');

obs_prob_A_open = n_A_open_correct / n_A_open_total;
obs_prob_A_closed = n_A_closed_correct / n_A_closed_total;
obs_prob_B_open = n_B_open_correct / n_B_open_total;
obs_prob_B_closed = n_B_closed_correct / n_B_closed_total;

fprintf('Observed proportions of success:\n');
fprintf('--------------------------------\n');
fprintf('Condition A - Eyes Open:    %.3f\n', obs_prob_A_open);
fprintf('Condition A - Eyes Closed:  %.3f\n', obs_prob_A_closed);
fprintf('Condition B - Eyes Open:    %.3f\n', obs_prob_B_open);
fprintf('Condition B - Eyes Closed:  %.3f\n\n', obs_prob_B_closed);

fprintf('Observed performance drops (open - closed):\n');
fprintf('-------------------------------------------\n');
fprintf('Condition A: %.3f\n', obs_prob_A_open - obs_prob_A_closed);
fprintf('Condition B: %.3f\n', obs_prob_B_open - obs_prob_B_closed);
fprintf('Difference in drops (A - B): %.3f\n\n', ...
        (obs_prob_A_open - obs_prob_A_closed) - (obs_prob_B_open - obs_prob_B_closed));

%% Step 9: Create visualization
figure('Position', [100, 100, 800, 600]);

% Plot observed proportions
subplot(1,2,1);
x_pos = [1, 2, 4, 5];
observed = [obs_prob_A_closed, obs_prob_A_open, obs_prob_B_closed, obs_prob_B_open];
bar(x_pos, observed, 'FaceColor', [0.7 0.7 0.7]);
hold on;

% Add error bars (using binomial standard errors)
se_obs = sqrt(observed .* (1 - observed) ./ [n_A_closed_total, n_A_open_total, ...
                                              n_B_closed_total, n_B_open_total]);
errorbar(x_pos, observed, se_obs, 'k.', 'LineWidth', 1.5);

% Formatting
set(gca, 'XTick', x_pos);
set(gca, 'XTickLabel', {'Closed', 'Open', 'Closed', 'Open'});
ylabel('Proportion Correct');
title('Observed Proportions');
% ylim([0 1]);
text(1.5, -0.1, 'Condition A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(4.5, -0.1, 'Condition B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
grid on;

% Plot model predictions
subplot(1,2,2);
predicted = [prob_A_closed, prob_A_open, prob_B_closed, prob_B_open];
bar(x_pos, predicted, 'FaceColor', [0.3 0.5 0.7]);
hold on;

% Connect with lines to show interaction
plot([1, 2], [prob_A_closed, prob_A_open], 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
plot([4, 5], [prob_B_closed, prob_B_open], 'b-o', 'LineWidth', 2, 'MarkerSize', 8);

% Formatting
set(gca, 'XTick', x_pos);
set(gca, 'XTickLabel', {'Closed', 'Open', 'Closed', 'Open'});
ylabel('Proportion Correct');
title('Model Predictions');
% ylim([0 1]);
text(1.5, -0.1, 'Condition A', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(4.5, -0.1, 'Condition B', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
legend({'Predictions', 'Condition A', 'Condition B'}, 'Location', 'best');
grid on;

sgtitle('Logistic Regression Analysis: Eyes ? Condition Interaction');

fprintf('========================================\n');
fprintf('Analysis complete. See figure for visualization.\n');
fprintf('========================================\n');