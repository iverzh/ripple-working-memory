function [sttc, p_value] = compute_sttc(spikes1, spikes2, dt, total_time, varargin)
% COMPUTE_STTC - Calculate Spike Time Tiling Coefficient (Cutts & Eglen, 2014)
%
% The STTC is a rate-corrected measure of neural synchrony that ranges from
% -1 to +1 and is independent of firing rate.
%
% INPUTS:
%   spikes1     - Vector of spike times for neuron 1 (in ms or seconds)
%   spikes2     - Vector of spike times for neuron 2 (in ms or seconds)
%   dt          - Coincidence window (same units as spike times)
%   total_time  - Total recording duration (same units as spike times)
%
% OPTIONAL INPUTS (as name-value pairs):
%   'n_surrogates' - Number of surrogates for significance testing (default: 0)
%   'jitter_window' - Window for jitter-based surrogates in ms (default: 20)
%   'method'        - Surrogate method: 'jitter' or 'shift' (default: 'jitter')
%   'verbose'       - Display progress (default: false)
%
% OUTPUTS:
%   sttc        - Spike Time Tiling Coefficient [-1 to 1]
%   p_value     - Significance level from surrogate testing (if requested)
%
% EXAMPLE:
%   spikes1 = [10, 25, 40, 55, 70, 85];
%   spikes2 = [12, 27, 42, 60, 75, 90];
%   sttc = compute_sttc(spikes1, spikes2, 5, 100);
%
% Author: Based on Cutts & Eglen (2014) J Neurosci Methods
% Reference: https://doi.org/10.1016/j.jneumeth.2014.01.002

% Parse optional inputs
p = inputParser;
addParameter(p, 'n_surrogates', 0, @isnumeric);
addParameter(p, 'jitter_window', 20, @isnumeric);
addParameter(p, 'method', 'jitter', @ischar);
addParameter(p, 'verbose', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

% Input validation
if isempty(spikes1) || isempty(spikes2)
    sttc = NaN;
    p_value = NaN;
    if opts.verbose
        warning('Empty spike train(s) provided');
    end
    return;
end

% Ensure column vectors
spikes1 = spikes1(:);
spikes2 = spikes2(:);

% Sort spike times (in case they're not sorted)
% spikes1 = sort(spikes1);
% spikes2 = sort(spikes2);

% Calculate the four key components
% TA: proportion of total time within ?dt of any spike in train A
% TB: proportion of total time within ?dt of any spike in train B  
% PA: proportion of spikes from A that fall within ?dt of any spike in B
% PB: proportion of spikes from B that fall within ?dt of any spike in A

TA = calculate_time_coverage(spikes1, dt, total_time);
TB = calculate_time_coverage(spikes2, dt, total_time);
PA = calculate_spike_proportion(spikes1, spikes2, dt);
PB = calculate_spike_proportion(spikes2, spikes1, dt);

% Calculate STTC using the formula from Cutts & Eglen (2014)
% STTC = 0.5 * [(PA - TB)/(1 - PA*TB) + (PB - TA)/(1 - PB*TA)]

% Handle edge cases to avoid division by zero
if PA * TB == 1
    term1 = 1;  % When PA*TB = 1, the limit approaches 1
else
    term1 = (PA - TB) / (1 - PA * TB);
end

if PB * TA == 1
    term2 = 1;  % When PB*TA = 1, the limit approaches 1
else
    term2 = (PB - TA) / (1 - PB * TA);
end

sttc = 0.5 * (term1 + term2);

% Significance testing with surrogates if requested
if opts.n_surrogates > 0
    if opts.verbose
        fprintf('Computing %d surrogates for significance testing...\n', opts.n_surrogates);
    end
    
    sttc_surrogates = zeros(opts.n_surrogates, 1);
    
    for i = 1:opts.n_surrogates
        if opts.verbose && mod(i, 100) == 0
            fprintf('  Surrogate %d/%d\n', i, opts.n_surrogates);
        end
        
        % Generate surrogate spike trains
        if strcmpi(opts.method, 'jitter')
            surr_spikes1 = jitter_spikes(spikes1, opts.jitter_window, total_time);
            surr_spikes2 = jitter_spikes(spikes2, opts.jitter_window, total_time);
        elseif strcmpi(opts.method, 'shift')
            % Circular shift method
            shift = rand() * total_time;
            surr_spikes1 = mod(spikes1 + shift, total_time);
            surr_spikes2 = spikes2;  % Keep one train fixed
        else
            error('Unknown surrogate method: %s', opts.method);
        end
        
        TA = calculate_time_coverage(surr_spikes1, dt, total_time);
        TB = calculate_time_coverage(surr_spikes2, dt, total_time);
        PA = calculate_spike_proportion(surr_spikes1, surr_spikes2, dt);
        PB = calculate_spike_proportion(surr_spikes2, surr_spikes1, dt);
                   % Handle edge cases to avoid division by zero
        if PA * TB == 1
            term1 = 1;  % When PA*TB = 1, the limit approaches 1
        else
            term1 = (PA - TB) / (1 - PA * TB);
        end

        if PB * TA == 1
            term2 = 1;  % When PB*TA = 1, the limit approaches 1
        else
            term2 = (PB - TA) / (1 - PB * TA);
        end

        
        % Calculate STTC for surrogate
        sttc_surrogates(i) = 0.5 * (term1 + term2);
    end
    
    % Calculate p-value (two-tailed test)
    p_value = sum(abs(sttc_surrogates) >= abs(sttc)) / opts.n_surrogates;
    
    if opts.verbose
        fprintf('STTC = %.4f, p-value = %.4f\n', sttc, p_value);
    end
else
    p_value = NaN;
end

end

%% Helper Functions

function T = calculate_time_coverage(spikes, dt, total_time)
% Calculate the proportion of total time covered by tiles around spikes
%
% This function merges overlapping intervals efficiently

if isempty(spikes)
    T = 0;
    return;
end

% Create intervals [spike - dt, spike + dt]
intervals = [spikes - dt, spikes + dt];

% Clip to recording boundaries
intervals(:,1) = max(0, intervals(:,1));
% intervals(:,2) = min(total_time, intervals(:,2));

% Sort intervals by start time
intervals = sortrows(intervals);

% Merge overlapping intervals
merged = intervals(1,:);
for i = 2:size(intervals, 1)
    if intervals(i,1) <= merged(end,2)
        % Overlapping or touching intervals - merge them
        merged(end,2) = max(merged(end,2), intervals(i,2));
    else
        % Non-overlapping interval - add as new row
        merged = [merged; intervals(i,:)];
    end
end

% Calculate total covered time
covered_time = sum(merged(:,2) - merged(:,1));

T = covered_time / total_time;

end

function P = calculate_spike_proportion(spikes_ref, spikes_target, dt)
% Calculate proportion of spikes_ref that fall within ?dt of any spike in spikes_target

if isempty(spikes_ref) || isempty(spikes_target)
    P = 0;
    return;
end

n_ref = length(spikes_ref);
n_within = 0;


% s = spikes_ref - spikes_target';
% closestSpikeTime = min(abs(s), [], 2);
% n_within = sum(closestSpikeTime <= dt);
% 
% For each spike in reference train
for i = 1:n_ref
    spike_time = spikes_ref(i);
    
    % Check if any target spike is within ?dt
    % Use binary search for efficiency with sorted spike trains
    min_time = spike_time - dt;
    max_time = spike_time + dt;
    
    % Find indices of target spikes in the window
    idx_start = find(spikes_target >= min_time, 1, 'first');
    idx_end = find(spikes_target <= max_time, 1, 'last');
    
    if ~isempty(idx_start) && ~isempty(idx_end) && idx_start <= idx_end
        n_within = n_within + 1;
    end
end

P = n_within / n_ref;

end

function jittered_spikes = jitter_spikes(spikes, jitter_window, total_time)
% Jitter spikes within specified window (interval jitter method)
%
% This implements the conditional inference framework from Amarasingham et al. (2012)

if isempty(spikes)
    jittered_spikes = spikes;
    return;
end

% Define jitter intervals
n_intervals = ceil(total_time / jitter_window);
interval_edges = (0:n_intervals) * jitter_window;
interval_edges(end) = total_time;

% Assign spikes to intervals
[~, interval_idx] = histc(spikes, interval_edges);

% Jitter each spike within its interval
jittered_spikes = zeros(size(spikes));
for i = 1:length(spikes)
    if interval_idx(i) > 0 && interval_idx(i) <= n_intervals
        interval_start = interval_edges(interval_idx(i));
        interval_end = interval_edges(interval_idx(i) + 1);
        
        % Uniform jitter within interval
        jittered_spikes(i) = interval_start + ...
            rand() * (interval_end - interval_start);
    end
end

jittered_spikes = sort(jittered_spikes);

end

%% Additional Analysis Functions

% function [sttc_matrix, pairs] = compute_sttc_population(spike_trains, dt, total_time, varargin)
% % COMPUTE_STTC_POPULATION - Calculate STTC for all pairs in a population
% %
% % INPUTS:
% %   spike_trains - Cell array of spike time vectors
% %   dt          - Coincidence window
% %   total_time  - Total recording duration
% %
% % OUTPUTS:
% %   sttc_matrix - N x N matrix of STTC values
% %   pairs       - Structure with pairwise information
% 
% n_neurons = length(spike_trains);
% sttc_matrix = NaN(n_neurons, n_neurons);
% 
% % Parse options
% p = inputParser;
% addParameter(p, 'parallel', false, @islogical);
% addParameter(p, 'n_surrogates', 0, @isnumeric);
% parse(p, varargin{:});
% opts = p.Results;
% 
% % Initialize pairs structure
% n_pairs = n_neurons * (n_neurons - 1) / 2;
% pairs = struct('neuron1', zeros(n_pairs, 1), ...
%                'neuron2', zeros(n_pairs, 1), ...
%                'sttc', zeros(n_pairs, 1), ...
%                'p_value', zeros(n_pairs, 1));
% 
% pair_idx = 0;
% 
% if opts.parallel && license('test', 'Distrib_Computing_Toolbox')
%     % Parallel computation
%     parfor i = 1:n_neurons
%         for j = i+1:n_neurons
%             [sttc_val, p_val] = compute_sttc(spike_trains{i}, spike_trains{j}, ...
%                                               dt, total_time, ...
%                                               'n_surrogates', opts.n_surrogates);
%             sttc_matrix(i,j) = sttc_val;
%             sttc_matrix(j,i) = sttc_val;
%         end
%     end
% else
%     % Serial computation
%     for i = 1:n_neurons
%         for j = i+1:n_neurons
%             pair_idx = pair_idx + 1;
%             
%             [sttc_val, p_val] = compute_sttc(spike_trains{i}, spike_trains{j}, ...
%                                               dt, total_time, ...
%                                               'n_surrogates', opts.n_surrogates);
%             
%             sttc_matrix(i,j) = sttc_val;
%             sttc_matrix(j,i) = sttc_val;
%             
%             pairs.neuron1(pair_idx) = i;
%             pairs.neuron2(pair_idx) = j;
%             pairs.sttc(pair_idx) = sttc_val;
%             pairs.p_value(pair_idx) = p_val;
%         end
%         
%         % Diagonal is 1 by definition
%         sttc_matrix(i,i) = 1;
%     end
% end
% 
% % Apply FDR correction if p-values were calculated
% if opts.n_surrogates > 0
%     pairs.p_value_fdr = fdr_correction(pairs.p_value);
% end
% 
% end
% 
% function p_adj = fdr_correction(p_values, alpha)
% % Benjamini-Hochberg FDR correction
% if nargin < 2
%     alpha = 0.05;
% end
% 
% n = length(p_values);
% [p_sorted, sort_idx] = sort(p_values);
% 
% % Find largest i such that P(i) <= (i/m) * alpha
% thresh_idx = find(p_sorted <= (1:n)'/n * alpha, 1, 'last');
% 
% if isempty(thresh_idx)
%     p_adj = ones(size(p_values));
% else
%     p_adj = p_values;
%     p_adj(p_values <= p_sorted(thresh_idx)) = p_values(p_values <= p_sorted(thresh_idx));
%     p_adj(p_values > p_sorted(thresh_idx)) = 1;
% end
% 
% end
% 
% %% Visualization Functions
% 
% function plot_sttc_analysis(spikes1, spikes2, dt, total_time, sttc, varargin)
% % PLOT_STTC_ANALYSIS - Visualize spike trains and STTC calculation
% %
% % Creates a comprehensive figure showing:
% %   - Raster plot of both spike trains
% %   - Time tiles around each spike
% %   - STTC value and components
% 
% p = inputParser;
% addParameter(p, 'figure_handle', [], @ishandle);
% parse(p, varargin{:});
% opts = p.Results;
% 
% if isempty(opts.figure_handle)
%     figure('Position', [100, 100, 1200, 600]);
% else
%     figure(opts.figure_handle);
% end
% 
% % Calculate STTC components
% TA = calculate_time_coverage(spikes1, dt, total_time);
% TB = calculate_time_coverage(spikes2, dt, total_time);
% PA = calculate_spike_proportion(spikes1, spikes2, dt);
% PB = calculate_spike_proportion(spikes2, spikes1, dt);
% 
% % Subplot 1: Raster plot with tiles
% subplot(3, 1, 1);
% hold on;
% 
% % Plot tiles for neuron 1
% for i = 1:length(spikes1)
%     rectangle('Position', [spikes1(i)-dt, 0.8, 2*dt, 0.4], ...
%               'FaceColor', [0.8, 0.8, 1], 'EdgeColor', 'none');
% end
% 
% % Plot tiles for neuron 2
% for i = 1:length(spikes2)
%     rectangle('Position', [spikes2(i)-dt, 1.8, 2*dt, 0.4], ...
%               'FaceColor', [1, 0.8, 0.8], 'EdgeColor', 'none');
% end
% 
% % Plot spikes
% plot(spikes1, ones(size(spikes1)), 'b|', 'MarkerSize', 20, 'LineWidth', 2);
% plot(spikes2, 2*ones(size(spikes2)), 'r|', 'MarkerSize', 20, 'LineWidth', 2);
% 
% xlim([0, total_time]);
% ylim([0.5, 2.5]);
% ylabel('Neuron');
% title(sprintf('Spike Trains with Tiles (dt = %.1f)', dt));
% set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
% grid on;
% 
% % Subplot 2: Coincident spikes
% subplot(3, 1, 2);
% hold on;
% 
% % Find coincident spikes
% coincident_mask1 = false(size(spikes1));
% coincident_mask2 = false(size(spikes2));
% 
% for i = 1:length(spikes1)
%     if any(abs(spikes2 - spikes1(i)) <= dt)
%         coincident_mask1(i) = true;
%     end
% end
% 
% for i = 1:length(spikes2)
%     if any(abs(spikes1 - spikes2(i)) <= dt)
%         coincident_mask2(i) = true;
%     end
% end
% 
% % Plot all spikes
% plot(spikes1, ones(size(spikes1)), 'bo', 'MarkerSize', 8);
% plot(spikes2, 2*ones(size(spikes2)), 'ro', 'MarkerSize', 8);
% 
% % Highlight coincident spikes
% plot(spikes1(coincident_mask1), ones(sum(coincident_mask1), 1), 'b*', ...
%      'MarkerSize', 12, 'LineWidth', 2);
% plot(spikes2(coincident_mask2), 2*ones(sum(coincident_mask2), 1), 'r*', ...
%      'MarkerSize', 12, 'LineWidth', 2);
% 
% xlim([0, total_time]);
% ylim([0.5, 2.5]);
% ylabel('Neuron');
% title('Coincident Spikes (starred)');
% set(gca, 'YTick', [1, 2], 'YTickLabel', {'Neuron 1', 'Neuron 2'});
% grid on;
% 
% % Subplot 3: STTC components
% subplot(3, 1, 3);
% 
% % Create bar plot of components
% components = [TA, TB, PA, PB, sttc];
% component_names = {'T_A', 'T_B', 'P_A', 'P_B', 'STTC'};
% colors = [0.5, 0.5, 1;    % Light blue for TA
%           1, 0.5, 0.5;    % Light red for TB
%           0, 0, 1;        % Blue for PA
%           1, 0, 0;        % Red for PB
%           0, 0.8, 0];     % Green for STTC
% 
% b = bar(components);
% b.FaceColor = 'flat';
% for i = 1:length(components)
%     b.CData(i,:) = colors(i,:);
% end
% 
% set(gca, 'XTickLabel', component_names);
% ylabel('Value');
% title(sprintf('STTC Components and Final Value'));
% ylim([-1, 1]);
% grid on;
% 
% % Add value labels on bars
% for i = 1:length(components)
%     text(i, components(i) + sign(components(i))*0.05, ...
%          sprintf('%.3f', components(i)), ...
%          'HorizontalAlignment', 'center', ...
%          'FontWeight', 'bold');
% end
% 
% % Add horizontal line at 0
% hold on;
% plot(xlim, [0, 0], 'k--', 'LineWidth', 1);
% 
% end
% 
% %% Example Usage Script
% 
% function run_sttc_example()
% % RUN_STTC_EXAMPLE - Demonstrate STTC analysis with simulated data
% 
% % Set random seed for reproducibility
% rng(42);
% 
% % Parameters
% total_time = 1000;  % ms
% dt = 30;  % ms coincidence window for 30-50ms oscillations
% 
% % Generate example spike trains with some synchrony
% % Neuron 1: Regular firing with some jitter
% base_times1 = 50:100:950;
% spikes1 = base_times1 + randn(size(base_times1)) * 5;
% 
% % Neuron 2: Similar pattern with some synchronized spikes
% base_times2 = 50:100:950;
% spikes2 = base_times2 + randn(size(base_times2)) * 5;
% 
% % Add some precisely synchronized spikes
% sync_indices = [2, 4, 6, 8];
% spikes2(sync_indices) = spikes1(sync_indices) + randn(size(sync_indices)) * 2;
% 
% % Sort spikes
% spikes1 = sort(spikes1);
% spikes2 = sort(spikes2);
% 
% fprintf('Example STTC Analysis\n');
% fprintf('=====================\n');
% fprintf('Total recording time: %.1f ms\n', total_time);
% fprintf('Coincidence window: %.1f ms\n', dt);
% fprintf('Neuron 1: %d spikes (%.1f Hz)\n', length(spikes1), ...
%         length(spikes1)/total_time*1000);
% fprintf('Neuron 2: %d spikes (%.1f Hz)\n', length(spikes2), ...
%         length(spikes2)/total_time*1000);
% 
% % Calculate STTC with significance testing
% [sttc, p_value] = compute_sttc(spikes1, spikes2, dt, total_time, ...
%                                 'n_surrogates', 1000, ...
%                                 'jitter_window', 20, ...
%                                 'method', 'jitter', ...
%                                 'verbose', true);
% 
% fprintf('\nResults:\n');
% fprintf('STTC = %.4f\n', sttc);
% fprintf('p-value = %.4f\n', p_value);
% 
% if p_value < 0.05
%     fprintf('Significant synchrony detected!\n');
% else
%     fprintf('No significant synchrony detected.\n');
% end
% 
% % Visualize
% plot_sttc_analysis(spikes1, spikes2, dt, total_time, sttc);
% 
% % Test with different coincidence windows
% fprintf('\nTesting different coincidence windows:\n');
% dt_values = [10, 20, 30, 40, 50];
% sttc_values = zeros(size(dt_values));
% 
% for i = 1:length(dt_values)
%     sttc_values(i) = compute_sttc(spikes1, spikes2, dt_values(i), total_time);
%     fprintf('  dt = %.0f ms: STTC = %.4f\n', dt_values(i), sttc_values(i));
% end
% 
% % Plot STTC vs dt
% figure;
% plot(dt_values, sttc_values, 'b.-', 'LineWidth', 2, 'MarkerSize', 20);
% xlabel('Coincidence Window (ms)');
% ylabel('STTC');
% title('STTC as a Function of Coincidence Window');
% grid on;
% 
% end
% 
% %% Validation Functions
% 
% function validate_sttc_implementation()
% % VALIDATE_STTC_IMPLEMENTATION - Test cases to validate STTC calculation
% %
% % Tests against known cases from Cutts & Eglen (2014)
% 
% fprintf('Validating STTC Implementation\n');
% fprintf('==============================\n');
% 
% % Test 1: Identical spike trains should give STTC = 1
% spikes = [10, 20, 30, 40, 50];
% sttc = compute_sttc(spikes, spikes, 5, 100);
% assert(abs(sttc - 1) < 1e-10, 'Test 1 failed: Identical trains should give STTC = 1');
% fprintf('Test 1 passed: Identical spike trains\n');
% 
% % Test 2: Non-overlapping spike trains should give STTC < 0
% spikes1 = [10, 30, 50, 70, 90];
% spikes2 = [20, 40, 60, 80, 100];
% sttc = compute_sttc(spikes1, spikes2, 5, 110);
% assert(sttc < 0, 'Test 2 failed: Non-overlapping trains should give STTC < 0');
% fprintf('Test 2 passed: Non-overlapping spike trains\n');
% 
% % Test 3: Empty spike trains
% sttc = compute_sttc([], [10, 20], 5, 100);
% assert(isnan(sttc), 'Test 3 failed: Empty train should give NaN');
% fprintf('Test 3 passed: Empty spike train handling\n');
% 
% % Test 4: Rate independence
% % Double the rate of both trains - STTC should remain similar
% spikes1 = sort(rand(20, 1) * 1000);
% spikes2 = sort(rand(20, 1) * 1000);
% sttc1 = compute_sttc(spikes1, spikes2, 50, 1000);
% 
% % Double the rate by adding more spikes
% spikes1_double = sort([spikes1; rand(20, 1) * 1000]);
% spikes2_double = sort([spikes2; rand(20, 1) * 1000]);
% sttc2 = compute_sttc(spikes1_double, spikes2_double, 50, 1000);
% 
% fprintf('Test 4: Rate independence check\n');
% fprintf('  Original rate STTC: %.4f\n', sttc1);
% fprintf('  Double rate STTC: %.4f\n', sttc2);
% fprintf('  Difference: %.4f (should be small)\n', abs(sttc1 - sttc2));
% 
% % Test 5: Perfectly synchronized trains at low rate
% spikes1 = [100, 500, 900];
% spikes2 = [100, 500, 900];
% sttc = compute_sttc(spikes1, spikes2, 10, 1000);
% assert(abs(sttc - 1) < 1e-10, 'Test 5 failed: Perfect sync should give STTC = 1');
% fprintf('Test 5 passed: Perfect synchrony at low rate\n');
% 
% fprintf('\nAll validation tests passed!\n');
% 
% end