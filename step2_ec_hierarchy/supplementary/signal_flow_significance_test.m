%% Signal Flow Significance Testing for 27 Modules
% This script tests the significance of signal flow patterns in 27 modules
% by comparing empirical iEC signal flow against null model distributions.
%
% Author: Younghyun Oh
% Date: 2025-12-19
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
clear; close all; clc;

% Define constants
NUM_MODULES_27 = 27;
NUM_NULL_ITERATIONS = 1000;
SIGNIFICANCE_ALPHA = 0.05;
MMP360_REGIONS = 360;

% Set directories
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
data_dir = fullfile(main_dir, 'data');
result_dir = fullfile(main_dir, 'step2_ec_hierarchy', 'supplementary');

if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

%--------------------------------------------------------------------------%
%% Load data
%--------------------------------------------------------------------------%
% Load iEC matrix
iEC_file = fullfile(data_dir, 'MMP360_resting_iEC.mat');
iEC = importdata(iEC_file);

% Load 27-module parcellation
data_file_modules = fullfile(data_dir, 'module27(ordered).mat');
module_27_labels = importdata(data_file_modules);

%--------------------------------------------------------------------------%
%% Compute empirical signal flow for 27 modules
%--------------------------------------------------------------------------%
% Initialize matrices for 27-module signal flow analysis
empirical_signal_flow_pos = zeros(NUM_MODULES_27, NUM_MODULES_27);
empirical_signal_flow_neg = zeros(NUM_MODULES_27, NUM_MODULES_27);

% Compute signal flow matrices between 27 modules
parfor module_idx = 1:NUM_MODULES_27
    [empirical_signal_flow_neg(:, module_idx), empirical_signal_flow_pos(:, module_idx)] = ...
        signalflow_modules(iEC, module_idx, '27modules','both');
end

% Calculate outgoing signal flow for each module
empirical_outgoing_positive = sum(empirical_signal_flow_pos, 1);
empirical_outgoing_negative = abs(sum(empirical_signal_flow_neg, 1));

%--------------------------------------------------------------------------%
%% Generate null models and compute null signal flow distributions
%--------------------------------------------------------------------------%
% Initialize storage for null distributions
null_outgoing_positive = zeros(NUM_NULL_ITERATIONS, NUM_MODULES_27);
null_outgoing_negative = zeros(NUM_NULL_ITERATIONS, NUM_MODULES_27);

parfor null_iter = 1:NUM_NULL_ITERATIONS
    % Generate null model using null_model_dir_sign
    iEC_null = null_model_dir_sign(iEC);

    % Compute signal flow for this null model
    null_signal_flow_pos = zeros(NUM_MODULES_27, NUM_MODULES_27);
    null_signal_flow_neg = zeros(NUM_MODULES_27, NUM_MODULES_27);

    % Compute signal flow matrices for null model
    for module_idx = 1:NUM_MODULES_27
        [null_signal_flow_neg(:, module_idx), null_signal_flow_pos(:, module_idx)] = ...
            signalflow_modules(iEC_null, module_idx, '27modules','both');
    end

    % Calculate outgoing signal flow for each module
    null_outgoing_positive(null_iter, :) = sum(null_signal_flow_pos, 1);
    null_outgoing_negative(null_iter, :) = abs(sum(null_signal_flow_neg, 1));
end

%--------------------------------------------------------------------------%
%% 1. Statistical Comparison: Empirical vs Null (With Inf Handling)
%--------------------------------------------------------------------------%

% Calculate null distribution statistics
% Note: Using std(..., 0, 1) for sample standard deviation (N-1 normalization)
null_pos_mean = mean(null_outgoing_positive, 1);
null_pos_std  = std(null_outgoing_positive, 0, 1);
null_neg_mean = mean(null_outgoing_negative, 1);
null_neg_std  = std(null_outgoing_negative, 0, 1);

% Calculate Z-scores
% Note: This can generate Inf if null_std is 0 (i.e., null models have 0 variance)
z_scores_pos = (empirical_outgoing_positive - null_pos_mean) ./ null_pos_std;
z_scores_neg = (empirical_outgoing_negative - null_neg_mean) ./ null_neg_std;

% --- CRITICAL FIX FOR INF / NAN VALUES ---
% If the null model has zero variance (std=0):
% Case A: Empirical != Null Mean -> Infinite Deviation (Real Signal)
% Case B: Empirical == Null Mean -> Zero Deviation (No Signal)

% Fix Positive Z-scores
z_scores_pos(isnan(z_scores_pos)) = 0; % Case B (0/0 -> 0)
max_finite_pos = max(z_scores_pos(~isinf(z_scores_pos)));
if isempty(max_finite_pos), max_finite_pos = 10; end 
% Cap Inf at 1.1x the max finite value to show it's "huge" without breaking plot
z_scores_pos(isinf(z_scores_pos)) = max_finite_pos * 1.1; 

% Fix Negative Z-scores
z_scores_neg(isnan(z_scores_neg)) = 0; % Case B (0/0 -> 0)
max_finite_neg = max(z_scores_neg(~isinf(z_scores_neg)));
if isempty(max_finite_neg), max_finite_neg = 10; end
% Cap Inf at 1.1x the max finite value
z_scores_neg(isinf(z_scores_neg)) = max_finite_neg * 1.1; 

% Calculate P-values (Two-tailed percentile method)
p_values_pos = zeros(1, NUM_MODULES_27);
p_values_neg = zeros(1, NUM_MODULES_27);

for m = 1:NUM_MODULES_27
    % Count how many nulls are more extreme than empirical
    greater_eq = sum(null_outgoing_positive(:, m) >= empirical_outgoing_positive(m));
    less_eq    = sum(null_outgoing_positive(:, m) <= empirical_outgoing_positive(m));
    p_values_pos(m) = 2 * min(greater_eq, less_eq) / NUM_NULL_ITERATIONS;

    greater_eq = sum(null_outgoing_negative(:, m) >= empirical_outgoing_negative(m));
    less_eq    = sum(null_outgoing_negative(:, m) <= empirical_outgoing_negative(m));
    p_values_neg(m) = 2 * min(greater_eq, less_eq) / NUM_NULL_ITERATIONS;
end

%--------------------------------------------------------------------------%
%% 2. Visualization: Aesthetic Lollipop Plot (Minimalist)
%--------------------------------------------------------------------------%

% Set figure size
figure('Position', [100, 100, 1000, 600], 'Color', 'w');

% Define Colors (MATLAB R2019a+ supports hex codes)
color_pos = '#bc0202'; % Muted Red
color_neg = '#2121cf'; % Muted Blue
color_gray = [0.6 0.6 0.6]; 

%-----------------------------------%
% Subplot 1: Positive Flow Z-Scores
%-----------------------------------%
subplot(2, 1, 1);
hold on;

% 1. Add Significance Threshold Lines (Behind data)
yline(1.96, '--', 'Color', color_gray, 'LineWidth', 1.5);
yline(-1.96, '--', 'Color', color_gray, 'LineWidth', 1.5);
yline(0, '-', 'Color', 'k', 'LineWidth', 2);

% 2. Create Lollipop Plot
h1 = stem(z_scores_pos, 'filled');

% 3. Style
h1.Color = color_pos;          
h1.LineWidth = 1.5;            
h1.MarkerSize = 6;             
h1.MarkerFaceColor = color_pos; 
h1.MarkerEdgeColor = color_pos;


% 4. Sparse Ticks
xlim([0, NUM_MODULES_27 + 1]);
box off;

% X-Axis
xticks(ceil(linspace(1,27,5))); 
set(gca, 'XTickLabel', []); % Uncomment to hide numbers

% Y-Axis (Min, Zero, Max)
ylim_curr = ylim;
yticks([floor(ylim_curr(1)), 0, ceil(ylim_curr(2))]); 
% set(gca, 'YTickLabel', []); % Uncomment to hide numbers

% Polish
set(gca, 'LineWidth', 1.5, 'TickDir', 'in');

%-----------------------------------%
% Subplot 2: Negative Flow Z-Scores
%-----------------------------------%
subplot(2, 1, 2);
hold on;

% 1. Significance Lines
yline(1.96, '--', 'Color', color_gray, 'LineWidth', 1.5);
yline(-1.96, '--', 'Color', color_gray, 'LineWidth', 1.5);
yline(0, '-', 'Color', 'k', 'LineWidth', 2);

% 2. Create Lollipop Plot
h2 = stem(z_scores_neg, 'filled');

% 3. Style
h2.Color = color_neg;
h2.LineWidth = 1.5;
h2.MarkerSize = 6;
h2.MarkerFaceColor = color_neg;
h2.MarkerEdgeColor = color_neg;

% 4. Sparse Ticks
xlim([0, NUM_MODULES_27 + 1]);
box off;

% X-Axis
xticks(ceil(linspace(1,27,5))); 
set(gca, 'XTickLabel', []); 

% Y-Axis (Zero, Mid, Max)
max_val = max(z_scores_neg);
yticks([0, round(max_val/2), round(max_val)]); 
ytickformat('%.1e'); % Scientific notation for huge Z-scores
% set(gca, 'YTickLabel', []); 

% Polish
set(gca, 'LineWidth', 1.5, 'TickDir', 'in');
exportgraphics(gcf, fullfile(result_dir, 'signal_flow_sig_test.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Save results
%--------------------------------------------------------------------------%
save(fullfile(result_dir, 'signal_flow_significance_results.mat'));

