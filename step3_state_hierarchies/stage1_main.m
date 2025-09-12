%% Stage 1: State-dependent cortical hierarchy analysis for Figure 6
% This script analyzes how cortical hierarchy changes across different brain states
% (resting, movie-watching, pain) and generates Figure 6 components.
%
% Main analyses:
% 1. Hierarchy computation: Derives hierarchy levels for each brain state
% 2. State comparison: Visualizes hierarchies across cortical zones (Figure 6a-b)
% 3. Signal flow analysis: Computes 27-module signal flow for each state (Figure 6c)
% 4. Delta analysis: Statistical comparison of state differences with significance testing
% 5. Joy plot visualization: Shows state-dependent signal flow changes
%
% Output: Figure 6 components demonstrating state-dependent hierarchy reorganization
%
% Author: Younghyun Oh
% Date: 2025-06-09
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Define analysis parameters
HIERARCHY_THRESHOLD = 0.15;       % Threshold for hierarchy computation
MMP360_REGIONS = 360;             % Number of MMP360 regions
NUM_CORTICAL_ZONES = 4;           % Number of cortical zones
NUM_MODULES_27 = 27;              % Number of 27-module parcellation
JOY_PLOT_SPACING = 7;             % Joy plot vertical spacing
SIGNIFICANCE_THRESHOLD_1 = 1.645; % Z-score threshold for p<0.1
SIGNIFICANCE_THRESHOLD_2 = 1.96;  % Z-score threshold for p<0.05
FDR_ALPHA = 0.05;                 % FDR significance level
FIGURE_RESOLUTION = 1200;         % Export resolution

% Define paths
main_path = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step3_state_hierarchies';
data_path = fullfile(main_path, 'data');
result_path = fullfile(main_path, 'results');

% Load ECs
EC_rest = importdata(fullfile(data_path, 'MMP360_resting_iEC.mat'));
EC_movie = importdata(fullfile(data_path,'MMP360_movie_iEC.mat'));
EC_pain = importdata(fullfile(data_path,'MMP360_pain_iEC.mat'));

%--------------------------------------------------------------------------%
%% Figure 6a: Hierarchy computation across brain states
%--------------------------------------------------------------------------%
% Compute hierarchical levels for each brain state
hierarchyRest = computeHierarchyLevels(EC_rest, HIERARCHY_THRESHOLD);
hierarchyMovie = computeHierarchyLevels(EC_movie, HIERARCHY_THRESHOLD);
hierarchyPain = computeHierarchyLevels(EC_pain, HIERARCHY_THRESHOLD);

% Plot hierarchy results for each state
surfaceplot(hierarchyMovie, 'MMP', 'both', 'viridis');
exportgraphics(gcf, fullfile(result_path, 'Figure6aMovie.png'), 'Resolution', FIGURE_RESOLUTION);

surfaceplot(hierarchyPain, 'MMP', 'both', 'viridis');
exportgraphics(gcf, fullfile(result_path, 'Figure6aPain.png'), 'Resolution', FIGURE_RESOLUTION);

%--------------------------------------------------------------------------%
%% Figure 6b: Cortical zone analysis
%--------------------------------------------------------------------------%
% Import Cortical Zones Atlas
CortZones = importdata(fullfile(data_path, 'cortical_zones.mat'));

% Uncomment to generate cortical zones brain map
% surfaceplot(CortZones, 'MMP', 'both', 'zone4');
% exportgraphics(gcf, fullfile(result_path, 'Figure6b_brainmap.png'), 'Resolution', 2000);

% Initialize arrays for hierarchy statistics
rest_ratio = zeros(NUM_CORTICAL_ZONES, 1);
pain_ratio = zeros(NUM_CORTICAL_ZONES, 1);
movie_ratio = zeros(NUM_CORTICAL_ZONES, 1);

rest_sem = zeros(NUM_CORTICAL_ZONES, 1);
pain_sem = zeros(NUM_CORTICAL_ZONES, 1);
movie_sem = zeros(NUM_CORTICAL_ZONES, 1);

% Calculate median hierarchy levels and SEM for each cortical zone
for cyto_type = 1:NUM_CORTICAL_ZONES
    indx = CortZones == cyto_type;
    rest_data = hierarchyRest(indx);
    pain_data = hierarchyPain(indx);
    movie_data = hierarchyMovie(indx);
    
    rest_ratio(cyto_type) = median(rest_data);
    pain_ratio(cyto_type) = median(pain_data);
    movie_ratio(cyto_type) = median(movie_data);
    
    % Calculate SEM
    rest_sem(cyto_type) = std(rest_data) / sqrt(numel(rest_data));
    pain_sem(cyto_type) = std(pain_data) / sqrt(numel(pain_data));
    movie_sem(cyto_type) = std(movie_data) / sqrt(numel(movie_data));
end

% Create hierarchy comparison plot across cortical zones
figure;
errorbar(1:NUM_CORTICAL_ZONES, rest_ratio, rest_sem, '-o', 'Color', '#ED6A5A', 'LineWidth', 2); hold on;
errorbar(1:NUM_CORTICAL_ZONES, movie_ratio, movie_sem, '-o', 'Color', '#9BC1BC', 'LineWidth', 2); hold on;
errorbar(1:NUM_CORTICAL_ZONES, pain_ratio, pain_sem, '-o', 'Color', '#B8B8FF', 'LineWidth', 2);

% Customize plot appearance
box off;
set(gca, 'LineWidth', 1.5, 'TickDir', 'in');
set(gca, 'YTick', [4, 5, 6, 7]);
set(gca, 'YTickLabel', [], 'XTickLabel', []);
savefig(fullfile(result_path, 'Figure6b.fig'));

%--------------------------------------------------------------------------%
%% Figure 6c: 27-module signal flow analysis
%--------------------------------------------------------------------------%
% ----------------------------- Pain Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_pain = zeros(NUM_MODULES_27, NUM_MODULES_27);
signal_flow_neg_pain = zeros(NUM_MODULES_27, NUM_MODULES_27);

% Compute signal flow matrices in parallel
parfor i = 1:NUM_MODULES_27
    [signal_flow_neg_pain(:, i), signal_flow_pos_pain(:, i)] = signalflow_modules(EC_pain, i, '27modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(NUM_MODULES_27, signal_flow_pos_pain, signal_flow_neg_pain);

% Save the figure
savefig(gcf, fullfile(result_path, 'energyflow_27module_pain.fig'));

% ----------------------------- Movie Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow
signal_flow_pos_movie = zeros(NUM_MODULES_27, NUM_MODULES_27);
signal_flow_neg_movie = zeros(NUM_MODULES_27, NUM_MODULES_27);

% Compute signal flow matrices in parallel
parfor i = 1:NUM_MODULES_27
    [signal_flow_neg_movie(:, i), signal_flow_pos_movie(:, i)] = signalflow_modules(EC_movie, i, '27modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(NUM_MODULES_27, signal_flow_pos_movie, signal_flow_neg_movie);

% Save the figure
savefig(gcf, fullfile(result_path, 'energyflow_27module_movie.fig'));

%--------------------------------------------------------------------------%
%% Delta analysis: Signal flow state differences
%--------------------------------------------------------------------------%

% Initialize matrices for positive and negative signal flow (region-level)
signal_flow_pos_rest = zeros(MMP360_REGIONS, MMP360_REGIONS);
signal_flow_neg_rest = zeros(MMP360_REGIONS, MMP360_REGIONS);

% Compute signal flow matrices in parallel
parfor i = 1:MMP360_REGIONS
    [signal_flow_neg_rest(:, i), signal_flow_pos_rest(:, i)] = signalflow_modules(EC_rest, i, 'nomodule');
end

% ----------------------------- Pain Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow (region-level)
signal_flow_pos_pain = zeros(MMP360_REGIONS, MMP360_REGIONS);
signal_flow_neg_pain = zeros(MMP360_REGIONS, MMP360_REGIONS);

% Compute signal flow matrices in parallel
parfor i = 1:MMP360_REGIONS
    [signal_flow_neg_pain(:, i), signal_flow_pos_pain(:, i)] = signalflow_modules(EC_pain, i, 'nomodule');
end


% ----------------------------- Movie Signal Flow ----------------------------- %

% Initialize matrices for positive and negative signal flow (region-level)
signal_flow_pos_movie = zeros(MMP360_REGIONS, MMP360_REGIONS);
signal_flow_neg_movie = zeros(MMP360_REGIONS, MMP360_REGIONS);

% Compute signal flow matrices in parallel
parfor i = 1:MMP360_REGIONS
    [signal_flow_neg_movie(:, i), signal_flow_pos_movie(:, i)] = signalflow_modules(EC_movie, i, 'nomodule');
end
%--------------------------------------------------------------------------%
%% Module averaging and statistical analysis
%--------------------------------------------------------------------------%
% Load the 27 modules atlas
modules = importdata(fullfile(data_path, 'module27(ordered).mat'));
unique_modules = unique(modules);

% Initialize variables to store average and raw signal flow values for each module
avg_signal_flow_pos_rest = zeros(length(unique_modules), 1);
avg_signal_flow_neg_rest = zeros(length(unique_modules), 1);
avg_signal_flow_pos_pain = zeros(length(unique_modules), 1);
avg_signal_flow_neg_pain = zeros(length(unique_modules), 1);
avg_signal_flow_pos_movie = zeros(length(unique_modules), 1);
avg_signal_flow_neg_movie = zeros(length(unique_modules), 1);

raw_signal_flow_pos_rest = cell(length(unique_modules), 1);
raw_signal_flow_neg_rest = cell(length(unique_modules), 1);
raw_signal_flow_pos_pain = cell(length(unique_modules), 1);
raw_signal_flow_neg_pain = cell(length(unique_modules), 1);
raw_signal_flow_pos_movie = cell(length(unique_modules), 1);
raw_signal_flow_neg_movie = cell(length(unique_modules), 1);


% Calculate the average and store raw signal flow values for each module
for module = 1:length(unique_modules)
    
    % Extract indices of each module
    moduleIndx = (modules == unique_modules(module));
    
    % Extract raw signal flow values for each state
    raw_signal_flow_pos_rest{module} = signal_flow_pos_rest(:, moduleIndx);
    raw_signal_flow_neg_rest{module} = signal_flow_neg_rest(:, moduleIndx);
    
    raw_signal_flow_pos_pain{module} = signal_flow_pos_pain(:, moduleIndx);
    raw_signal_flow_neg_pain{module} = signal_flow_neg_pain(:, moduleIndx);
    
    raw_signal_flow_pos_movie{module} = signal_flow_pos_movie(:, moduleIndx);
    raw_signal_flow_neg_movie{module} = signal_flow_neg_movie(:, moduleIndx);
    
    % Compute average signal flow values for each state
    avg_signal_flow_pos_rest(module) = mean(raw_signal_flow_pos_rest{module}(:));
    avg_signal_flow_neg_rest(module) = mean(raw_signal_flow_neg_rest{module}(:));
    
    avg_signal_flow_pos_pain(module) = mean(raw_signal_flow_pos_pain{module}(:));
    avg_signal_flow_neg_pain(module) = mean(raw_signal_flow_neg_pain{module}(:));
    
    avg_signal_flow_pos_movie(module) = mean(raw_signal_flow_pos_movie{module}(:));
    avg_signal_flow_neg_movie(module) = mean(raw_signal_flow_neg_movie{module}(:));
    
end

% Calculate state differences for statistical analysis
% Rest vs Movie comparisons
diff_avg_signal_flow_pos_rest_movie = avg_signal_flow_pos_movie - avg_signal_flow_pos_rest;
diff_avg_signal_flow_neg_rest_movie = abs(avg_signal_flow_neg_movie) - abs(avg_signal_flow_neg_rest);

% Rest vs Pain comparisons
diff_avg_signal_flow_pos_rest_pain = avg_signal_flow_pos_pain - avg_signal_flow_pos_rest;
diff_avg_signal_flow_neg_rest_pain = abs(avg_signal_flow_neg_pain) - abs(avg_signal_flow_neg_rest);

%--------------------------------------------------------------------------%
%% Joy plot visualization for state differences
%--------------------------------------------------------------------------%
figure;

% Extract signal flow differences for rest vs pain comparison
vector1 = diff_avg_signal_flow_neg_rest_pain; 
vector2 = diff_avg_signal_flow_pos_rest_pain; 

% Normalize the data to retain negative values
normalized_neg = (vector1 - mean(vector1)) / std(vector1);

% Compute statistical thresholds for negative signal flow
mu_neg = mean(normalized_neg);
sigma_neg = std(normalized_neg);
threshold_low_neg1 = mu_neg - SIGNIFICANCE_THRESHOLD_1 * sigma_neg;  % Lower threshold for p<0.1
threshold_high_neg1 = mu_neg + SIGNIFICANCE_THRESHOLD_1 * sigma_neg; % Upper threshold for p<0.1
threshold_low_neg2 = mu_neg - SIGNIFICANCE_THRESHOLD_2 * sigma_neg;  % Lower threshold for p<0.05
threshold_high_neg2 = mu_neg + SIGNIFICANCE_THRESHOLD_2 * sigma_neg; % Upper threshold for p<0.05

% Normalize positive signal flow data
normalized_pos = (vector2 - mean(vector2)) / std(vector2);

% Compute statistical thresholds for positive signal flow
mu_pos = mean(normalized_pos);
sigma_pos = std(normalized_pos);
threshold_low_pos1 = mu_pos - SIGNIFICANCE_THRESHOLD_1 * sigma_pos;  % Lower threshold for p<0.1
threshold_high_pos1 = mu_pos + SIGNIFICANCE_THRESHOLD_1 * sigma_pos; % Upper threshold for p<0.1
threshold_low_pos2 = mu_pos - SIGNIFICANCE_THRESHOLD_2 * sigma_pos;  % Lower threshold for p<0.05
threshold_high_pos2 = mu_pos + SIGNIFICANCE_THRESHOLD_2 * sigma_pos; % Upper threshold for p<0.05

% Statistical significance testing
% Convert z-scores to p-values (two-tailed test)
p_neg = 2 * (1 - normcdf(abs(normalized_neg)));
p_pos = 2 * (1 - normcdf(abs(normalized_pos)));

% Apply False Discovery Rate (FDR) correction using Benjamini-Hochberg procedure
[~, ~, ~, fdr_p_neg] = fdr_bh(p_neg);
[~, ~, ~, fdr_p_pos] = fdr_bh(p_pos);

% Add spacing between curves for joy plot visualization
normalized_pos = normalized_pos + JOY_PLOT_SPACING;

% Define x-axis data points
xData = 1:NUM_MODULES_27;

% Create interpolated data for smooth visualization
num_interp_points = 500;
x_interp = linspace(min(xData), max(xData), num_interp_points);
normalized_neg_interp = interp1(xData, normalized_neg, x_interp, 'linear');
normalized_pos_interp = interp1(xData, normalized_pos, x_interp, 'linear');

% Define color scheme for visualization
color_negative = [0.0824, 0.3765, 0.5098];  % Blue for negative flow
color_positive = [0.7373, 0.0078, 0.0078];  % Red for positive flow


% Create the joy plot
hold on;
plot(xData, normalized_neg, '-o', 'Color', color_negative, 'LineWidth', 1.5);
plot(xData, normalized_pos, '-o', 'Color', color_positive, 'LineWidth', 1.5);

% Add reference lines at zero for each distribution
plot(x_interp, zeros(size(x_interp)), '--k', 'LineWidth', 1);
plot(x_interp, JOY_PLOT_SPACING * ones(size(x_interp)), '--k', 'LineWidth', 1); 

hold off;

% Add symbols based on thresholds for significance
% Add significance markers based on statistical thresholds
for i = 1:length(xData)
    % Significance level p<0.1 (*)
    if (normalized_neg(i) < threshold_low_neg1 || normalized_neg(i) > threshold_high_neg1) && ...
            ~(normalized_neg(i) < threshold_low_neg2 || normalized_neg(i) > threshold_high_neg2)
        text(xData(i), normalized_neg(i) + 0.1, '*', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if (normalized_pos(i) < threshold_low_pos1 + JOY_PLOT_SPACING || normalized_pos(i) > threshold_high_pos1 + JOY_PLOT_SPACING) && ...
            ~(normalized_pos(i) < threshold_low_pos2 + JOY_PLOT_SPACING || normalized_pos(i) > threshold_high_pos2 + JOY_PLOT_SPACING)
        text(xData(i), normalized_pos(i) + 0.1, '*', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end

    % Significance level p<0.05 (**)
    if (normalized_neg(i) < threshold_low_neg2 || normalized_neg(i) > threshold_high_neg2) && ...
            fdr_p_neg(i) >= FDR_ALPHA
        text(xData(i), normalized_neg(i) + 0.1, '**', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if (normalized_pos(i) < threshold_low_pos2 + JOY_PLOT_SPACING || normalized_pos(i) > threshold_high_pos2 + JOY_PLOT_SPACING) && ...
            fdr_p_pos(i) >= FDR_ALPHA
        text(xData(i), normalized_pos(i) + 0.1, '**', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    
    % FDR corrected significance (***)
    if fdr_p_neg(i) < FDR_ALPHA
        text(xData(i), normalized_neg(i) + 0.1, '***', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
    if fdr_p_pos(i) < FDR_ALPHA
        text(xData(i), normalized_pos(i) + 0.1, '***', 'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
    end
end

% Finalize joy plot appearance
set(gca, 'XColor', 'none', 'YColor', 'none');
box off;
set(gca, 'Visible', 'off');
hold off;

% Save the joy plot
savefig(fullfile(result_path, 'Figure6_joy_plot_rest_vs_pain.fig'));
exportgraphics(gcf, fullfile(result_path, 'Figure6_joy_plot_rest_vs_pain.png'), 'Resolution', FIGURE_RESOLUTION);

%--------------------------------------------------------------------------%
%% Summary of outputs
%--------------------------------------------------------------------------%
fprintf('\n=== STATE HIERARCHIES ANALYSIS COMPLETED ===\n');
fprintf('Generated files:\n');
fprintf('  - Figure6aMovie.png: Movie state hierarchy surface plot\n');
fprintf('  - Figure6aPain.png: Pain state hierarchy surface plot\n');
fprintf('  - Figure6b.fig: Cortical zone hierarchy comparisons\n');
fprintf('  - energyflow_27module_pain.fig: Pain state signal flow\n');
fprintf('  - energyflow_27module_movie.fig: Movie state signal flow\n');
fprintf('  - Figure6_joy_plot_rest_vs_pain.[fig,png]: State difference joy plot\n');
fprintf('\nAnalyses demonstrate state-dependent reorganization of cortical hierarchy\n');
fprintf('and signal flow patterns across brain states.\n');
