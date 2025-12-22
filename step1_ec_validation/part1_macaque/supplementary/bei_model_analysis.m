%% Stage 5: Bei Model Analysis - Control vs Timeseries Comparison
% This script analyzes control vs timeseries conditions from bei_model data.
%
% Main workflow:
% 1. Load 50 iterations of control and timeseries BOLD data
% 2. Preprocess: Discard initial 10 timepoints and z-score across time
% 3. Run VAR, FASK, LiNGAM algorithms on preprocessed data
% 4. Create iEC_vfl using existing VFL betas for each iteration
% 5. Generate group-averaged iECs and correlation scatter plot
% 6. Perform edge-level statistical testing with FDR correction
% 7. Create histogram of p-values to demonstrate similarity between conditions
%
% Author: Younghyun Oh
% Date: 2025-11-18
%
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc

% Set up directories
base_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
main_dir = fullfile(base_dir, 'step1_ec_validation/part1_macaque');
bei_data_dir = fullfile(main_dir, 'bei_model/results');

% Load existing VFL betas
result_dir = fullfile(main_dir, 'results');
EC_betas_vfl = load(fullfile(result_dir, 'EC_betas.mat')).EC_betas_vfl;
median_betas_vfl = median(EC_betas_vfl);

% Load FLN
FLN_dir = fullfile(main_dir, 'data/FLN40.mat');
FLN = importdata(FLN_dir);
FLN = 1.2*FLN.^0.3; % log-linear transformation

fprintf('Loaded VFL betas: [%.4f, %.4f, %.4f]\n', median_betas_vfl);
fprintf('Starting Bei model analysis...\n');

%--------------------------------------------------------------------------%
%% Load and preprocess BOLD data
%--------------------------------------------------------------------------%

% Parameters
num_iterations = 20;
num_regions = 40;
discard_initial = 10;  % Discard first 10 timepoints

% Algorithm parameters (from stage1_run_algorithms.m)
lambda = 1e+8;
threshold = 0.1;
iter_alg = 50;

% Initialize storage for iEC results
iec_control = zeros(num_regions, num_regions, num_iterations);
iec_timeseries = zeros(num_regions, num_regions, num_iterations);

fprintf('Processing %d iterations of control and timeseries data...\n', num_iterations);

%--------------------------------------------------------------------------%
%% Process each iteration
%--------------------------------------------------------------------------%

for iter_idx = 1:num_iterations
    fprintf('Processing iteration %d/%d...\n', iter_idx, num_iterations);

    % Load control data
    control_file = fullfile(bei_data_dir, sprintf('iter_%d_BOLD_control.txt', iter_idx));
    control_bold = readmatrix(control_file);

    % Remove header line (dimensions) and discard initial timepoints
    control_bold = control_bold(:,(discard_initial+1):end);  % Discard initial timepoints

    % Z-score across time dimension
    control_bold = zscore(control_bold, 0, 2);

    % Load timeseries data
    timeseries_file = fullfile(bei_data_dir, sprintf('iter_%d_BOLD_timeseries.txt', iter_idx));
    timeseries_bold = readmatrix(timeseries_file);

    % Remove header line and discard initial timepoints
    timeseries_bold = timeseries_bold(:,(discard_initial+1):end);  % Discard initial timepoints

    % Z-score across time dimension
    timeseries_bold = zscore(timeseries_bold, 0, 2);

    %% Run EC algorithms for control condition
    fprintf('  Running EC algorithms for control condition...\n');

    % VAR
    var_control = run_var(control_bold, lambda);
    var_control_norm = var_control / max(var_control(:));

    % FASK
    fask_control = run_fask(control_bold', 1e-6, threshold, iter_alg);
    fask_control_norm = fask_control / max(fask_control(:));

    % LiNGAM
    lingam_control = run_lingam(control_bold', threshold, iter_alg);
    lingam_control_norm = lingam_control / max(lingam_control(:));

    % Create iEC_vfl for control
    iec_control(:, :, iter_idx) = (var_control_norm * median_betas_vfl(1) + ...
                                   fask_control_norm * median_betas_vfl(2) + ...
                                   lingam_control_norm * median_betas_vfl(3)) / 3;

    %% Run EC algorithms for timeseries condition
    fprintf('  Running EC algorithms for timeseries condition...\n');

    % VAR
    var_timeseries = run_var(timeseries_bold, lambda);
    var_timeseries_norm = var_timeseries / max(var_timeseries(:));

    % FASK
    fask_timeseries = run_fask(timeseries_bold', 1e-6, threshold, iter_alg);
    fask_timeseries_norm = fask_timeseries / max(fask_timeseries(:));

    % LiNGAM
    lingam_timeseries = run_lingam(timeseries_bold', threshold, iter_alg);
    lingam_timeseries_norm = lingam_timeseries / max(lingam_timeseries(:));

    % Create iEC_vfl for timeseries
    iec_timeseries(:, :, iter_idx) = (var_timeseries_norm * median_betas_vfl(1) + ...
                                      fask_timeseries_norm * median_betas_vfl(2) + ...
                                      lingam_timeseries_norm * median_betas_vfl(3)) / 3;
end

fprintf('EC analysis completed for all iterations.\n');

%--------------------------------------------------------------------------%
%% Create group-averaged iECs and correlation analysis
%--------------------------------------------------------------------------%

fprintf('Creating group-averaged iECs and correlation analysis...\n');

% Calculate group averages
group_avg_control = mean(iec_control(:,:,1:20), 3);
group_avg_timeseries = mean(iec_timeseries(:,:,1:20), 3);


% Create figure
figure('Position', [100 100 400 270]);
% Create scatter plot
scatter(group_avg_control(:), FLN(:), 40, [0.5 0.5 0.5], 'filled', ...
'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
hold on
% Add regression line
p = polyfit(group_avg_control(:), FLN(:), 1);
x_range = [min(group_avg_control(:)) max(group_avg_control(:))];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, '-r', 'LineWidth', 1.5);
% Calculate correlation
corr_val = corr(FLN(:), group_avg_control(:));
% Customize appearance
box on
grid on
set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
% Set axis limits
xlim([min(group_avg_control(:)) max(group_avg_control(:))])
ylim([min(FLN(:)) max(FLN(:))])
% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/bei_model_control_vs_FLN.png'), 'Resolution', 1200);


% Calculate correlation between group averages
correlation_group = corr(group_avg_control(:), group_avg_timeseries(:));

fprintf('Correlation between group-averaged iECs: %.4f\n', correlation_group);

% Create scatter plot of group averages
figure('Position', [100 100 400 270]);

% Create scatter plot
scatter(group_avg_control(:), group_avg_timeseries(:), 40, [0.5 0.5 0.5], 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'white', 'MarkerEdgeAlpha', 0.5);
hold on

% Add regression line
p = polyfit(group_avg_control(:), group_avg_timeseries(:), 1);
x_range = [min(group_avg_control(:)) max(group_avg_control(:))];
y_fit = polyval(p, x_range);
plot(x_range, y_fit, '-r', 'LineWidth', 1.5);


% Customize appearance
box on
grid on
% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))
xlim([min(group_avg_control(:)) max(group_avg_control(:))])
ylim([min(group_avg_timeseries(:)) max(group_avg_timeseries(:))])

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/bei_model_group_correlation.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Edge-level statistical testing
%--------------------------------------------------------------------------%

fprintf('Performing edge-level statistical testing...\n');

% Get total number of edges
total_edges = num_regions * num_regions;
fprintf('Testing %d edges for statistical differences...\n', total_edges);

% Initialize p-value storage
p_values = zeros(num_regions, num_regions);

% Perform Wilcoxon signed-rank test for each edge
for i = 1:num_regions
    for j = 1:num_regions
        % Extract data for this edge across all iterations
        control_edge_data = squeeze(iec_control(i, j, 1:20));
        timeseries_edge_data = squeeze(iec_timeseries(i, j,1:20));

        % Perform two-tailed Wilcoxon signed-rank test (signrank)
        p_values(i, j) = signrank(control_edge_data, timeseries_edge_data, 'tail', 'both');
    end
end

% Apply FDR correction
p_values_vector = p_values(:);
alpha = 0.05;
[h_corrected, crit_p, adj_ci_cvrg, p_values_corrected] = fdr_bh(p_values_vector, alpha, 'pdep', 'yes');

fprintf('FDR correction applied:\n');
fprintf('  Original significant edges (p < 0.05): %d/%d (%.2f%%)\n', ...
    sum(p_values_vector < 0.05), total_edges, 100*sum(p_values_vector < 0.05)/total_edges);
fprintf('  FDR-corrected significant edges: %d/%d (%.2f%%)\n', ...
    sum(h_corrected), total_edges, 100*sum(h_corrected)/total_edges);
fprintf('  Critical p-value: %.6f\n', crit_p);

%--------------------------------------------------------------------------%
%% Create histogram of p-values
%--------------------------------------------------------------------------%

close all

fprintf('Creating p-value heatmap...\n');

% Reshape the corrected p-values back to matrix form if necessary
if isvector(p_values_corrected)
    p_values_heatmap = reshape(p_values_corrected, num_regions, num_regions);
else
    p_values_heatmap = p_values_corrected;
end

% Fill the diagonal with 1
p_values_heatmap(eye(num_regions) == 1) = 1;

figure('Position', [100 100 270 270]);
imagesc(p_values_heatmap);
colormap(hot); % low p = white/yellow, high p = dark
% colorbar;

set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', []);
axis square;

% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/bei_model_pvalue_heatmap.png'), 'Resolution', 1200);



%--------------------------------------------------------------------------%
%% Summary statistics and visualization of edge differences
%--------------------------------------------------------------------------%

% Calculate effect sizes (Cohen's d) for each edge
effect_sizes = zeros(num_regions, num_regions);
for i = 1:num_regions
    for j = 1:num_regions
        control_data = squeeze(iec_control(i, j, :));
        timeseries_data = squeeze(iec_timeseries(i, j, :));

        pooled_std = sqrt(((var(control_data) + var(timeseries_data)) / 2));
        effect_sizes(i, j) = abs(mean(control_data) - mean(timeseries_data)) / pooled_std;
    end
end

% Create histogram of Cohen's d values
figure('Position', [100 100 600 270]);

effect_sizes_vector = effect_sizes(:);
histogram(effect_sizes_vector, 50, 'Normalization', 'probability', ...
    'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on

% set(gca, 'LineWidth', 1.5, 'TickDir', 'in', 'XTickLabel', [], 'YTickLabel', []);
box on
grid on
set(gca, 'XTick', linspace(min(get(gca, 'XTick')), max(get(gca, 'XTick')), 5))
set(gca, 'YTick', linspace(min(get(gca, 'YTick')), max(get(gca, 'YTick')), 5))


% Save figure
exportgraphics(gcf, fullfile(main_dir, 'results/bei_model_cohens_d_histogram.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Save results
%--------------------------------------------------------------------------%

% Save all results
save(fullfile(result_dir, 'bei_model_analysis_results.mat'));