%% Cortical Hierarchy T-statistic Analysis
% This script compares t-statistics across individual EC algorithms and iEC
% for their relationship with cortical types using GLM analysis.
%
% Main workflow:
% 1. Load group-level EC matrices for test set (individual algorithms)
% 2. Load final iEC matrix for test set
% 3. Load cortical type data
% 4. Compute cortical hierarchy for each algorithm using computeHierarchyLevels
% 5. Fit GLM models for each algorithm against cortical types
% 6. Extract t-statistics and create lollipop plot
%
% Author: Younghyun Oh
% Date: 2025-12-19
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Define constants
parcellation = 'MMP360';
rois = 360;
HIERARCHY_THRESHOLD = 0.15; % From stage2_cortical_hierarchy_main.m

% Set directories
main_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy';
ec_dir = fullfile(main_dir, 'step1_ec_validation', 'part3_empirical_human');
data_dir = fullfile(main_dir, 'data');

% Selected algorithms for comparison (VAR and FASK only)
algorithm_names = {'var','fask'};

%--------------------------------------------------------------------------%
%% Load data
%--------------------------------------------------------------------------%
fprintf('Loading data...\n');

% Load group-level EC matrices for test set (individual algorithms)
ec_alg_groups_test_full = importdata(fullfile(ec_dir, 'data', [parcellation,'_ec_alg_groups_test.mat']));

% Extract only VAR and FASK algorithms
full_algorithm_names = {'rdcm','var','fask','ccd','boss','lingam','grasp','patel'};
ec_alg_groups_test = cell(length(algorithm_names), 1);
for i = 1:length(algorithm_names)
    alg_idx = find(strcmp(full_algorithm_names, algorithm_names{i}));
    ec_alg_groups_test{i} = ec_alg_groups_test_full{alg_idx};
end

% Load final iEC matrices for test set
iEC_full= importdata(fullfile(data_dir, 'MMP360_resting_iEC.mat'));

% Load cortical type data
corticalTypes = importdata(fullfile(data_dir, 'cortical_types.mat'));

fprintf('Data loaded successfully.\n');
fprintf('Number of algorithms: %d\n', length(algorithm_names));
fprintf('Number of regions: %d\n', rois);
fprintf('Cortical types range: %d to %d\n', min(corticalTypes), max(corticalTypes));

%--------------------------------------------------------------------------%
%% Compute cortical hierarchy for each algorithm
%--------------------------------------------------------------------------%
fprintf('\nComputing cortical hierarchy levels...\n');

% Initialize storage for hierarchy levels
hierarchy_levels = cell(length(algorithm_names) + 1, 1); % +1 for iEC
algorithm_labels = [algorithm_names, {'iEC'}]; % Add iEC to the list

% Compute hierarchy for each individual algorithm
for a = 1:length(algorithm_names)
    fprintf('Computing hierarchy for %s...\n', algorithm_names{a});
    ec_matrix = ec_alg_groups_test{a};
    hierarchy_levels{a} = computeHierarchyLevels(ec_matrix, HIERARCHY_THRESHOLD);
end

% Compute hierarchy for iEC
fprintf('Computing hierarchy for iEC...\n');
hierarchy_levels{end} = computeHierarchyLevels(iEC_full, HIERARCHY_THRESHOLD);

%--------------------------------------------------------------------------%
%% Fit GLM models and extract t-statistics
%--------------------------------------------------------------------------%
fprintf('\nFitting GLM models and extracting t-statistics...\n');

% Initialize storage for t-statistics
t_statistics = zeros(length(algorithm_labels), 1);
p_values = zeros(length(algorithm_labels), 1);

% Fit GLM for each algorithm/iEC
for a = 1:length(algorithm_labels)
    fprintf('Fitting GLM for %s...\n', algorithm_labels{a});

    % Prepare data for GLM analysis (following stage2_cortical_hierarchy_main.m approach)
    hierarchy_data = [];
    cyto_group_labels = [];

    % Extract hierarchy values for each cytoarchitectonic type
    num_cyto_types = max(corticalTypes);
    for cyto_type = 1:num_cyto_types
        regions_in_type = (corticalTypes == cyto_type);
        hierarchy_values = hierarchy_levels{a}(regions_in_type);

        % Accumulate data for GLM analysis
        hierarchy_data = [hierarchy_data; hierarchy_values];
        cyto_group_labels = [cyto_group_labels; repmat(cyto_type, length(hierarchy_values), 1)];
    end

    % Fit linear model to test hierarchy-cytoarchitecture relationship
    glm_model = fitglm(cyto_group_labels, hierarchy_data, 'linear');

    % Extract t-statistic and p-value (coefficient for cytoarchitectonic type)
    t_statistics(a) = glm_model.Coefficients.tStat(2); % Second coefficient (slope)
    p_values(a) = glm_model.Coefficients.pValue(2);

    fprintf('  %s: t-stat = %.3f, p-value = %.3e\n', ...
        algorithm_labels{a}, t_statistics(a), p_values(a));
end

%--------------------------------------------------------------------------%
%% Create lollipop plot
%--------------------------------------------------------------------------%
fprintf('\nCreating lollipop plot...\n');

% Sort algorithms by t-statistic
[t_sorted, sort_idx] = sort(t_statistics, 'descend');
labels_sorted = algorithm_labels(sort_idx);

% Define colors
color_iEC = [0.6, 0.3, 0.8];     % Purple for iEC
color_VAR = [0.4, 0.4, 0.4];     % Dark grey for VAR
color_others = [0.8, 0.8, 0.8];  % Light grey for others

% Create figure
figure('Position', [100, 100, 800, 500]);

% Create lollipop plot
for i = 1:length(t_sorted)
    % Determine color based on algorithm
    if strcmp(labels_sorted{i}, 'iEC')
        plot_color = color_iEC;
        line_width = 3;
        marker_size = 8;
    elseif strcmp(labels_sorted{i}, 'var')
        plot_color = color_VAR;
        line_width = 2;
        marker_size = 6;
    else
        plot_color = color_others;
        line_width = 1.5;
        marker_size = 5;
    end

    % Plot stem (lollipop stick)
    line([i i], [0 t_sorted(i)], 'Color', plot_color, 'LineWidth', line_width);
    hold on;

    % Plot circle (lollipop head)
    scatter(i, t_sorted(i), 50, plot_color, 'filled', 'MarkerEdgeColor', 'k', ...
        'MarkerEdgeAlpha', 0.5, 'SizeData', marker_size^2);
end

% Add horizontal line at t = 1.96 (significance threshold)
yline(1.96, '--k', 'LineWidth', 2, 'Alpha', 0.7);
yline(-1.96, '--k', 'LineWidth', 2, 'Alpha', 0.7);

% Customize plot
xlabel('EC Algorithms (sorted by t-statistic)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('t-statistic', 'FontSize', 12, 'FontWeight', 'bold');
title('Cortical Hierarchy vs. Cytoarchitectonic Types: GLM t-statistics', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Set x-axis
xticks(1:length(labels_sorted));
xticklabels(upper(labels_sorted));
xtickangle(45);

% Customize appearance
grid on;
% grid minor;
set(gca, 'GridAlpha', 0.3);
set(gca, 'LineWidth', 1.5, 'FontSize', 11);
box on;

% % Add significance threshold annotation
% text(0.5, 1.96 + 0.1*range(ylim), 't = 1.96', 'FontSize', 10, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
% text(0.5, -1.96 - 0.1*range(ylim), 't = -1.96', 'FontSize', 10, ...
%     'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Set y-axis limits with some padding
y_range = max(t_sorted) - min(t_sorted);
ylim([min(t_sorted) - 0.2*y_range, max(t_sorted) + 0.2*y_range]);

% Ensure horizontal line for significance is visible
ylim_current = ylim;
if ylim_current(1) > -2.5
    ylim([min(-2.5, ylim_current(1)), ylim_current(2)]);
end
if ylim_current(2) < 2.5
    ylim([ylim_current(1), max(2.5, ylim_current(2))]);
end

hold off;

%--------------------------------------------------------------------------%
%% Save results
%--------------------------------------------------------------------------%
% Save the figure
result_dir = fullfile(main_dir, 'step2_ec_hierarchy', 'results');
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end

exportgraphics(gcf, fullfile(result_dir, 'cortical_hierarchy_tstat_lollipop.png'), 'Resolution', 1200);
savefig(gcf, fullfile(result_dir, 'cortical_hierarchy_tstat_lollipop.fig'));

% Save numerical results
results_table = table(labels_sorted', t_sorted, p_values(sort_idx), ...
    'VariableNames', {'Algorithm', 'T_Statistic', 'P_Value'});

fprintf('\n=== RESULTS SUMMARY ===\n');
disp(results_table);

% Save results to file
save(fullfile(result_dir, 'cortical_hierarchy_tstat_results.mat'), ...
    'results_table', 't_statistics', 'p_values', 'algorithm_labels', ...
    'hierarchy_levels', 'sort_idx');

fprintf('\n=== ANALYSIS COMPLETED ===\n');
fprintf('Generated files:\n');
fprintf('  - cortical_hierarchy_tstat_lollipop.png: Lollipop plot\n');
fprintf('  - cortical_hierarchy_tstat_lollipop.fig: MATLAB figure\n');
fprintf('  - cortical_hierarchy_tstat_results.mat: Numerical results\n');