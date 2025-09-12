%% Stage 2: Cortical hierarchy analysis and signal flow patterns
% This script performs cortical hierarchy analysis of the iEC matrix and generates
% components for Figure 5, examining hierarchy relationships across multiple features.
%
% Main analyses:
% 1. Hierarchy computation: Derives cortical hierarchy levels from iEC matrix
% 2. Principal gradient comparison: Correlates hierarchy with functional gradients (PC1)
% 3. Cytoarchitectonic analysis: Examines hierarchy across 6 cortical types
% 4. Modular signal flow: Computes 27-module signal flow patterns
% 5. Joy plot visualization: Creates signal flow distribution plots
%
% Output: Figure 5 components and hierarchy validation analyses
%
% Author: Younghyun Oh
% Date: 2025-06-10
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Define analysis parameters
HIERARCHY_THRESHOLD = 0.15;     % Threshold for hierarchy computation
PC1_SCALING_FACTOR = 10;        % Scaling factor for PC1 normalization
MMP360_REGIONS = 360;           % Number of MMP360 regions
NUM_CYTO_TYPES = 6;             % Number of cytoarchitectonic types
NUM_MODULES_27 = 27;            % Number of 27-module parcellation
JOY_PLOT_SPACING = 0.7;         % Spacing for joy plot visualization
SCATTER_POINT_SIZE = 50;        % Size of scatter plot points

% Define paths
data_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/data';
result_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step2_ec_hierarchy/results';
iEC_file = fullfile(data_dir, 'MMP360_resting_iEC.mat');

% Import iEC matrix
iEC = importdata(iEC_file);

%--------------------------------------------------------------------------%
%% Cortical hierarchy analysis
%--------------------------------------------------------------------------%
% Compute hierarchical levels
hierarchyLevels = computeHierarchyLevels(iEC, HIERARCHY_THRESHOLD);

% Plot the result
surfaceplot(hierarchyLevels, 'MMP360', 'both', 'viridis');
% exportgraphics(gcf, fullfile(result_dir, 'MMP360_cortical_hierarchy.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Compare against principal gradient
%--------------------------------------------------------------------------%
% Import principal gradient
data_file_PC = fullfile(data_dir, 'MMP360_PC1.mat');
pc1 = importdata(data_file_PC);

% Plot the gradient 
% surfaceplot(pc1, 'MMP', 'both', 'pc1');
% exportgraphics(gcf, fullfile(result_dir, 'Figure5c_pc1.png'), 'Resolution', 1200);

% Normalize the data for the scatter plot
pc1_range = max(pc1(:)) - min(pc1(:));  
pc1_normalized = ((pc1 - min(pc1(:))) / pc1_range) * PC1_SCALING_FACTOR;  

% Prepare data for scatter plot correlation analysis
pc1_data = pc1_normalized;  
hierarchy_data = hierarchyLevels;  

% Calculate the correlation between PC1 and hierarchy
[pc1_hierarchy_corr, p_value] = corr(pc1_data, hierarchy_data); 
fprintf('PC1-Hierarchy correlation: r = %.3f, p = %.3e\n', pc1_hierarchy_corr, p_value);

% Create the scatter plot
figure('Position', [100, 100, 350, 350]);
scatter(pc1_data, hierarchy_data, 'o', 'SizeData', SCATTER_POINT_SIZE, ...
    'MarkerFaceColor', [0.15, 0.15, 0.15], 'MarkerEdgeColor', [1, 1, 1], ...
    'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8); 
hold on;

% Calculate threshold lines (mean + standard deviation)
pc1_mean = mean(pc1_data);  
pc1_std = std(pc1_data); 
pc1_threshold = round(pc1_mean + pc1_std, 1);  

hierarchy_mean = mean(hierarchy_data); 
hierarchy_std = std(hierarchy_data);  
hierarchy_threshold = round(hierarchy_mean + hierarchy_std, 1);  

% Display threshold lines on the scatter plot
xline(pc1_threshold, '--', 'LineWidth', 1.5);  
yline(hierarchy_threshold, '--', 'LineWidth', 1.5);  

% Identify regions with divergent PC1-hierarchy relationships
higher_hierarchy_regions = find(hierarchy_data > hierarchy_threshold & pc1_data < pc1_threshold);  
scatter(pc1_data(higher_hierarchy_regions), hierarchy_data(higher_hierarchy_regions), ...
    'r', 'SizeData', SCATTER_POINT_SIZE, 'filled');  

lower_hierarchy_regions = find(hierarchy_data < hierarchy_threshold & pc1_data > pc1_threshold);  
scatter(pc1_data(lower_hierarchy_regions), hierarchy_data(lower_hierarchy_regions), ...
    'b', 'SizeData', SCATTER_POINT_SIZE, 'filled');  

% Fit and plot linear regression line
linear_model = fitlm(pc1_data, hierarchy_data);  
plot(pc1_data, linear_model.Fitted, 'k-', 'LineWidth', 2);  
set(gca, 'LineWidth', 1.5, 'TickDir', 'in')
set(gca, 'YTick', [], 'XTick', [], 'YTickLabel', [], 'XTickLabel', [])
exportgraphics(gcf, fullfile(result_dir, 'Figure5c.png'), 'Resolution', 1200);

% Create surface plots highlighting divergent regions
% Plot 1: Higher hierarchy than expected from PC1
divergent_map = zeros(MMP360_REGIONS, 1);  
divergent_map(higher_hierarchy_regions) = 1;  
surfaceplot(divergent_map, 'MMP360', 'both', 'viridis2');  
exportgraphics(gcf, fullfile(result_dir, 'Figure5c_pc1_hierarchy_difference1.png'), 'Resolution', 1200);  

% Plot 2: Lower hierarchy than expected from PC1
divergent_map = zeros(MMP360_REGIONS, 1); 
divergent_map(lower_hierarchy_regions) = 1;
surfaceplot(divergent_map, 'MMP360', 'both', 'viridis2'); 
exportgraphics(gcf, fullfile(result_dir, 'Figure5c_pc1_hierarchy_difference2.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Compare against cortical type
%--------------------------------------------------------------------------%

% Define the file path for the cortical types data
data_file_corTypes = fullfile(data_dir, 'cortical_types.mat');

% Import cortical types data from the specified file
corticalTypes = importdata(data_file_corTypes);

% Create a surface plot of the cortical types data
surfaceplot(corticalTypes, 'MMP360', 'both', 'tpl');

% Export the surface plot as a high-resolution image
exportgraphics(gcf, 'Figure5e.png', 'Resolution', 2000);

% Cytoarchitectonic type analysis
cyto_labels = {'Konicortex', 'Eulaminate-III', 'Eulaminate-II', ...
               'Eulaminate-I', 'Dysgranular', 'Agranular'};

% Generate boxplots for cytoarchitectonic groups
plotCytoGroups(hierarchyLevels, corticalTypes, cyto_labels, 'tpl', NUM_CYTO_TYPES);
savefig(fullfile(result_dir, 'cyto_hierarchy_boxplot.fig'));

% Prepare data for statistical analysis
cyto_hierarchy_data = [];
cyto_group_labels = [];

% Extract hierarchy values for each cytoarchitectonic type
for cyto_type = 1:NUM_CYTO_TYPES
    regions_in_type = (corticalTypes == cyto_type);
    hierarchy_values = hierarchyLevels(regions_in_type);
    
    fprintf('Cytoarchitectonic type %d (%s): mean hierarchy = %.3f\n', ...
        cyto_type, cyto_labels{cyto_type}, mean(hierarchy_values));
    
    % Accumulate data for GLM analysis
    cyto_hierarchy_data = [cyto_hierarchy_data; hierarchy_values];
    cyto_group_labels = [cyto_group_labels; repmat(cyto_type, length(hierarchy_values), 1)];
end

% Fit linear model to test hierarchy-cytoarchitecture relationship
cyto_glm_model = fitglm(cyto_group_labels, cyto_hierarchy_data, 'linear');
fprintf('\nCytoarchitectonic-Hierarchy GLM Results:\n');
fprintf('R-squared: %.3f\n', cyto_glm_model.Rsquared.Ordinary);
fprintf('p-value: %.3e\n', cyto_glm_model.Coefficients.pValue(2));
disp(cyto_glm_model);

%--------------------------------------------------------------------------%
%% Signal flow analysis with 27 modules
%--------------------------------------------------------------------------%
% Load 27-module parcellation
data_file_modules = fullfile(data_dir, 'module27(ordered).mat');
module_27_labels = importdata(data_file_modules);

% Visualize 27-module parcellation
surfaceplot(module_27_labels, 'MMP360', 'both', 'module27');
exportgraphics(gcf, fullfile(result_dir, 'Figure5f_brainmap.png'), 'Resolution', 2000);

% Initialize matrices for 27-module signal flow analysis
signal_flow_pos_27 = zeros(NUM_MODULES_27, NUM_MODULES_27);
signal_flow_neg_27 = zeros(NUM_MODULES_27, NUM_MODULES_27);

% Compute signal flow matrices between 27 modules in parallel
parfor module_idx = 1:NUM_MODULES_27
    [signal_flow_neg_27(:, module_idx), signal_flow_pos_27(:, module_idx)] = ...
        signalflow_modules(iEC, module_idx, '27modules');
end

% Plot the 27-module signal flow using edge bundling
plot_signal_flow(NUM_MODULES_27, signal_flow_pos_27, signal_flow_neg_27);

% Save the figure
savefig(gcf, fullfile(result_dir, 'energyflow_27module_rest.fig'));
exportgraphics(gcf, fullfile(result_dir, 'energyflow_27module_rest.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Joy plot analysis for signal flow distributions
%--------------------------------------------------------------------------%
% Compute region-level signal flow for joy plot visualization
signal_flow_pos_regions = zeros(MMP360_REGIONS, MMP360_REGIONS);
signal_flow_neg_regions = zeros(MMP360_REGIONS, MMP360_REGIONS);

% Compute signal flow matrices for all regions in parallel
parfor region_idx = 1:MMP360_REGIONS
    [signal_flow_neg_regions(:, region_idx), signal_flow_pos_regions(:, region_idx)] = ...
        signalflow_modules(iEC, region_idx, 'nomodule');
end

% Calculate outgoing signal flow for each region
outgoing_negative = abs(sum(signal_flow_neg_regions, 1));
outgoing_positive = sum(signal_flow_pos_regions, 1);

% Aggregate signal flow by 27 modules
avg_outgoing_negative = zeros(1, NUM_MODULES_27);
avg_outgoing_positive = zeros(1, NUM_MODULES_27);

for module_idx = 1:NUM_MODULES_27
    regions_in_module = (module_27_labels == module_idx);
    
    % Calculate average outgoing flow for this module
    avg_outgoing_positive(module_idx) = mean(outgoing_positive(regions_in_module));
    avg_outgoing_negative(module_idx) = mean(outgoing_negative(regions_in_module));
end

% Create joy plot visualization
figure('Position', [100, 100, 600, 300]);

% Normalize the data for visualization
neg_min = min(avg_outgoing_negative);
neg_max = max(avg_outgoing_negative);
normalized_negative = (avg_outgoing_negative - neg_min) / (neg_max - neg_min);

pos_min = min(avg_outgoing_positive);
pos_max = max(avg_outgoing_positive);
normalized_positive = (avg_outgoing_positive - pos_min) / (pos_max - pos_min);

% Create x-axis for plotting
x_axis = 1:NUM_MODULES_27;

% Add spacing for joy plot effect
normalized_positive_spaced = normalized_positive + JOY_PLOT_SPACING;

% Plot signal flow distributions
hold on;
plot(x_axis, normalized_negative, '-o', 'Color', [0.0824, 0.3765, 0.5098], ...
     'LineWidth', 1.5, 'DisplayName', 'Negative Flow');
plot(x_axis, normalized_positive_spaced, '-o', 'Color', [0.7373, 0.0078, 0.0078], ...
     'LineWidth', 1.5, 'DisplayName', 'Positive Flow');

% Fill areas under the curves
fill([x_axis fliplr(x_axis)], [zeros(1, NUM_MODULES_27) fliplr(normalized_negative)], ...
     [0.0824, 0.3765, 0.5098], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([x_axis fliplr(x_axis)], [JOY_PLOT_SPACING * ones(1, NUM_MODULES_27) fliplr(normalized_positive_spaced)], ...
     [0.7373, 0.0078, 0.0078], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Customize plot appearance
set(gca, 'XColor', 'none', 'YColor', 'none');
box off;
set(gca, 'Visible', 'off');
hold off;

% Save joy plot
savefig(gcf, fullfile(result_dir, 'Figure5f_joy.fig'));
exportgraphics(gcf, fullfile(result_dir, 'Figure5f_joy.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Summary of outputs
%--------------------------------------------------------------------------%
fprintf('\n=== CORTICAL HIERARCHY ANALYSIS COMPLETED ===\n');
fprintf('Generated files:\n');
fprintf('  - MMP360_cortical_hierarchy.png: Hierarchy surface plot\n');
fprintf('  - Figure5c.png: PC1-hierarchy correlation scatter plot\n');
fprintf('  - Figure5c_pc1_hierarchy_difference[1,2].png: Divergent region maps\n');
fprintf('  - Figure5e.png: Cytoarchitectonic type surface plot\n');
fprintf('  - cyto_hierarchy_boxplot.fig: Cytoarchitecture boxplots\n');
fprintf('  - Figure5f_brainmap.png: 27-module parcellation\n');
fprintf('  - energyflow_27module_rest.[fig,png]: Signal flow visualization\n');
fprintf('  - Figure5f_joy.[fig,png]: Joy plot of signal flow distributions\n');  
