%% Stage 1: Network profiling and signal flow analysis for Figure 4
% This script performs comprehensive network analysis of the iEC matrix to 
% generate Figure 4 components showing network structure and signal flow patterns.
%
% Main analyses:
% 1. Matrix visualization: Displays full iEC matrix with sorted network organization
% 2. Edge distribution analysis: Examines strength distribution with tail analysis
% 3. Network degree analysis: Computes in/out-flow and positive/negative connection ratios
% 4. Modular boxplot analysis: Analyzes connection patterns across 7 functional networks
% 5. Signal flow analysis: Computes directed signal flow between 22 cortical modules
% 6. Cross-module flow analysis: Examines unimodal vs heteromodal signal flow patterns
%
% Output: Figure 4a-e components and supplementary network visualizations
%
% Author: Younghyun Oh  
% Date: 2025-06-10
%--------------------------------------------------------------------------%
%% Set up
%--------------------------------------------------------------------------%
% Clear workspace
clear; close all; clc;

% Define analysis parameters
TAIL_ANALYSIS_PERCENT = 0.05;  % Top 5% for tail analysis
NUM_HISTOGRAM_BINS = 100;      % Number of bins for histogram
LEFT_CORTEX_REGIONS = 180;     % Number of left cortex regions
NUM_YEONETWORKS = 7;           % Number of Yeo networks
NUM_MODULES_22 = 22;           % Number of 22-module parcellation

% Define paths
data_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/data';
iEC_file = fullfile(data_dir, 'MMP360_resting_iEC.mat');
module_file = fullfile(data_dir, 'Yeo7MMP.mat');
result_dir = '/combinelab/03_user/younghyun/01_project/01_Signalflow_Hierarchy/step2_ec_hierarchy/results';

% Import iEC matrix
iEC = importdata(iEC_file);

% Import module data
modules = importdata(module_file);

%--------------------------------------------------------------------------%
%% Matrix visualization
%--------------------------------------------------------------------------%
% Flatten matrix into a vector
X = iEC(:);

% Uniform weights (since no specific weighting is given)
weights = ones(size(iEC));

% Sort X in descending order
[X_sorted, sort_idx] = sort(X, 'descend');

% Flatten the weight matrix into a vector
W = weights(:);

% Sort weights to match the sorted X values
W_sorted = W(sort_idx);

% Choose the number of largest observations (k)
n = length(X_sorted);
k = round(n * TAIL_ANALYSIS_PERCENT);  % Select top percentage for tail analysis

% Call the TailWHill function
[Hill, Hillsd, DJV1, DJV2, AM, AMsd, T1, T1sd, T2, T3, D, Dsd] = TailWHill(X_sorted, W_sorted, k);

% Extract left cortex of iEC matrix
iEC_left = iEC(1:LEFT_CORTEX_REGIONS, 1:LEFT_CORTEX_REGIONS);
modules_left = modules(1:LEFT_CORTEX_REGIONS);

% Sort out the left cortex according the the 7 networks
[~, sortedIndices] = sort(modules_left);
iEC_sorted = iEC_left(sortedIndices, sortedIndices);

% Plot sorted iEC matrix
figure('Position', [100, 100, 250, 250]);
imagesc(iEC);
colormap(generateColorMap(iEC(:), 1000));
axis square;
set(gca, 'XTick', [], 'YTick', []);
exportgraphics(gcf, fullfile(result_dir, 'Figure4a.png'), 'Resolution', 1200);

%-------------------------------------------------------------------------%
%% Histogram showing edge strength distribution
%-------------------------------------------------------------------------%

figure('Position', [100, 100, 300, 270]);
histogram(iEC(:), NUM_HISTOGRAM_BINS, 'FaceColor', '#808080', 'Normalization', 'probability');
set(gca, 'XTick', [0 0.2 0.4 0.6], 'YTick', [0 0.2 0.4 0.6]);
set(gca, 'LineWidth', 1.5, 'XTickLabel', [], 'YTickLabel', []);
box off

% Main axis limits
mainAx = gca;
mainAxLimits = axis(mainAx);

% Create an inset axes on the bottom part of the current axes for histogram
% Position inset relative to the new figure size
insetAxHist = axes('Position', [0.55 0.22 0.28 0.1]); 
% Histogram for the inset, focusing on a specific range
histogram(insetAxHist, iEC_left(:), 'BinLimits', [0.4, 0.8], 'FaceColor', '#808080', 'Normalization', 'probability');

% Set the axis limits and ticks for the histogram inset
set(insetAxHist, 'XLim', [0.4, 0.8], 'YLim', [0, 0.002]);
set(insetAxHist, 'XTick', [0.4, 0.8], 'YTick', [0, 0.002]);
set(insetAxHist, 'XTickLabel', [], 'YTickLabel', []);
box off

% Create an inset axes on the upper part of the current axes for bar plot
% Position inset relative to the new figure size
insetAxBar = axes('Position', [0.55 0.5 0.3 0.3]); 

% Calculate proportions for the bar plot
totalConnections = numel(nonzeros(iEC(:)));
numNegativeConnections = sum(nonzeros(iEC(:)) < 0);
numPositiveConnections = sum(nonzeros(iEC(:)) > 0);

% Bar plot for the inset
bar(insetAxBar, [1, 2], [numNegativeConnections, numPositiveConnections] / totalConnections, 'FaceColor', 'flat', ...
    'CData', [0 0 1; 1 0 0]); % Blue for negative, red for positive
set(insetAxBar, 'XTickLabel', [], 'XTick', [1 2], 'ylim', [0 0.75], 'YTick', [0.3 0.6], 'YTickLabel', []);
grid on;
box on

% Export the figure
exportgraphics(gcf, fullfile(result_dir, 'Figure4b.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Network degree and ratio of positive and negative connections
%--------------------------------------------------------------------------%
% Calculate outflow and inflow
outflow = sum(abs(iEC), 1);
inflow = sum(abs(iEC), 2);

% Plot and export outflow
surfaceplot(outflow, 'MMP360', 'both', 'viridis');
exportgraphics(gcf, fullfile(result_dir, 'Supple_Outflow.png'), 'Resolution', 2000);

% Plot and export inflow
surfaceplot(inflow, 'MMP360', 'both', 'viridis');
exportgraphics(gcf, fullfile(result_dir, 'Supple_Inflow.png'), 'Resolution', 2000);

% Number of nodes in the connectivity matrix 'EC'
num_nodes = size(iEC, 2); 

% Initialize arrays to store the ratios for outgoing and incoming connections
ratios_out = zeros(2, num_nodes); % Ratio of positive outgoing connections
ratios_in = zeros(2, num_nodes); % Ratio of positive incoming connections

% Calculate the ratios for each node
for i = 1:num_nodes
    % Outgoing connections for node i (column-wise)
    outgoing_connections = iEC(:, i);
    ratios_out(1, i) = sum(outgoing_connections > 0) / sum(outgoing_connections ~= 0);
    ratios_out(2, i) = sum(outgoing_connections < 0) / sum(outgoing_connections ~= 0);
    
    % Incoming connections for node i (row-wise)
    incoming_connections = iEC(i, :);
    ratios_in(1, i) = sum(incoming_connections > 0) / sum(incoming_connections ~= 0);
    ratios_in(2, i) = sum(incoming_connections < 0) / sum(incoming_connections ~= 0);
end

% Aggregate and plot ratio differences
ratio_outgoing_diff = ratios_out(1,:) - ratios_out(2,:);
surfaceplot(ratio_outgoing_diff, 'MMP360', 'both', 'BlRd');
exportgraphics(gcf, fullfile(result_dir, 'Supple_OutgoingRatio.png'), 'Resolution', 2000);

ratio_incoming_diff = ratios_in(1,:) - ratios_in(2,:);
surfaceplot(ratio_incoming_diff, 'MMP360', 'both', 'BlRd');
exportgraphics(gcf, fullfile(result_dir, 'Supple_IncomingRatio.png'), 'Resolution', 2000);

% Organize data for modular analysis
module_order = [2 1 4 3 5 6 7]; % SM, Vis, DAN, VAN, Limbic, FPN, DMN
num_modules = length(module_order);

% Pre-allocate for efficiency
total_regions = length(modules);
module_values_pos_out = zeros(total_regions, 1);
module_values_neg_out = zeros(total_regions, 1);
module_nums = zeros(total_regions, 1);

% Collect data efficiently
start_idx = 1;
for m = 1:num_modules
    current_module = module_order(m);
    module_mask = (modules == current_module);
    num_regions_in_module = sum(module_mask);
    end_idx = start_idx + num_regions_in_module - 1;
    
    % Extract values for current module
    module_values_pos_out(start_idx:end_idx) = ratios_out(1, module_mask)';
    module_values_neg_out(start_idx:end_idx) = ratios_out(2, module_mask)';
    module_nums(start_idx:end_idx) = m;
    
    start_idx = end_idx + 1;
end

% Create structured data table for Figure 4c
values = [module_values_pos_out; module_values_neg_out];
connection_type = [repmat({'posOut'}, total_regions, 1); repmat({'negOut'}, total_regions, 1)];
module_labels = [module_nums; module_nums];

ratio_array = table(values, module_labels, connection_type, ...
    'VariableNames', {'Var1', 'Var2', 'reshaped_labels'});

% Save the resulting table to a .mat file
save(fullfile(result_dir, 'Figure4c.mat'), 'ratio_array');

% Prepare data for boxplot visualization
var2_names = {'SM', 'Vis', 'DAN', 'VAN', 'Limbic', 'FPN', 'DMN'};
connection_types = {'posOut', 'negOut'};

% Pre-allocate for efficiency  
data = [];
labels = {};

% Create organized data for boxplot
for module_idx = 1:length(var2_names)
    for conn_idx = 1:length(connection_types)
        % Find matching data
        matching_rows = (ratio_array.Var2 == module_idx) & ...
                       strcmp(ratio_array.reshaped_labels, connection_types{conn_idx});
        
        % Extract data and create labels
        current_data = ratio_array.Var1(matching_rows);
        current_label = sprintf('%s_%s', var2_names{module_idx}, connection_types{conn_idx});
        
        % Append to arrays
        data = [data; current_data];
        labels = [labels; repmat({current_label}, length(current_data), 1)];
    end
end

% Create the boxplot and set box colors
figure('Position', [100, 100, 600, 270]);
h = boxplot(data, labels, 'Colors', 'k', 'Symbol', ''); % 'k' sets the box color to black

% Define colors for each box
col_vals = [...
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00';
    '21'; '21'; 'cf';
    'bc'; '08'; '00'];
col_vals = reshape(hex2dec(col_vals), [3 numel(col_vals)/6])' ./ 255;
P = size(col_vals, 1);
colors = interp1(1:size(col_vals, 1), col_vals, linspace(1, P, 14), 'linear');

% Set the colors for each box
hGroups = findobj(gca, 'Tag', 'Box');
for j = 1:length(hGroups)
    patch(get(hGroups(j), 'XData'), get(hGroups(j), 'YData'), colors(j, :), 'FaceAlpha', 0.8); % Fill the box with color
end

% Hold the current plot
hold on;

% Add scatter plot points in grey color
unique_groups = {...
    'SM_posOut', 'SM_negOut', ...
    'Vis_posOut', 'Vis_negOut', ...
    'DAN_posOut', 'DAN_negOut', ...
    'VAN_posOut', 'VAN_negOut', ...
    'Limbic_posOut', 'Limbic_negOut', ...
    'FPN_posOut', 'FPN_negOut', ...
    'DMN_posOut', 'DMN_negOut'};
num_groups = length(unique_groups);

% Adjust jitter for better visibility
jitterAmount = 0.1;

for i = 1:num_groups
    % Find the x-position for each group
    x_pos = repmat(i, sum(strcmp(labels, unique_groups{i})), 1) + (rand(sum(strcmp(labels, unique_groups{i})), 1) - 0.5) * jitterAmount;
    y_pos = data(strcmp(labels, unique_groups{i}));
    scatter(x_pos, y_pos, 5, 'filled', 'MarkerFaceColor', [0.5, 0.5, 0.5], 'MarkerEdgeColor', [0.5, 0.5, 0.5]); % Grey color
end

box on
ylim([0 1]);
set(gca, 'YTick', [0.2 0.5,0.8]);
set(gca, 'LineWidth', 1.5, 'XTickLabel', [], 'YTickLabel', []);

% Export the figure
exportgraphics(gcf, fullfile(result_dir, 'Figure4c.png'), 'Resolution', 1200);

%--------------------------------------------------------------------------%
%% Signal flow analysis
%--------------------------------------------------------------------------%
% Number of modules
nModules = NUM_MODULES_22;

% Initialize matrices for positive and negative signal flow
signal_flow_pos = zeros(nModules, nModules);
signal_flow_neg = zeros(nModules, nModules);

% Compute signal flow matrices in parallel
parfor i = 1:nModules
    [signal_flow_neg(:, i), signal_flow_pos(:, i)] = signalflow_modules(iEC, i, '22modules');
end

% Plot the signal flow using edge bundling
plot_signal_flow(nModules, signal_flow_pos, signal_flow_neg)

% Export figure
exportgraphics(gcf, fullfile(result_dir, 'Figure4d.png'), 'Resolution', 1200);
savefig(gcf,fullfile(result_dir, 'Figure4d.fig'))

%-----------------------------------%
%% Figure 3e
%-----------------------------------%
% Cross module analysis
% Define the module information with indices for Unimodal and Heteromodal categories
ModuleInfo = {'Unimodal', 1:11, 'Heteromodal', 12:22};

% Extract indices for unimodal and heteromodal modules
unimodal_idx = ModuleInfo{2};
heteromodal_idx = ModuleInfo{4};

% Reshape the positive signal flow data and filter out zero values
within_unimodal_pos = nonzeros(reshape(signal_flow_pos(unimodal_idx, unimodal_idx), [], 1));
unimodal_to_heteromodal_pos = nonzeros(reshape(signal_flow_pos(heteromodal_idx, unimodal_idx), [], 1));
heteromodal_to_unimodal_pos = nonzeros(reshape(signal_flow_pos(unimodal_idx, heteromodal_idx), [], 1));
within_heteromodal_pos = nonzeros(reshape(signal_flow_pos(heteromodal_idx, heteromodal_idx), [], 1));

% Reshape the negative signal flow data, take the absolute value, and filter out zeros
within_unimodal_neg = nonzeros(abs(reshape(signal_flow_neg(unimodal_idx, unimodal_idx), [], 1)));
unimodal_to_heteromodal_neg = nonzeros(abs(reshape(signal_flow_neg(heteromodal_idx, unimodal_idx), [], 1)));
heteromodal_to_unimodal_neg = nonzeros(abs(reshape(signal_flow_neg(unimodal_idx, heteromodal_idx), [], 1)));
within_heteromodal_neg = nonzeros(abs(reshape(signal_flow_neg(heteromodal_idx, heteromodal_idx), [], 1)));

% Combine positive and negative flow data for bar plot preparation
pos_data = [within_unimodal_pos; unimodal_to_heteromodal_pos; heteromodal_to_unimodal_pos; within_heteromodal_pos];
neg_data = [within_unimodal_neg; unimodal_to_heteromodal_neg; heteromodal_to_unimodal_neg; within_heteromodal_neg];

% Combine all data into one array
Data = [pos_data; neg_data];

% Create group labels for the data (positive and negative)
xGroup1 = repmat({'positive'}, length(pos_data), 1);  % Group for positive data
xGroup2 = repmat({'negative'}, length(neg_data), 1);  % Group for negative data
xGroup = [xGroup1; xGroup2];  % Combine group labels

% Create labels for the different categories within the groups
temp1 = repmat({'within_unimodal'}, length(within_unimodal_pos), 1);
temp2 = repmat({'unimodal_to_heteromodal'}, length(unimodal_to_heteromodal_pos), 1);
temp3 = repmat({'heteromodal_to_unimodal'}, length(heteromodal_to_unimodal_pos), 1);
temp4 = repmat({'within_heteromodal'}, length(within_heteromodal_pos), 1);
temp5 = repmat({'within_unimodal'}, length(within_unimodal_neg), 1);
temp6 = repmat({'unimodal_to_heteromodal'}, length(unimodal_to_heteromodal_neg), 1);
temp7 = repmat({'heteromodal_to_unimodal'}, length(heteromodal_to_unimodal_neg), 1);
temp8 = repmat({'within_heteromodal'}, length(within_heteromodal_neg), 1);
Legend = [temp1; temp2; temp3; temp4; temp5; temp6; temp7; temp8];  % Combine all labels

% Combine the data and labels into a table
result = [array2table(Data), cell2table([xGroup Legend])];

% Save the result table to a .mat file
save(fullfile(result_dir, 'Figure4e.mat'), 'result');

%--------------------------------------------------------------------------%
%% Summary of outputs
%--------------------------------------------------------------------------%
fprintf('Analysis completed successfully!\n');
fprintf('Generated files:\n');
fprintf('  - Figure4a.png: iEC matrix visualization\n');
fprintf('  - Figure4b.png: Edge strength distribution with insets\n');
fprintf('  - Figure4c.png: Network ratio boxplots\n');
fprintf('  - Figure4c.mat: Network ratio data table\n'); 
fprintf('  - Figure4d.png: Signal flow visualization\n');
fprintf('  - Figure4d.fig: Signal flow figure file\n');
fprintf('  - Figure4e.mat: Cross-module flow analysis data\n');
fprintf('  - Supplementary surface plots: Outflow, Inflow, Ratios\n');




